// -*- mode:groovy -*-

/*
 Author: Jan Kässens <j.kaessens@ikmb.uni-kiel.de>

 TODO:
*/

// Helper closure to check input files
fileExists = { fn ->
              if (fn.exists())
                  return fn;
              else
                  error("File not found: $fn")
}

// initialize configuration
params.output = "."
params.PCA_SNPList = ""

if (params.PCA_SNPList == "nofileexists") {
    params.PCA_SNPList = ""
}

script_dir = file(SCRIPT_DIR)

input_ch = Channel.create()
Channel.fromFilePairs("${params.input}.{bed,bim,fam}", size: 3, flat: true) { file -> file.baseName } \
    .ifEmpty{ error "Could not find Plink input dataset" } \
    .map { a -> [fileExists(a[1]), fileExists(a[2]), fileExists(a[3])] }
    .separate (input_ch) { a -> [a] }

//input_bim = file(params.input + ".bim")
//input_bed = file(params.input + ".bed")
//input_fam = file(params.input + ".fam")

sampleqci_variant_filter = file("bin/SampleQCI_variant_filter.py")
sampleqci_pca_convert = file("bin/SampleQCI_pca_convert.py")



/*
 *
 *
 * Imported from SampleQCI_parallel_part1
 *
 *
 *
 */
/* Apply pre-calculated remove list for indidividuals */
process apply_precalc_remove_list {
    input:
    file plink from input_ch
//    file precalculated_remove_list

    output:
    file "manually-removed.{bed,bim,fam}" into for_det_miss_het, for_calc_pi_hat

    module "IKMB"
    module "Plink/1.9b4.5"

//    def basename = new File(input_bim.toString()).getBaseName()
    // plink --bfile ${basename} --remove ${precalculated_remove_list} --make-bed --out manually-removed

    script:
        base = plink[0].baseName
"""
# generates manually-removed.{bim,bed,fam,log,nosex,hh}
plink --bfile ${base} --make-bed --out manually-removed --allow-no-sex
"""
}


process determine_miss_het {
    input:
    file dataset from for_det_miss_het

    output:
    file "het.het.{1.png,2.png,logscale.1.png,logscale.2.png}"
    file "miss.lmiss.{1.png,2.png,logscale.1.png,logscale.2.png}"
    file "miss.outlier.txt" into for_calc_pi_hat_outliers
    file "het.het.outlier.txt"

    publishDir params.output ?: '.', mode: 'copy', overwrite: true

    //cpus {4 * task.attempt}

    module "IKMB"
    module "Plink/1.9b4.5"

shell:
'''
# generates  miss.{hh,imiss,lmiss,log,nosex}
plink --bfile "!{new File(dataset[0].toString()).getBaseName()}" --out miss --missing --allow-no-sex&
# generates het.{het,hh,log,nosex}
plink --bfile "!{new File(dataset[0].toString()).getBaseName()}" --out het --het --allow-no-sex&
wait

R --slave --args het.het miss.imiss < "!{script_dir + "/heterozygosity_logimiss_withoutthresh.r"}"& # generates het.het.logscale.1.png
R --slave --args het.het miss.imiss < "!{script_dir + "/heterozygosity_imiss_withoutthresh.r"}"&    # generates het.het.1.png
R --slave --args miss.lmiss < "!{script_dir + "/lmiss_withoutthresh.r"}"&    # generates miss.lmiss.1.png
R --slave --args miss.lmiss < "!{script_dir + "/loglmiss_withoutthresh.r"}"& # generates miss.lmiss.logscale.1.png
wait

R --slave --args het.het miss.imiss !{params.mind} < "!{script_dir + "/heterozygosity_logimiss.r"}"& # generates het.het.logscale.2.png
R --slave --args het.het miss.imiss !{params.mind} < "!{script_dir + "/heterozygosity_imiss.r"}"&    # generates het.het.2.png
R --slave --args miss.lmiss !{params.geno_batch} < "!{script_dir + "/lmiss.r"}"&    # generates miss.lmiss.2.png
R --slave --args miss.lmiss !{params.geno_batch} < "!{script_dir + "/loglmiss.r"}"& # generates miss.lmiss.logscale.2.png
wait

R --slave --args het.het < "!{script_dir + "/heterozygosity.r"}" # generates het.het.outlier.txt

# generate individual outliers
perl -ne 'chomp;next if $.==1; @s=split /\\A\\s+|\\s+\\z/;print "$s[0]\t$s[1]\n" if $s[6]>!{params.mind};' <miss.imiss >miss.outlier.txt
'''
}

process calc_pi_hat {

    input:
    file dataset from for_calc_pi_hat
    file outliers from for_calc_pi_hat_outliers
    file sampleqci_variant_filter

    output:
    set file('pruned.bed'), file('pruned.bim'), file('pruned.fam') into for_calc_imiss,for_merge_hapmap,for_pca_convert_pruned

    module "IKMB"
    module "Plink/1.9b4.5"

    script:
        base = dataset[0].baseName
        bim = dataset[1]

    if (params.PCA_SNPList != "") {
"""
echo Using PCA SNP List file for variant selection
plink --bfile "${base}" --extract "$params.PCA_SNPList" --remove  "$outliers" --make-bed --out pruned --allow-no-sex
"""
    } else {
"""
echo Generating PCA SNP List file for variant selection
plink --bfile "${base}" --indep-pairwise 50 5 0.2 --out _prune --allow-no-sex
plink --bfile "${base}" --extract _prune.prune.in --maf 0.05 --remove "$outliers" --make-bed --out intermediate --allow-no-sex
python ${sampleqci_variant_filter} "${bim}" include_variants
plink --bfile "${base}" --extract include_variants --make-bed --out pruned --allow-no-sex
"""
    }
}

process calc_imiss {
    module "IKMB"
    module "Plink/1.9b4.5"

    input:
    file dataset from for_calc_imiss

    output:
    file "IBS.genome"
    set file("miss.imiss"), file ("miss.lmiss")

    script:
        base = dataset[0].baseName
"""
plink --bfile "${base}" --genome --out IBS --allow-no-sex&
plink --bfile "${base}" --missing --out miss --allow-no-sex&
wait
"""
}

/*
 *
 *
 * Imported from SampleQCI_parallel_part2
 *
 *
 *
 */

process merge_dataset_with_hapmap {

    cpus { 8 * 1 }

    module "IKMB"
    module "Plink/1.9b4.5"

    input:
    file pruned from for_merge_hapmap


    def hapmap = params.preQCIMDS_HapMap2


    output:
        set file('pruned_hapmap.bed'), file('pruned_hapmap.bim'), file('pruned_hapmap.fam') into for_pca_convert_pruned_hapmap

    script:
        base_pruned = pruned[0].baseName
        bim_pruned = pruned[1]

    if (params.PCA_SNPList != "") {
        """
        plink --bfile "${base_pruned}" --extract "${hapmap}.bim" --exclude "${params.PCA_SNPList}" --make-bed --out pruned_tmp --allow-no-sex
        plink --bfile "${hapmap}" --extract "${bim_pruned}" --exclude "${params.PCA_SNPList}" --make-bed --out hapmap_tmp --allow-no-sex
        plink --bfile pruned_tmp --bmerge hapmap_tmp --out pruned_hapmap --allow-no-sex
        """
    } else {
        """
        plink --bfile "${base_pruned}" --extract "${hapmap}.bim" --make-bed --out pruned_tmp --allow-no-sex
        plink --bfile "${hapmap}" --extract "${bim_pruned}" --make-bed --out hapmap_tmp --allow-no-sex
        plink --bfile pruned_tmp --bmerge hapmap_tmp --out pruned_hapmap --allow-no-sex
        """
    }
}


process pca_convert {
  module "IKMB"
  module "Plink/1.9b4.5"
  module "Eigensoft"

  input:
  file pruned_hapmap from for_pca_convert_pruned_hapmap
  file pruned from for_pca_convert_pruned

  output:
  file 'eigenstrat-parameters'
  file 'pruned.eigenstratgeno'
  file 'pruned.ind'
  file 'pruned.snp'

  def annotations = BATCH_DIR + "/" + params.individuals_annotation_hapmap2

    script:
        base_pruned = pruned[0].baseName
        base_hapmap = pruned_hapmap[0].baseName
"""
  SampleQCI_pca_convert.py "${base_pruned}" eigenstrat-parameters ${annotations}
  SampleQCI_pca_convert.py "${base_hapmap}" eigenstrat-parameters-all ${annotations}
"""
}

process pca_eigenstrat {
    when params.program_for_second_PCA == "EIGENSTRAT"


process pca_flashpca2 {
    when params.program_for_second_PCA == "FLASHPCA2"
}
