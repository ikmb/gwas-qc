// -*- mode:groovy -*-

/*
 Author: Jan Kässens <j.kaessens@ikmb.uni-kiel.de>
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

// match auto-generated "no file exists" to actual not-existing files
if (params.PCA_SNPList == "nofileexists") {
    params.PCA_SNPList = ""
}

script_dir = file(SCRIPT_DIR)
batch_dir = file(BATCH_DIR)

input_ch = Channel.create()
Channel.fromFilePairs("${params.input}.{bed,bim,fam}", size: 3, flat: true) { file -> file.baseName } \
    .ifEmpty{ error "Could not find Plink input dataset" } \
    .map { a -> [fileExists(a[1]), fileExists(a[2]), fileExists(a[3])] }
    .separate (input_ch) { a -> [a] }

sampleqci_variant_filter = file("bin/SampleQCI_variant_filter.py")
sampleqci_pca_convert = file("bin/SampleQCI_pca_convert.py")
sampleqci_pca_run = file("bin/SampleQCI_pca_run.py")

sampleqci_helpers = file("bin/SampleQCI_helpers.py")


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
    module "Plink/1.9"

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
    module "Plink/1.9"

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
    file sampleqci_helpers

    output:
    set file('pruned.bed'), file('pruned.bim'), file('pruned.fam') into for_calc_imiss,for_merge_hapmap,for_pca_convert_pruned,for_second_pca_eigen,for_second_pca_flashpca

    module "IKMB"
    module "Plink/1.9"

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
#python ${sampleqci_variant_filter} "${bim}" include_variants
python -c 'from SampleQCI_helpers import *; write_snps_autosomes_noLDRegions_noATandGC_noIndels("${bim}", "include_variants")'
# ${sampleqci_helpers} 
plink --bfile "${base}" --extract include_variants --make-bed --out pruned --allow-no-sex
"""
    }
}

process calc_imiss {
    module "IKMB"
    module "Plink/1.9"

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

//  cpus { 8 * 1 }

    module "IKMB"
    module "Plink/1.9"

    input:
    file pruned from for_merge_hapmap


    def hapmap = params.preQCIMDS_HapMap2


    output:
        set file('pruned_hapmap.bed'), file('pruned_hapmap.bim'), file('pruned_hapmap.fam') into for_pca_convert_pruned_hapmap, for_pca_run_pruned_hapmap

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
  module "Plink/1.9"
  module "Eigensoft"

  input:
  file pruned_hapmap from for_pca_convert_pruned_hapmap
  file pruned from for_pca_convert_pruned

  output:
  file 'eigenstrat-parameters'
  file 'pruned.{eigenstratgeno,ind,snp}'
  file 'pruned_hapmap.{eigenstratgeno,ind,snp}'


  def annotations = BATCH_DIR + "/" + params.individuals_annotation_hapmap2

  script:
  base_pruned = pruned[0].baseName
  base_hapmap = pruned_hapmap[0].baseName
  plink_pca = pruned[0].baseName + "_" + String.valueOf(params.numof_pc) + "PC"
"""
  python -c 'from SampleQCI_helpers import *; pca_convert ("${base_pruned}", "eigenstrat-parameters", "${annotations}", "${plink_pca}")'
  python -c 'from SampleQCI_helpers import *; pca_convert ("${base_hapmap}", "eigenstrat-parameters-all", "${annotations}", "${plink_pca}")'
#  SampleQCI_pca_convert.py "${base_pruned}" eigenstrat-parameters ${annotations}
#  SampleQCI_pca_convert.py "${base_hapmap}" eigenstrat-parameters-all ${annotations}
"""
}

projection_on_populations_controls =          script_dir + '/' + params.projection_on_populations_controls
projection_on_populations_hapmap  =     script_dir + '/' + params.projection_on_populations_hapmap 

process pca_run {
    module "IKMB"
    module "Plink/1.9"
    module "Eigensoft/4.2"

    input:
    file pruned_hapmap from for_pca_run_pruned_hapmap
    file projection_on_populations_hapmap

    script:
    base_pruned = pruned_hapmap[0].baseName
    sigma_threshold = 100.0

"""
# SampleQCI_pca_run.py "${base_pruned}" ${sigma_threshold} "${projection_on_populations_hapmap}" ${params.numof_pc}"
python -c 'from SampleQCI_helpers import *; pca_run("${base_pruned}", ${sigma_threshold}, "${projection_on_populations_hapmap}", ${params.numof_pc})'
"""
}

process second_pca_eigenstrat {
    when params.program_for_second_PCA == "EIGENSTRAT"

    input:
    file pruned from for_second_pca_eigen
    file projection_on_populations_controls

    output:
    file '*_flashpca2'

    script:
    sigma_threshold = 6.0

"""
python -c 'from SampleQCI_helpers import *; pca_run("${pruned}", ${sigma_threshold}, "${projection_on_populations_controls}", ${params.numof_pc})'
# SampleQCI_pca_run.py "${pruned}" ${sigma_threshold} "${projection_on_populations_controls}" ${params.numof_pc}
# determine PCA outlier
"""
}



process second_pca_flashpca2 {
    when params.program_for_second_PCA == "FLASHPCA2"

    module "IKMB"
    module "FlashPCA"

    input:
    file pruned from for_second_pca_flashpca

    script:
    base_pruned = pruned[0].baseName
    plink_pca = pruned[0].baseName + "_" + String.valueOf(params.numof_pc) + "PC"
    draw_evec_FLASHPCA2 = file(script_dir+"/draw_evec_FLASHPCA2.r")
    pcaplot_1KG = file(script_dir+"/pcaplot_1KG_v2.R")

"""
flashpca2 -d ${params.numof_pc} --bfile "${base_pruned}" --outval "${base_pruned}_eigenvalues_flashpca2" --outvec "${base_pruned}_eigenvectors_flashpca2" --outpc "${base_pruned}_pcs_flashpca2" --outpve "${base_pruned}_pve_flashpca2" --numthreads ${params.numof_threads} --outload "${base_pruned}_loadings_flashpca2" --outmeansd "${base_pruned}_meansd_flashpca2"

# addbatchinfo_10PCs(evec_file=${base_pruned}_pcs_flashpca2, eval_file=
R --slave --args "${plink_pca}" <"${draw_evec_FLASHPCA2}"
# addcountryinfo
R --slave --args "${plink_pca}.country" <"${draw_evec_FLASHPCA2}"

# merge__new_plink_collection_pruned__1kG
#flashpca2 -d ${params.numof_pc} 
"""

} 
