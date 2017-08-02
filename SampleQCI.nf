// -*- mode:groovy -*-

/*
 Author: Jan Kässens <j.kaessens@ikmb.uni-kiel.de>

 TODO:
*/



// initialize configuration
params.output = "."
input_basename = params.input

script_dir = file(SCRIPT_DIR)

input_bim = file(params.input + ".bim")
input_bed = file(params.input + ".bed")
input_fam = file(params.input + ".fam")


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
    file input_bim
    file input_bed
    file input_fam
//    file precalculated_remove_list

    output:
    file "manually-removed.{bim,bed,fam}" into for_det_miss_het, for_calc_pi_hat

    module "IKMB"
    module "Plink/1.9b4.5"

    def basename = new File(input_bim.toString()).getBaseName()
// plink --bfile ${basename} --remove ${precalculated_remove_list} --make-bed --out manually-removed

"""
# generates manually-removed.{bim,bed,fam,log,nosex,hh}
plink --bfile ${basename} --make-bed --out manually-removed
"""
}


process determine_miss_het {
    cpus 2

    input:
    file dataset from for_det_miss_het

    output:
    file "het.het.{1.png,2.png,logscale.1.png,logscale.2.png}"
    file "miss.lmiss.{1.png,2.png,logscale.1.png,logscale.2.png}"
    file "{het.het,miss}.outlier.txt" into for_calc_pi_hat_outliers

    publishDir params.output ?: '.', mode: 'copy', overwrite: true

    //cpus {4 * task.attempt}

    module "IKMB"
    module "Plink/1.9b4.5"

shell:
'''
# generates  miss.{hh,imiss,lmiss,log,nosex}
plink --bfile "!{new File(dataset[0].toString()).getBaseName()}" --out miss --missing&
# generates het.{het,hh,log,nosex}
plink --bfile "!{new File(dataset[0].toString()).getBaseName()}" --out het --het&
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
perl -ne 'chomp;next if $.==1; @s=split /\\s+/;print "$s[6]\n" if $s[6]>!{params.mind};' <miss.imiss >miss.outlier.txt
'''
}

process calc_pi_hat {

    input:
    file dataset from for_calc_pi_hat
    file outliers from for_calc_pi_hat_outliers

    output:
    file "pruned.{bim,bed,fam}" into for_calc_imiss,for_merge_hapmap,for_pca_convert_pruned

    module "IKMB"
    module "Plink/1.9b4.5"

    script:
    if (binding.variables.containsKey("params.PCA_SNPList") or params.PCA_SNPList != "nofileexists") {
"""
echo Using PCA SNP List file for variant selection
plink --bfile "${new File(dataset[0].toString()).getBaseName()}" --extract "$params.PCA_SNPList" --remove  "$outliers" --out pruned
"""
    } else {
"""
echo Generating PCA SNP List file for variant selection
plink --bfile "${new File(dataset[0].toString()).getBaseName()}" --indep-pairwise 50 5 0.2 --out _prune
plink --bfile "${new File(dataset[0].toString()).getBaseName()}" --extract _prune.prune.in --maf 0.05 --remove "$outliers" --make-bed --out intermediate
python SampleQCI_variant_filter.py "${dataset[0].toString()}" include_variants
plink --bfile "${new File(dataset[0].toString()).getBaseName()}" --extract include_variants --make-bed --out pruned
"""
    }
}

process calc_imiss {
    module "IKMB"
    module "Plink/1.9b4.5"

    input:
    file dataset from for_calc_imiss

    output:
    file "IBS.{bim,bed,fam}"
    file "miss.{bim,bed,fam}"

"""
plink --bfile "${new File(dataset[0].toString()).getBaseName()}" --genome --out IBS&
plink --bfile "${new File(dataset[0].toString()).getBaseName()}" --missing --out miss&
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
    module "IKMB"
    module "Plink/1.9b4.5"

    input:
    file pruned from for_merge_hapmap
    file hapmap

    output:
    file "pruned_hapmap.{bim,bed,fam}" into for_pca_convert_pruned_hapmap

    script:
    if (binding.variables.containsKey("params.PCA_SNPList") or params.PCA_SNPList != "nofileexists") {
        """
        plink --bfile "${new File(pruned[0].toString()).getBaseName()}" --extract "${hapmap[0]}" --exclude "${params.PCA_SNPList}" --make-bed --out pruned_tmp
        plink --bfile "${new File(hapmap[0].toString()).getBaseName()}" --extract "${pruned[0]}" --exclude "${params.PCA_SNPList}" --make-bed --out hapmap_tmp
        plink --bfile pruned_tmp --bmerge hapmap_tmp --out pruned_hapmap
        """
    } else {
        """
        plink --bfile "${new File(pruned[0].toString()).getBaseName()}" --extract "${hapmap[0]}" --make-bed --out pruned_tmp
        plink --bfile "${new File(hapmap[0].toString()).getBaseName()}" --extract "${pruned[0]}" --make-bed --out hapmap_tmp
        plink --bfile pruned_tmp --bmerge hapmap_tmp --out pruned_hapmap
        """
    }
}

process pca_convert {
    module "IKMB"
    module "Plink/1.9b4.5"

    input:
    file pruned_hapmap from for_pca_convert_pruned_hapmap
    file pruned from for_pca_convert_pruned

    
}
