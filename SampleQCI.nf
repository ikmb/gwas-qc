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
params.sampleqci_dir = "."
params.PCA_SNPList = ""

// match auto-generated "no file exists" to actual not-existing files
if (params.PCA_SNPList == "nofileexists") {
    params.PCA_SNPList = ""
}
if (params.PCA_SNPexcludeList == "nofileexists") {
    params.PCA_SNPexcludeList = ""
}

// may be overwritten by parameters
params.individuals_remove_manually = "/dev/null"

// default value, will be overwritten by config
params.projection_on_populations_CON_only = "False"

individuals_annotation = file(params.individuals_annotation)

script_dir = file(SCRIPT_DIR)
batch_dir = file(BATCH_DIR)

input_ch = Channel.create()
Channel.fromFilePairs("${params.snpqci_dir}/${params.collection_name}_QCI.{bed,bim,fam}", size: 3, flat: true) { file -> file.baseName } \
    .ifEmpty{ error "Could not find Plink input dataset" } \
    .map { a -> [fileExists(a[1]), fileExists(a[2]), fileExists(a[3])] }
    .separate (input_ch) { a -> [a] }

for_remove_bad_samples = Channel.create()
Channel.fromFilePairs("${params.snpqci_dir}/${params.collection_name}_QCI.{bed,bim,fam}", size: 3, flat: true) { file -> file.baseName } \
    .ifEmpty{ error "Could not find Plink input dataset" } \
    .map { a -> [fileExists(a[1]), fileExists(a[2]), fileExists(a[3])] }
    .separate (for_remove_bad_samples) { a -> [a] }

sampleqci_helpers = file("${workflow.projectDir}/bin/SampleQCI_helpers.py")


prefix = params.collection_name + "_SampleQCI_"

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
    tag "${params.collection_name}"
    input:
    file plink from input_ch
//    file precalculated_remove_list

    output:
    file "manually-removed.{bed,bim,fam}" into for_det_miss_het, for_calc_pi_hat


    script:
    base = plink[0].baseName
    remove_list = "/${params.individuals_remove_manually}"

    if(params.individuals_remove_manually == '' || params.individuals_remove_manually == 'nothing.txt' || params.individuals_remove_manually == '/dev/null') {
        """
        module load "IKMB"
        module load "Plink/1.9"

        # generates manually-removed.{bim,bed,fam,log,nosex,hh}
        plink --bfile ${base} --make-bed --out manually-removed --allow-no-sex --memory ${task.memory.toMega()} 
        """
    } else {
        """
        module load "IKMB"
        module load "Plink/1.9"

        plink --bfile ${base} --remove ${remove_list} --make-bed --out manually-removed --allow-no-sex --memory ${task.memory.toMega()} 
        """
    }
}


process determine_miss_het {
//    cpus 4
    tag "${params.collection_name}"
	memory { 12.GB * task.attempt }
	time { 4.h * task.attempt }
	errorStrategy 'retry'

    input:
        file dataset from for_det_miss_het
//    val prefix

    output:
    file prefix+"het.het.{1.png,2.png,logscale.1.png,logscale.2.png}"
    file prefix+"miss.lmiss.{1.png,2.png,logscale.1.png,logscale.2.png}"
    file prefix+"miss.outlier.txt" into for_calc_pi_hat_outliers, for_remove_bad_samples_miss, for_prune_wr_miss
    file prefix+"het.het.outlier.txt" into for_remove_bad_samples_het

    publishDir params.sampleqci_dir ?: '.', mode: 'copy'

    //cpus {4 * task.attempt}



shell:
'''
    module load "IKMB"
    module load "Plink/1.7"
    
# generates  miss.{hh,imiss,lmiss,log,nosex}
plink --noweb --bfile "!{new File(dataset[0].toString()).getBaseName()}" --out !{prefix}miss --missing&
# generates het.{het,hh,log,nosex}
plink --noweb --bfile "!{new File(dataset[0].toString()).getBaseName()}" --out !{prefix}het --het&
wait

R --slave --args !{prefix}het.het !{prefix}miss.imiss < "!{script_dir + "/heterozygosity_logimiss_withoutthresh.r"}"& # generates het.het.logscale.1.png
R --slave --args !{prefix}het.het !{prefix}miss.imiss < "!{script_dir + "/heterozygosity_imiss_withoutthresh.r"}"&    # generates het.het.1.png
R --slave --args !{prefix}miss.lmiss < "!{script_dir + "/lmiss_withoutthresh.r"}"&    # generates miss.lmiss.1.png
R --slave --args !{prefix}miss.lmiss < "!{script_dir + "/loglmiss_withoutthresh.r"}"& # generates miss.lmiss.logscale.1.png
wait

R --slave --args !{prefix}het.het !{prefix}miss.imiss !{params.mind} < "!{script_dir + "/heterozygosity_logimiss.r"}"& # generates het.het.logscale.2.png
R --slave --args !{prefix}het.het !{prefix}miss.imiss !{params.mind} < "!{script_dir + "/heterozygosity_imiss.r"}"&    # generates het.het.2.png
R --slave --args !{prefix}miss.lmiss !{params.geno_batch} < "!{script_dir + "/lmiss.r"}"&    # generates miss.lmiss.2.png
R --slave --args !{prefix}miss.lmiss !{params.geno_batch} < "!{script_dir + "/loglmiss.r"}"& # generates miss.lmiss.logscale.2.png
wait

R --slave --args !{prefix}het.het < "!{script_dir + "/heterozygosity.r"}" # generates het.het.outlier.txt

# generate individual outliers
perl -ne 'chomp;next if $.==1; @s=split /\\s+/;print "$s[1]\\t$s[2]\\n" if $s[6]>!{params.mind};' <!{prefix}miss.imiss >!{prefix}miss.outlier.txt
'''
}

process prune {
    tag "${params.collection_name}"
    publishDir params.sampleqci_dir ?: '.', mode: 'copy'
    time 3.h

    input:
    file dataset from for_calc_pi_hat
    file outliers from for_calc_pi_hat_outliers

    output:
        set file(prefix+'pruned.bed'), file(prefix+'pruned.bim'), file(prefix+'pruned.fam') into for_calc_imiss,for_detect_duplicates,for_calc_imiss_ibs,for_ibs_merge_and_verify_ds,for_merge_hapmap,for_pca_convert_pruned,for_second_pca_eigen,for_second_pca_flashpca,for_second_pca_flashpca_1kg



    script:
        base = dataset[0].baseName
        bim = dataset[1]

    if (params.PCA_SNPList != "" && params.PCA_SNPList != "nofileexists") {
"""
    module load "IKMB"
    module load "Plink/1.7"
echo Using PCA SNP List file and sample outliers for variant selection
plink --noweb --bfile "${base}" --extract "$params.PCA_SNPList" --remove  "$outliers" --make-bed --out pruned
"""
    } else {
"""
    module load "IKMB"
    module load "Plink/1.9"
echo Generating PCA SNP List file for variant selection
plink --bfile "${base}" --indep-pairwise 50 5 0.2 --out _prune --memory ${task.memory.toMega()} 
plink --bfile "${base}" --extract _prune.prune.in --maf 0.05 --remove "$outliers" --make-bed --out intermediate --memory ${task.memory.toMega()} 
python -c 'from SampleQCI_helpers import *; write_snps_autosomes_noLDRegions_noATandGC_noIndels("${bim}", "include_variants")'
plink --noweb --bfile intermediate --extract include_variants --make-bed --out "${prefix}pruned"
"""
    }
}

process calc_imiss {
    tag "${params.collection_name}"
    input:
    file dataset from for_calc_imiss

    output:
    file(prefix+"miss.imiss") into for_detect_duplicates_miss
    file (prefix+"miss.lmiss")

    script:
        base = dataset[0].baseName
"""
    module load "IKMB"
    module load "Plink/1.7"
plink --noweb --bfile "${base}" --missing --out ${prefix}miss
"""
}

final calc_imiss_job_count = 5

calc_imiss_job_ids = Channel.from(1..calc_imiss_job_count) // plink expects 1-based job indices
process calc_imiss_IBS {
    input:
    tag "${params.collection_name}/$job"
    file dataset from for_calc_imiss_ibs
    each job from calc_imiss_job_ids

    output:
    file "*.part.genome.*" into for_ibs_merge_and_verify

"""
module load "IKMB"
module load "Plink/1.9"

plink --bfile ${dataset[0].baseName} --genome --parallel ${job} ${calc_imiss_job_count} --threads 1 --memory ${task.memory.toMega()} --out ${dataset[0].baseName}.part
"""
}

process ibs_merge_and_verify {

    publishDir params.sampleqci_dir ?: '.', mode: 'copy'  
    tag "${params.collection_name}/$job"
    maxRetries 5
    memory { 20.GB * (task.attempt+1)/2 }

    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    input:
    file chunks from for_ibs_merge_and_verify.collect()
    file dataset from for_ibs_merge_and_verify_ds

    output:
    file "${dataset[0].baseName}_IBS.genome" into for_detect_duplicates_genome
    file "*.png"



    shell:
//        sorted_chunks = chunks.sort({ a, b -> a.toString() <=> b.toString() })
    ibdplot = SCRIPT_DIR + "/IBD-plot-genomefile.r"
'''
    module load "IKMB"
    module load "Plink/1.9"
#IFS=' '
#UNSORTED_CHUNKS='!{chunks}'

CHUNKS=$(ls --sort=extension *.genome.*)

echo Chunks: $CHUNKS

OUTFILE="!{dataset[0].baseName}_IBS.genome"
#FIRST=1
#echo Got ${#CHUNKS} chunks.
echo "Merging..."

for chunk in $CHUNKS; do
#    echo "Processing $chunk"
#    if [ $FIRST -eq 1 ]; then
#        ((FIRST--))
        cat $chunk >>$OUTFILE
#    else
#        tail -n +2 $chunk >>$OUTFILE
#    fi
done

echo "Verifying..."

# Count samples in merged set (without header line)
LINES_MERGED=$(wc -l $OUTFILE | cut -d " " -f 1 | tr -d '[:space:]')
let LINES_MERGED--

# Count samples in original dataset (need only "half of the matrix")
LINES_MATRIX=$(wc -l !{dataset[2]} | cut -d " " -f 1 | tr -d '[:space:]')
let LINES_ORIG="$LINES_MATRIX * ($LINES_MATRIX - 1) / 2"

if [ $LINES_MERGED -ne $LINES_ORIG ]; then
    echo Warning: merged genome files are broken. >&2
    echo   FAM file !{dataset[2]} has $LINES_MATRIX samples >&2
    echo   Expected pairwise estimates: $LINES_ORIG >&2
    echo   Observed pairwise estimates: $LINES_MERGED >&2
    exit 1
else
    echo Check passed, genome files do have the expected number of pairwise estimates.
fi

gawk '{ print $7, $8 }' $OUTFILE >$OUTFILE.Z0.Z1
R --slave --args $OUTFILE.Z0.Z1 < !{ibdplot}

'''
}

process merge_dataset_with_hapmap {
    tag "${params.collection_name}"
    input:
    file pruned from for_merge_hapmap

    output:
        set file(prefix+'pruned_hapmap.bed'), file(prefix+'pruned_hapmap.bim'), file(prefix+'pruned_hapmap.fam') into for_pca_convert_pruned_hapmap, for_pca_run_pruned_hapmap

    shell:
        base_pruned = pruned[0].baseName
        bim_pruned = pruned[1]
     
     exclude = BATCH_DIR + "/" + params.PCA_SNPexcludeList
    hapmap = params.preQCIMDS_HapMap2

'''
module load IKMB
module load Plink/1.9

if [ -e "!{exclude}" ]; then
    plink --bfile "!{base_pruned}" --extract "!{hapmap}.bim" --exclude "!{exclude}" --make-bed --out pruned_tmp --allow-no-sex --memory !{task.memory.toMega()} 
    plink --bfile "!{hapmap}" --extract "!{bim_pruned}" --exclude "!{exclude}" --make-bed --out hapmap_tmp --allow-no-sex --memory !{task.memory.toMega()} 
else
    plink --bfile "!{base_pruned}" --extract "!{hapmap}.bim" --make-bed --out pruned_tmp --allow-no-sex --memory !{task.memory.toMega()}
    plink --bfile "!{hapmap}" --extract "!{bim_pruned}" --make-bed --out hapmap_tmp --allow-no-sex --memory !{task.memory.toMega()} 
fi

DONE=0
MERGE_INPUT=pruned_tmp

while [ $DONE -lt 1 ]
do
    plink --bfile $MERGE_INPUT --bmerge hapmap_tmp --out !{prefix}pruned_hapmap --allow-no-sex --memory !{task.memory.toMega()} || true
    
    if [ -e "!{prefix}pruned_hapmap.missnp" ]; then
        rm -f removelist
        while read f
        do
            CHR=$(echo $f | cut -f1 -d:)
            POS=$(echo $f | cut -f2 -d:)
            grep -E "^$CHR\\\\s.*$POS" ${MERGE_INPUT}.bim | cut -f2 -d$'\\t' >>removelist
            echo "$f" >>removelist
        done <"!{prefix}pruned_hapmap.missnp"
        plink --bfile $MERGE_INPUT --exclude removelist --make-bed --out ${MERGE_INPUT}-clean
        MERGE_INPUT=${MERGE_INPUT}-clean
        rm "!{prefix}pruned_hapmap.missnp"
    else
        DONE=1
    fi
done

'''
}


process pca_convert {

    tag "${params.collection_name}"

  input:
  file pruned_hapmap from for_pca_convert_pruned_hapmap
  file pruned from for_pca_convert_pruned

  output:
  file prefix+'eigenstrat-parameters'
  file prefix+'pruned.{eigenstratgeno,ind,snp}'
  file prefix+'pruned_hapmap.{eigenstratgeno,ind,snp}' into for_pca_run


  def annotations = params.individuals_annotation_hapmap2

  script:
  base_pruned = pruned[0].baseName
  base_hapmap = pruned_hapmap[0].baseName
"""
  module load "IKMB"
  module load "Plink/1.9"
  module load "Eigensoft"
  python -c 'from SampleQCI_helpers import *; pca_convert ("${base_pruned}", "${prefix}eigenstrat-parameters", "${annotations}")'
  python -c 'from SampleQCI_helpers import *; pca_convert ("${base_hapmap}", "${prefix}eigenstrat-parameters-all", "${annotations}")'
"""
}


process pca_run {
    publishDir params.sampleqci_dir ?: '.', mode: 'copy'
    tag "${params.collection_name}"

    input:
    file pruned_hapmap from for_pca_run_pruned_hapmap

    set file(prefix+'pruned_hapmap.eigenstratgeno'), file(prefix+'pruned_hapmap.ind'), file(prefix+'pruned_hapmap.snp') from for_pca_run

    output:
        file "*.png"
        file "*.pdf"

    script:
    base_pruned = pruned_hapmap[0].baseName
    sigma_threshold = 100.0
    draw_eigenstrat = SCRIPT_DIR + "/draw_evec_EIGENSTRAT.r"
    draw_without = SCRIPT_DIR + "/draw_evec_withoutProjection.r"
    projection_on_populations_hapmap = ANNOTATION_DIR + "/" + params.projection_on_populations_hapmap
"""
    module load "IKMB"
    module load "Plink/1.9"
    module load "Eigensoft/4.2"
export TMPDIR="\$(pwd)"
export TEMPDIR="\$(pwd)"
export TEMP="\$(pwd)"
python -c 'from SampleQCI_helpers import *; pca_run("${base_pruned}", ${sigma_threshold}, "${projection_on_populations_hapmap}", ${params.numof_pc}, 1, "${draw_eigenstrat}", "${draw_without}")'
"""
}

process second_pca_eigenstrat {
    publishDir params.sampleqci_dir ?: '.', mode: 'copy' 
    tag "${params.collection_name}"


    input:
    file pruned from for_second_pca_eigen

    output:
    file "${pruned[0].baseName}_${params.numof_pc}PC.outlier.txt" into for_remove_bad_samples_eigenstrat

    script:
    sigma_threshold = 6.0
    projection_on_populations_controls = ANNOTATION_DIR + "/${params.projection_on_populations_controls}"
"""
if [ "${params.program_for_second_PCA}" == "EIGENSTRAT" ]; then
  python -c 'from SampleQCI_helpers import *;  pca_run("${pruned}", ${sigma_threshold}, "${projection_on_populations_controls}", ${params.numof_pc})'
  python -c 'from SampleQCI_helpers import *;  determine_pca_outlier(log="${pruned[0].baseName}_${params.numof_pc}PC.log", outlier_EIGENSTRAT_file="${pruned[0].baseName}_${params.numof_pc}PC.outlier.txt")'
else
  touch ${pruned[0].baseName}_${params.numof_pc}PC.outlier.txt
fi
"""
}

process flashpca2_pruned {
     publishDir params.sampleqci_dir ?: '.', mode: 'copy'   
    tag "${params.collection_name}"

    input:
    file pruned from for_second_pca_flashpca
    file individuals_annotation

    output:
    def target = "${params.collection_name}_SampleQCI_pruned_${params.numof_pc}PC"
    set file("${target}.pca.evec"), file( "${target}.country.pca.evec") into for_extract_qced_samples_evec
    set file("${target}.eval"), file("${target}.country.eval") into for_tracy_widom_stats_eval
    file "*.png"

    script:
    base_pruned = pruned[0].baseName
    base_pruned_1kG = base_pruned + "_1kG"
    plink_pca = pruned[0].baseName + "_" + String.valueOf(params.numof_pc) + "PC"
    draw_evec_FLASHPCA2 = SCRIPT_DIR + "/draw_evec_FLASHPCA2.r"
    pcaplot_1KG = SCRIPT_DIR + "/pcaplot_1KG_v2.R"

    // individuals_annotation = ANNOTATION_DIR + "/" + params.individuals_annotation
"""
    module load "IKMB"
    module load "FlashPCA"
    module load "Plink/1.7"
  flashpca2 -d ${params.numof_pc} --bfile "${base_pruned}" --outval "${base_pruned}_eigenvalues_flashpca2" --outvec "${base_pruned}_eigenvectors_flashpca2" --outpc "${base_pruned}_pcs_flashpca2" --outpve "${base_pruned}_pve_flashpca2" --numthreads 2 --outload "${base_pruned}_loadings_flashpca2" --outmeansd "${base_pruned}_meansd_flashpca2"
  echo Adding batch info | ts
  python -c 'from SampleQCI_helpers import *; addbatchinfo_10PCs("${base_pruned}_pcs_flashpca2", "${base_pruned}_eigenvalues_flashpca2", "${plink_pca}.pca.evec", "${plink_pca}.eval", "${individuals_annotation}", "${params.preQCIMDS_1kG_sample}")'
  echo Drawing FLASHPCA2 eigenvectors for set with batch info | ts
  R --slave --args "${plink_pca}" "${params.preQCIMDS_1kG_sample}" <"${draw_evec_FLASHPCA2}"
  echo Adding country info | ts
  python -c 'from SampleQCI_helpers import *; addcountryinfo_10PCs("${base_pruned}_pcs_flashpca2", "${base_pruned}_eigenvalues_flashpca2", "${plink_pca}.country.pca.evec", "${plink_pca}.country.eval", "${individuals_annotation}", "${params.preQCIMDS_1kG_sample}")'
  echo Drawing FLASHPCA2 eigenvectors for set with country info | ts
  R --slave --args "${plink_pca}.country" "${params.preQCIMDS_1kG_sample}" <"${draw_evec_FLASHPCA2}"
"""
} 

process flashpca2_pruned_1kG {
        publishDir params.sampleqci_dir ?: '.', mode: 'copy'
    tag "${params.collection_name}"

    input:
    file pruned from for_second_pca_flashpca_1kg
    file individuals_annotation

    output:
    file "${pruned[0].baseName}_1kG_${params.numof_pc}PC.fail-pca-1KG-qc.txt" into for_remove_bad_samples_flashpca
    set file("${pruned[0].baseName}_1kG.bed"), file("${pruned[0].baseName}_1kG.bim"), file("${pruned[0].baseName}_1kG.fam") into for_pca_without_projection
    file "*.pdf"

    script:
    base_pruned = pruned[0].baseName
    base_pruned_1kG = base_pruned + "_1kG"
    plink_pca = pruned[0].baseName + "_" + String.valueOf(params.numof_pc) + "PC"
    plink_pca_1kG = base_pruned_1kG + "_" + String.valueOf(params.numof_pc) + "PC"
    draw_evec_FLASHPCA2 = SCRIPT_DIR + "/draw_evec_FLASHPCA2.r"
    pcaplot_1KG = SCRIPT_DIR + "/pcaplot_1KG_v2.R"
    if (params.PCA_SNPexcludeList == "nofileexists") {
        PCA_SNPexcludeList = ""
    } else {
        PCA_SNPexcludeList = BATCH_DIR + "/" + params.PCA_SNPexcludeList
    }

    //individuals_annotation = ANNOTATION_DIR + "/" + params.individuals_annotation
"""
    module load "IKMB"
    module load "FlashPCA"
    module load "Plink/1.9"

DONE=0
BASE_PRUNED="${base_pruned}"

while [ \$DONE -lt 1 ]
do
    echo Merge with 1kG
    python -c 'from SampleQCI_helpers import *; merge__new_plink_collection_pruned__1kG("'\$BASE_PRUNED'", "${base_pruned_1kG}", "${PCA_SNPexcludeList}", "${params.preQCIMDS_1kG}")' || true
    if [ -e "${base_pruned_1kG}-merge.missnp" ]; then
        NEW_PRUNED=\${BASE_PRUNED}_clean
        rm -f removelist
        while read f
        do
            CHR=\$(echo \$f | cut -f1 -d:)
            POS=\$(echo \$f | cut -f2 -d:)
            grep -E "^\$CHR\\\\s.*\$POS" \${BASE_PRUNED}.bim | cut -f2 -d\$'\\t' >>removelist
        done <"${base_pruned_1kG}-merge.missnp"
        plink --bfile "\$BASE_PRUNED" --exclude removelist --make-bed --out "\$NEW_PRUNED"
        BASE_PRUNED=\$NEW_PRUNED
        mv "${base_pruned_1kG}-merge.missnp" "${base_pruned_1kG}.missnp.removed"
    else
        DONE=1
    fi
done

echo PCA of 1kG-merge
flashpca2 -d ${params.numof_pc} --bfile "${base_pruned_1kG}" --outval "${base_pruned_1kG}_eigenvalues_flashpca2" --outvec "${base_pruned_1kG}_eigenvectors_flashpca2" --outpc "${base_pruned_1kG}_pcs_flashpca2" --outpve "${base_pruned_1kG}_pve_flashpca2" --numthreads 2 --outload "${base_pruned_1kG}_loadings_flashpca2" --outmeansd "${base_pruned_1kG}_meansd_flashpca2"

echo Adding batch info
python -c 'from SampleQCI_helpers import *; addbatchinfo_10PCs("${base_pruned_1kG}_pcs_flashpca2", "${base_pruned_1kG}_eigenvalues_flashpca2", "${plink_pca_1kG}.pca.evec", "${plink_pca_1kG}.eval", "${individuals_annotation}", "${params.preQCIMDS_1kG_sample}")'
echo Drawing FLASHPCA2 eigenvectors for set with batch info
echo Calling R --slave --args "${plink_pca_1kG}" "${params.preQCIMDS_1kG_sample}"  lt "${pcaplot_1KG}"
R --slave --args "${plink_pca_1kG}" "${params.preQCIMDS_1kG_sample}" <"${pcaplot_1KG}"

echo Adding country info
python -c 'from SampleQCI_helpers import *; addcountryinfo_10PCs("${base_pruned_1kG}_pcs_flashpca2", "${base_pruned_1kG}_eigenvalues_flashpca2", "${plink_pca_1kG}.country.pca.evec", "${plink_pca_1kG}.country.eval", "${individuals_annotation}", "${params.preQCIMDS_1kG_sample}")'
echo Drawing FLASHPCA2 eigenvectors for set with country info
echo Calling R --slave --args "${plink_pca_1kG}.country" "${params.preQCIMDS_1kG_sample}" lt "${pcaplot_1KG}"
R --slave --args "${plink_pca_1kG}.country" "${params.preQCIMDS_1kG_sample}" <"${pcaplot_1KG}"


"""
}

process detect_duplicates_related {
        publishDir params.sampleqci_dir ?: '.', mode: 'copy'
        memory '12 GB'
    tag "${params.collection_name}"
    input:
    file pruned from for_detect_duplicates
    file ibs from for_detect_duplicates_genome
    file miss from for_detect_duplicates_miss

    output:
    file "${params.collection_name}_SampleQCI_final_flag_relatives.txt" into for_prune_related_rel
    file "${params.collection_name}_SampleQCI_final_duplicates.txt" into for_remove_bad_samples_duplicates

    shell:
    plot_script = SCRIPT_DIR + "/IBD-plot-genomefile.r"
    threshold_duplicates = params.max_ibd_threshold_duplicates
    threshold_relatives = params.max_ibd_threshold_relatives

'''
python -c 'from SampleQCI_helpers import *; detect_duplicates_related_individuals_by_pihat("!{pruned[0].baseName}", "!{ibs}", "!{miss}", "!{threshold_duplicates}", "!{threshold_relatives}", "!{plot_script}", "!{params.collection_name}_SampleQCI_final_flag_relatives.txt", "!{params.collection_name}_SampleQCI_final_duplicates.txt", "${params.batch_mode}")'
'''
}

process remove_bad_samples {
    publishDir params.sampleqci_dir ?: '.', mode: 'copy'
    tag "${params.collection_name}"

    input:
    file pruned from for_remove_bad_samples // TODO: naming issue. These are original files from SNPQCI, not pruned files.
    file het_outliers from for_remove_bad_samples_het
    file miss_outliers from for_remove_bad_samples_miss
    file duplicates from for_remove_bad_samples_duplicates
    file eigenstrat_outliers from for_remove_bad_samples_eigenstrat
    file flashpca_outliers from for_remove_bad_samples_flashpca

    output:
        file "${params.collection_name}_SampleQCI_final.{bed,bim,fam}" into for_extract_qced_samples_ds, for_draw_histograms, for_tracy_widom_stats, for_prune_related
    file "${params.collection_name}_SampleQCI_final.log" into for_tracy_widom_stats_log
    file "remove-samples" into for_pca_without_projection_removelist



    shell:
    target_basename = params.collection_name + "_SampleQCI_final"
    remove_manually = BATCH_DIR + "/${params.individuals_remove_manually}"
'''
    module load 'IKMB'
    module load 'Plink/1.7'
touch remove-samples

# pre-calculated remove list for individuals
if [ -e !{remove_manually} ]; then
    cut -f 1,2 !{remove_manually} | tr -s ' \t' ' '>>remove-samples
fi

# heterozygosity outliers
cat !{het_outliers} | tr -s ' \t' ' '>>remove-samples

# too high missingness
cat !{miss_outliers} | tr -s ' \t' ' '>>remove-samples

# duplicates
cut -d" " -f 1,2 !{duplicates} | tr -s ' \t' ' ' >>remove-samples

# eigenstrat outliers (empty if flashpca has been used)
cat !{eigenstrat_outliers} | tr -s ' \t' ' '>>remove-samples

# flashpca outliers (empty if eigenstrat has been used)
cat !{flashpca_outliers} | tr -s ' \t' ' ' >>remove-samples

cut -f1,2 -d' ' remove-samples | sort | uniq | TMPDIR=. sponge remove-samples

plink --noweb --bfile !{pruned[0].baseName} --remove remove-samples --make-bed --out !{target_basename} --allow-no-sex
'''
}

process extract_qced_samples {
	publishDir params.sampleqci_dir ?: '.', mode: 'copy'
    tag "${params.collection_name}"
    input:
    file dataset from for_extract_qced_samples_ds
    file evec from for_extract_qced_samples_evec // 0 is .evec, 1 is .country.evec
    file individuals_annotation

    output:
    set file ("${dataset[0].baseName}.evec"), file ("${dataset[0].baseName}.country.evec") into for_draw_histograms_evec, for_prune_related_evec
    file "${dataset[0].baseName}_annotation.txt" into for_draw_histograms_ann

    shell:
    //individuals_annotation = ANNOTATION_DIR + "/${params.individuals_annotation}"
    qced_annotations = "${dataset[0].baseName}_annotation.txt"
    newevec = "${dataset[0].baseName}.evec"
    newcountryevec = "${dataset[0].baseName}.country.evec"

'''
echo "Extract QCed samples from annotation file"
python -c 'from SampleQCI_helpers import *; extract_QCsamples_from_pc_file("!{evec[0]}", "!{newevec}", "!{dataset[2]}")'
python -c 'from SampleQCI_helpers import *; extract_QCsamples_from_annotationfile("!{newevec}", "!{individuals_annotation}", "!{dataset[0].baseName}_annotation.txt")'

python -c 'from SampleQCI_helpers import *; extract_QCsamples_from_pc_file("!{evec[1]}", "!{newcountryevec}", "!{dataset[2]}")'
'''
}

process draw_histograms {
        publishDir params.sampleqci_dir ?: '.', mode: 'copy'
    tag "${params.collection_name}"
    input:
    file dataset from for_draw_histograms
    file evec from for_draw_histograms_evec // same here
    file annotations from for_draw_histograms_ann

    output:
    file "*.png"

    shell:
    histos_con_flashpca = SCRIPT_DIR + "/draw_histos_CON_CASE_FLASHPCA2.r"
    histos_all_flashpca = SCRIPT_DIR + "/draw_histos_CON_PS_AS_CD_UC_PSC_FLASHPCA2.r"
    histos_con_country_flashpca = SCRIPT_DIR + "/draw_histos_CON_CASE_FLASHPCA2_country.r"

'''
# Currently, these R scripts expect the evec files to end in ".pca.evec"
cp !{evec[0]} !{evec[0]}.pca.evec
cp !{evec[1]} !{evec[1]}.pca.evec

if [ "!{params.hf_test_CON_only}" == "True" ]; then
    R --slave --args !{evec[0]} !{annotations} < !{histos_con_flashpca}
    R --slave --args !{evec[1]} !{annotations} < !{histos_con_country_flashpca}
else
    R --slave --args !{evec[0]} !{annotations} < !{histos_con_flashpca}
    R --slave --args !{evec[1]} !{annotations} < !{histos_all_flashpca}
fi
'''
}

process tracy_widom_stats {
    publishDir params.sampleqci_dir ?: '.', mode: 'copy'
    tag "${params.collection_name}"


    input:
    file logfile from for_tracy_widom_stats_log
    file dataset from for_tracy_widom_stats
    file eval from for_tracy_widom_stats_eval

    output:
    file "*.tracy_widom"

    shell:
'''
    module load 'IKMB'
    module load 'Eigensoft/4.2'
NUM_CASES=$(grep -P 'After.*\\d+.cases' !{logfile} | cut -d' ' -f3)
NUM_CONTROLS=$(grep -P 'After .*\\d+.controls' !{logfile} | cut -d' ' -f5)
echo Cases: $NUM_CASES, controls: $NUM_CONTROLS
NUM_SAMPLES=$(($NUM_CASES + $NUM_CONTROLS))
NUM_SNPS=$(grep -P 'After.*\\d+.SNPs' !{logfile} | cut -d' ' -f8)
echo Samples: $NUM_SAMPLES, markers: $NUM_SNPS

head -n 2000 !{eval[0]} >eval.tw

if [ "!{params.projection_on_populations_CON_only}" == "True" ]; then
    if [ "$NUM_CONTROLS" -gr "$NUM_SNPS" ]; then
        let NUM_CONTROLS--
        twstats -t !{params.twtable} -i eval.tw -o !{dataset[0].baseName}.eval.tracy_widom -n "$NUM_CONTROLS"
    else
        twstats -t !{params.twtable} -i eval.tw -o !{dataset[0].baseName}.eval.tracy_widom
    fi
else
    if [ "$NUM_SAMPLES" -gr "$NUM_SNPS" ]; then
        let NUM_SAMPLES--
        twstats -t !{params.twtable} -i eval.tw -o !{dataset[0].baseName}.eval.tracy_widom -n "$NUM_SAMPLES"
    else
        twstats -t !{params.twtable} -i eval.tw -o !{dataset[0].baseName}.eval.tracy_widom
    fi
fi
'''
}

process pca_without_projection {
    publishDir params.sampleqci_dir ?: '.', mode: 'copy'
    tag "${params.collection_name}"

    input:
    file pruned_1kg from for_pca_without_projection
    file remove_list from for_pca_without_projection_removelist
    file individuals_annotation
    output:
    file "${params.collection_name}_SampleQCI_pruned_1kG_${params.numof_pc}PC.fail-pca-1KG-qc.txt"
    file "${params.collection_name}_SampleQCI_pruned_1kG_${params.numof_pc}PC.country.fail-pca-1KG-qc.txt"

    shell:
    kg_out = "${params.collection_name}_SampleQCI_pruned_1kG"
    kg_pca = "${kg_out}_${params.numof_pc}PC"
    /*individuals_annotation = ANNOTATION_DIR + "/${params.individuals_annotation}"*/
    pcaplot_1KG = SCRIPT_DIR + "/pcaplot_1KG_v2.R"
'''
    module load 'IKMB'
    module load 'Plink/1.7'
    module load 'FlashPCA'
plink --noweb \
      --bfile "!{pruned_1kg[0].baseName}" \
      --remove "!{remove_list}" \
      --make-bed \
      --out "!{kg_out}" \
      --allow-no-sex

flashpca2 -d !{params.numof_pc} \
         --bfile !{kg_out} \
         --outval !{kg_out}_eigenvalues_flashpca2 \
         --outvec !{kg_out}_eigenvectors_flashpca2 \
         --outpc !{kg_out}_pcs_flashpca2 \
         --outpve !{kg_out}_pve_flashpca2 \
         --outload !{kg_out}_loadings_flashpca2 \
         --outmeansd !{kg_out}_meansd_flashpca2 \
         --numthreads !{task.cpus}


echo Adding batch info
python -c 'from SampleQCI_helpers import *; addbatchinfo_10PCs("!{kg_out}_pcs_flashpca2", "!{kg_out}_eigenvalues_flashpca2", "!{kg_pca}.pca.evec", "!{kg_pca}.eval", "!{individuals_annotation}", "!{params.preQCIMDS_1kG_sample}")'
echo Drawing FLASHPCA2 eigenvectors for set with batch info
touch !{kg_pca}.fail-pca-1KG-qc.txt
R --slave --args "!{kg_pca}" "!{params.preQCIMDS_1kG_sample}" <"!{pcaplot_1KG}"

echo Adding country info
python -c 'from SampleQCI_helpers import *; addcountryinfo_10PCs("!{kg_out}_pcs_flashpca2", "!{kg_out}_eigenvalues_flashpca2", "!{kg_pca}.country.pca.evec", "!{kg_pca}.country.eval", "!{individuals_annotation}", "!{params.preQCIMDS_1kG_sample}")'
echo Drawing FLASHPCA2 eigenvectors for set with country info
touch !{kg_pca}.country.fail-pca-1KG-qc.txt
R --slave --args "!{kg_pca}.country" "!{params.preQCIMDS_1kG_sample}" <"!{pcaplot_1KG}"
'''
}

// SampleQCI_parallel_part3 starts here
process remove_relatives {
    publishDir params.sampleqci_dir ?: '.', mode: 'copy'
    tag "${params.collection_name}"


    input:

    file dataset from for_prune_related
    file relatives from for_prune_related_rel
    file evec from for_prune_related_evec
    file individuals_annotation
    output:

    set file("${dataset[0].baseName}_withoutRelatives.bed"), file("${dataset[0].baseName}_withoutRelatives.bim"), file("${dataset[0].baseName}_withoutRelatives.fam") into for_prune_wr
    file "${dataset[0].baseName}_withoutRelatives.pca.evec"
    file "${dataset[0].baseName}_withoutRelatives.annotation.txt"

    shell:
    /*individuals_annotation = ANNOTATION_DIR + "/${params.individuals_annotation}"*/
'''
    module load "IKMB"
    module load "Plink/1.7"
plink --noweb --bfile "!{dataset[0].baseName}" --remove !{relatives} --make-bed --out "!{dataset[0].baseName}_withoutRelatives"

python -c 'from SampleQCI_helpers import *; extract_QCsamples_from_pc_file("!{evec[0]}", "!{dataset[0].baseName}_withoutRelatives.pca.evec", "!{dataset[2]}")'
python -c 'from SampleQCI_helpers import *; extract_QCsamples_from_annotationfile("!{dataset[0].baseName}_withoutRelatives.pca.evec", "!{individuals_annotation}", "!{dataset[0].baseName}_withoutRelatives.annotation.txt")'
'''

}

process prune_withoutRelatives {

    tag "${params.collection_name}"
    input:
    file dataset from for_prune_wr
    file miss_outliers from for_prune_wr_miss

    output:
    set file("${dataset[0].baseName}_pruned.bed"), file("${dataset[0].baseName}_pruned.bim"), file("${dataset[0].baseName}_pruned.fam") into for_calc_imiss_wr, for_calc_imiss_IBS_wr, for_ibs_merge_and_verify_wr_ds



    script:
    base = dataset[0].baseName
    bim = dataset[1]
    target = "${dataset[0].baseName}_pruned"

    if (params.PCA_SNPList != "" && params.PCA_SNPList != "nofileexists") {
"""
    module load "IKMB"
    module load "Plink/1.7"
echo Using PCA SNP List file and sample outliers for variant selection
plink --noweb --bfile "${base}" --extract "$params.PCA_SNPList" --remove  "$miss_outliers" --make-bed --out ${target}
"""
    } else {
"""
    module load "IKMB"
    module load "Plink/1.9"
echo Generating PCA SNP List file for variant selection
plink --bfile "${base}" --indep-pairwise 50 5 0.2 --out _prune
plink --bfile "${base}" --extract _prune.prune.in --maf 0.05 --remove "$miss_outliers" --make-bed --out intermediate
python -c 'from SampleQCI_helpers import *; write_snps_autosomes_noLDRegions_noATandGC_noIndels("${bim}", "include_variants")'
plink --noweb --bfile intermediate --extract include_variants --make-bed --out "${target}"
"""
    }
}


calc_imiss_job_ids_wr = Channel.from(1..calc_imiss_job_count) // plink expects 1-based job indices
process calc_imiss_IBS_wr {

    tag "${params.collection_name}/$job"
    input:
    file dataset from for_calc_imiss_IBS_wr
    each job from calc_imiss_job_ids_wr

    output:
    file "*.part.genome.*" into for_ibs_merge_and_verify_wr

"""
    module load "IKMB"
    module load "Plink/1.9"
plink --bfile ${dataset[0].baseName} --genome --parallel ${job} ${calc_imiss_job_count} --threads 1 --memory ${task.memory.toMega()} --out ${dataset[0].baseName}.part
"""
}

process calc_imiss_withoutRelatives {

    tag "${params.collection_name}"
    input:
        file dataset from for_calc_imiss_wr

    output:
        file("${dataset[0].baseName}_miss.imiss")

    script:
"""
    module load "IKMB"
    module load "Plink/1.7"
plink --noweb --bfile "${dataset[0].baseName}" --missing --out ${dataset[0].baseName}_miss
"""
}

process ibs_merge_and_verify_withoutRelatives {
    publishDir params.sampleqci_dir ?: '.', mode: 'copy'
    tag "${params.collection_name}"
    memory '16 GB'
    input:
    file chunks from for_ibs_merge_and_verify_wr.collect()
    file dataset from for_ibs_merge_and_verify_wr_ds

    output:
        file "*.png"


    shell:
    ibdplot = SCRIPT_DIR + "/IBD-plot-genomefile.r"
'''

#IFS=' '
#UNSORTED_CHUNKS='!{chunks}'

CHUNKS=$(ls --sort=extension *.genome.*)

echo Chunks: $CHUNKS

OUTFILE="!{dataset[0].baseName}_IBS.genome"
#FIRST=1
#echo Got ${#CHUNKS} chunks.
echo "Merging..."

for chunk in $CHUNKS; do
#    echo "Processing $chunk"
#    if [ $FIRST -eq 1 ]; then
#        ((FIRST--))
        cat $chunk >>$OUTFILE
#    else
#        tail -n +2 $chunk >>$OUTFILE
#    fi
done

echo "Verifying..."

# Count samples in merged set (without header line)
LINES_MERGED=$(wc -l $OUTFILE | cut -d " " -f 1 | tr -d '[:space:]')
let LINES_MERGED--

# Count samples in original dataset (need only "half of the matrix")
LINES_MATRIX=$(wc -l !{dataset[2]} | cut -d " " -f 1 | tr -d '[:space:]')
let LINES_ORIG="$LINES_MATRIX * ($LINES_MATRIX - 1) / 2"

if [ $LINES_MERGED -ne $LINES_ORIG ]; then
    echo Warning: merged genome files are broken. >&2
    echo   FAM file !{dataset[2]} has $LINES_MATRIX samples >&2
    echo   Expected pairwise estimates: $LINES_ORIG >&2
    echo   Observed pairwise estimates: $LINES_MERGED >&2
    exit 1
else
    echo Check passed, genome files do have the expected number of pairwise estimates.
    echo Plotting genome file result for comparison.
fi

gawk '{ print $7, $8 }' $OUTFILE >$OUTFILE.Z0.Z1
R --slave --args $OUTFILE.Z0.Z1 < !{ibdplot}
'''
}


