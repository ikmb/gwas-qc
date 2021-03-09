// -*- mode:groovy -*-
// vim: syntax=nextflow

/*
 Author: Jan Kässens <j.kaessens@ikmb.uni-kiel.de>
*/


// initialize configuration
params.snpqci_dir = "."
params.skip_snpqc = 0

// hwe_template_script = file("${workflow.projectDir}/"+params.hwe_template)
hwe_script = file("${workflow.projectDir}/bin/hwe.R")

// Lots of indirection layers require lots of backslash escaping
individuals_annotation = file(params.individuals_annotation)
definetti_r = file(SCRIPT_DIR + "/DeFinetti_hardy.r")
//autosomes = file(ANNOTATION_DIR + "/" + params.chip_build + "/" + ChipDefinitions.Producer(params.chip_producer) + "/" + ChipDefinitions.RsAutosomes(params.chip_version))
draw_fdr = file("${workflow.projectDir}/bin/SNP_QCI_draw_FDR.r")
draw_fdr_allbatches = file("${workflow.projectDir}/bin/SNP_QCI_draw_FDR_Fail_allbatches.r")

// set up channels
to_calc_hwe_script = Channel.create()
to_calc_hwe = Channel.create()

process merge_batches {
    label 'big_mem'
    tag "${params.collection_name}"
    file individuals_annotation

    output:
    file "${params.collection_name}_Rs.bim" into merged_bim, to_split_bim, to_hwe_bim, to_verify_bim, to_exclude_bim, to_miss_bim, to_miss_batch_bim
    file "${params.collection_name}_Rs.bed" into merged_bed, to_hwe_bed, to_verify_bed, to_exclude_bed, to_miss_bed, to_miss_batch_bed
    file "${params.collection_name}_Rs.fam" into merged_fam, to_hwe_fam, to_verify_fam, to_exclude_fam, to_miss_fam, to_miss_batch_fam
    file "${params.collection_name}_Rs.log"
    file "${params.collection_name}.indels" into preqc_hwe_indels,postqc_hwe_indels
    file "ethnicities.txt" into ethnicities

    shell:

'''
#!/usr/bin/env bash

module load IKMB
module load Plink/1.9

MEM=!{task.memory.toMega()-1000}

set -x

# Param input: space-separated file names, relative to params.rs_dir
BIMS=(!{params.rs_bims})
BEDS=(!{params.rs_beds})
FAMS=(!{params.rs_fams})
INDELS=(!{params.rs_indels})

IFS=' '

rm -f merge-list
for idx in ${!BIMS[@]}; do
    # All values except the first go into the merge list
    if [ "$idx" -gt 0 ]; then
        echo -n !{params.rs_dir}/${BEDS[$idx]} >>merge-list
        echo -n " " >>merge-list
        echo -n !{params.rs_dir}/${BIMS[$idx]} >>merge-list
        echo -n " " >>merge-list
        echo !{params.rs_dir}/${FAMS[$idx]} >>merge-list
    fi
done

if [ "${#BIMS[*]}" -gt 1 ]; then
    plink --memory $MEM --bed !{params.rs_dir}/${BEDS[0]} --bim !{params.rs_dir}/${BIMS[0]} --fam !{params.rs_dir}/${FAMS[0]} --merge-list merge-list --make-bed --out !{params.collection_name}_Rs --allow-no-sex
    if [ -s "!{params.collection_name}_Rs-merge.missnp" ]; then
        echo "Removing $(wc -l !{params.collection_name}_Rs.missnp) variants that would have 3 or more alleles after merge."

        # Remove from every entry in the merge list
        rm -f merge-list
        for idx in ${!BIMS[@]}; do
            plink --memory $MEM --bed !{params.rs_dir}/${BEDS[$idx]} --bim !{params.rs_dir}/${BIMS[$idx]} --fam !{params.rs_dir}/${FAMS[$idx]} --exclude !{params.collection_name}_Rs-merge.missnp --allow-no-sex --make-bed --out ds-$idx

            # All values except the first go into the merge list
            if [ "$idx" -gt 0 ]; then
                echo ds-$idx.bed ds-$idx.bim ds-$idx.fam >>merge-list
            fi
        done

        plink --memory $MEM --bfile ds-0 --merge-list merge-list --make-bed --out !{params.collection_name}_Rs --allow-no-sex
    fi
else
    echo "No merge necessary, only one batch found." | tee !{params.collection_name}_Rs.log
    ln -s !{params.rs_dir}/${BIMS[0]} !{params.collection_name}_Rs.bim
    ln -s !{params.rs_dir}/${BEDS[0]} !{params.collection_name}_Rs.bed
    ln -s !{params.rs_dir}/${FAMS[0]} !{params.collection_name}_Rs.fam
fi
ln -s !{params.rs_dir}/${INDELS[0]} !{params.collection_name}.indels
plinkinfo.pl !{params.collection_name}_Rs.bim !{params.collection_name}_Rs.fam >info.txt
gawk 'NR>1 { ethn[$8] += 1 } END { for(key in ethn) { print ethn[key] " " key } }' !{individuals_annotation} | sort --reverse --numeric-sort >ethnicities.txt

'''
}

/*
 Generate HWE tables and draw DeFinetti plots of the whole data set
 */
process hwe_definetti_preqc {
  tag "${params.collection_name}"
  publishDir params.snpqci_dir ?: '.', mode: 'copy', overwrite: true

  input:
    file definetti_r
    file merged_bim
    file merged_bed
    file merged_fam
    file indels from preqc_hwe_indels
  output:
    file "${params.collection_name}_controls_DeFinetti.jpg"
    file "${params.collection_name}_cases_DeFinetti.jpg"
    file "${params.collection_name}_cases_controls_DeFinetti.jpg"
    file "${params.collection_name}_hardy.hwe"
    file "${params.collection_name}_hardy.log"

"""
    module load IKMB
    module load Plink/1.9
MEM=${task.memory.toMega()-1000}

#plink --memory \$MEM --bfile ${merged_bim.baseName} --exclude ${indels} --make-bed --out no-indels
#plink --memory \$MEM --bfile no-indels --hardy --out ${params.collection_name}_hardy --hwe 0.0 --chr 1-22 --allow-no-sex
plink --memory \$MEM --bfile ${merged_bim.baseName} --exclude ${indels} --hardy --out ${params.collection_name}_hardy --hwe 0.0 --chr 1-22 --allow-no-sex
R --slave --args ${params.collection_name}_hardy.hwe ${params.collection_name}_controls_DeFinetti ${params.collection_name}_cases_DeFinetti ${params.collection_name}_cases_controls_DeFinetti <$definetti_r
"""
}

/*
 Split the data set into 1000-SNP chunks so they can be evaluated in parallel.
 */

process split_dataset {
  label 'small_mem'
  label 'short_running'

  tag "${params.collection_name}"
input:
  file input_bim from to_split_bim
  file individuals_annotation

output:
  file 'chunk_*' into to_calc_hwe

shell:
'''
  CONTROLS=$(awk '{if($9=="Control") print}' "!{individuals_annotation}"| wc -l)
  IS_QUANT=0
  cut -f6 "!{individuals_annotation}" | tail -n +2 | sort | uniq >phenotypes
  while read pheno; do
    case "$pheno" in
        -9|0|1|2) 
            ;;
        *) IS_QUANT=1
            ;;
    esac
  done <phenotypes

  echo $IS_QUANT >is_quantitative.txt

  if [ "$CONTROLS" -ge 1 ] || [ "$IS_QUANT" == "1" ]; then
    cut -f 2 !{input_bim} | split -l 1000 -a 5 -d - chunk_
  else
    touch chunk_00000
  fi
'''
}


/*
 Run the HWE calculation R script on every 1000-SNP chunk.

 HWE Rules:
 - Run HWE only on Controls with chr 1-22+25 (XY PAR regions), add 23 only for females
 */

process calculate_hwe {
  label 'long_running'
  label 'small_mem'
  // Should have been zero, but killall -q returns 1 if it didn't find anything
  errorStrategy 'retry'
  tag "${params.collection_name}/${chunk}"

  input:
    each file(chunk) from to_calc_hwe.flatten()
    file input_bim from to_hwe_bim
    file input_bed from to_hwe_bed
    file input_fam from to_hwe_fam
    file individuals_annotation
    file ethnicities

  output:
  file "${chunk}-out.auto.R" into from_calc_hwe

shell:
'''
#!/usr/bin/env bash

    module load IKMB
    module load Plink/1.9

MEM=!{task.memory.toMega()-1000}

# Filter annotations to include only controls
<!{individuals_annotation} mawk ' { if($9 == "Control" || $9 == "diagnosis") print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10 }' >individuals_annotation
NUM_CONTROLS=$(wc -l <individuals_annotation)

IS_QUANT=0
cut -f6 "!{individuals_annotation}" | tail -n +2 | sort | uniq >phenotypes
while read pheno; do
    case "$pheno" in
        -9|0|1|2) 
            ;;
        *) IS_QUANT=1
            ;;
    esac
done <phenotypes

# note that NUM_CONTROLS includes the header line so "1" means 0 controls
# Overwrite the local copy if we're checking a quantitative trait but we have
# no defined controls:
if [ "$NUM_CONTROLS" == "1" ] && [ "$IS_QUANT" == "1" ]; then
    cp !{individuals_annotation} individuals_annotation
    # if it's quantitative and no controls have been defined, use everything.
    FILTER_CONTROLS=""
elif [ "$NUM_CONTROLS" -gt "1" ] && [ "$IS_QUANT" == "1" ]; then
    # Plink can't handle control phenotype in a quantitative setting, so just
    # handle them explicitly
    cut -f1,2 | tail -n +2 individuals_annotation >keep.fam
    FILTER_CONTROLS="--keep keep.fam"
else
    FILTER_CONTROLS="--filter-controls"
fi

# sorted by count, take the ethnicity with the highest count
ETHNICITY=$(head -n1 !{ethnicities} | cut -f2 -d" ")

# annotation file always has at least one header line
if [ "$NUM_CONTROLS" -gt 1 ] || [ "$IS_QUANT" == "1" ]; then

    # Extract our chunk from the main dataset, keep only the controls.
    plink --memory $MEM --bfile "!{new File(input_bim.toString()).getBaseName()}" $FILTER_CONTROLS --extract !{chunk} --make-bed --out !{chunk}-controls.all

    # It might be, that this chunk lies entirely in chr23, so the next one will fail and
    # we don't need to merge later on.
    plink --memory $MEM --bfile !{chunk}-controls.all --chr 1-22,25 --make-bed --out 1-22and25 || true
    HAS_OTHERS=1
    if [ ! -e 1-22and25.bim ]; then
        HAS_OTHERS=0
    fi

    # If we do have X chromosomes, extract them for females only and merge into main dataset
    HAS_X=$(cut -f1 !{chunk}-controls.all.bim | grep -c 23)
    if [ $HAS_X -ne 0 ]; then
        if [ $HAS_OTHERS -eq 1 ]; then
            plink --memory $MEM --allow-no-sex --bfile !{chunk}-controls.all --filter-females --chr 23 --make-bed --out 23
            plink --memory $MEM --allow-no-sex --bfile 1-22and25 --bmerge 23 --make-bed --out !{chunk}-controls
        else
            plink --memory $MEM --allow-no-sex --bfile !{chunk}-controls.all --filter-females --chr 23 --make-bed --out !{chunk}-controls
        fi

    else
        ln -s !{chunk}-controls.all.bim !{chunk}-controls.bim
        ln -s !{chunk}-controls.all.bed !{chunk}-controls.bed
        ln -s !{chunk}-controls.all.fam !{chunk}-controls.fam
    fi

    # Calc HWE
    Rscript !{hwe_script} !{chunk}-controls individuals_annotation !{chunk}-out.auto.R $ETHNICITY
else
    touch !{chunk}-out.auto.R
fi

'''
}


/*
 Merge all chunked HWE tables into a single file and screen for obvious errors (i.e. N/A HWE values or wrong SNP counts)
 */

process hwe_fdr_filter {
    label 'long_running'
  tag "${params.collection_name}"
publishDir params.snpqci_dir ?: '.', mode: 'copy'
    input:
    file 'chunk' from from_calc_hwe.collect()
    file input_bim from to_verify_bim
    file input_bed from to_verify_bed
    file input_fam from to_verify_fam

    file individuals_annotation
    file draw_fdr

    output:
    file ("${params.collection_name}_exclude-whole-collection-worst-batch-removed") into excludes_whole
    file ("${params.collection_name}_exclude-per-batch-fail-in-2-plus-batches") into excludes_perbatch
    file ("${params.collection_name}_exclude-whole-collection-all-batches") into excludes_allbatches
    file ("${params.collection_name}_exclude-whole-collection-worst-batch-removed.FDRthresholds.SNPQCI.1.txt.png")
    file ("${params.collection_name}_exclude-per-batch-fail-in-2-plus-batches.FDRthresholds.SNPQCI.2.txt.png")
    file ("${params.collection_name}_exclude-whole-collection-worst-batch-removed.FDRthresholds.SNPQCI.1.txt")
    file ("${params.collection_name}_exclude-per-batch-fail-in-2-plus-batches.FDRthresholds.SNPQCI.2.txt")

    shell:
    hwe_result = file("${params.collection_name}_chunks_combined.hwe")
'''
#!/usr/bin/env bash

HWE_RESULTS="!{params.collection_name}_chunks_combined.hwe"
cat chunk* >$HWE_RESULTS


N_CONTROLS=$(awk '{if($9=="Control") print}' "!{individuals_annotation}"| wc -l)

# see if we really do not have controls or if some upstream process just failed silently
# ...then early abort. Create dummy outputs for the non-optional files
if [ "$N_CONTROLS" -eq 0 ] && [ ! -s "$HWE_RESULTS" ]; then
    touch "!{params.collection_name}_exclude-whole-collection-worst-batch-removed"
    touch "!{params.collection_name}_exclude-per-batch-fail-in-2-plus-batches"
    touch "!{params.collection_name}_exclude-whole-collection-all-batches"
    touch "!{params.collection_name}_exclude-whole-collection-worst-batch-removed.FDRthresholds.SNPQCI.1.txt.png"
    touch "!{params.collection_name}_exclude-per-batch-fail-in-2-plus-batches.FDRthresholds.SNPQCI.2.txt.png"
    touch "!{params.collection_name}_exclude-whole-collection-worst-batch-removed.FDRthresholds.SNPQCI.1.txt"
    touch "!{params.collection_name}_exclude-per-batch-fail-in-2-plus-batches.FDRthresholds.SNPQCI.2.txt"
    exit 0
fi

N_CONTROLS=$(<!{individuals_annotation} tr -s '\\t ' ' ' | cut -f7,9 -d' ' | uniq | tail -n +2 | cut -f2 -d' ' | grep -c Control)

echo Found $N_CONTROLS control batches

if [ "$N_CONTROLS" -gt 2 ]; then
    DRAWSCRIPT="!{draw_fdr}"
else
    DRAWSCRIPT="!{draw_fdr_allbatches}"
fi

fdrfilter.pl "$HWE_RESULTS" !{params.FDR_index_remove_variants} \
    !{params.collection_name}_exclude-whole-collection-worst-batch-removed.FDRthresholds.SNPQCI.1.txt \
    !{params.collection_name}_exclude-whole-collection-all-batches \
    !{params.collection_name}_exclude-whole-collection-worst-batch-removed \
    !{params.collection_name}_exclude-per-batch-fail-in-2-plus-batches.FDRthresholds.SNPQCI.2.txt \
    !{params.collection_name}_exclude-per-batch-fail-in-2-plus-batches

R --slave --args \
    !{params.collection_name}_exclude-whole-collection-worst-batch-removed.FDRthresholds.SNPQCI.1.txt \
    !{params.collection_name}_exclude-per-batch-fail-in-2-plus-batches.FDRthresholds.SNPQCI.2.txt \
    !{params.FDR_index_remove_variants} \
    <$DRAWSCRIPT
'''
}


/*
 Determine the missingness for the entire collection
 */
process determine_missingness_entire {
    publishDir params.snpqci_dir ?: '.', mode: 'copy'
    errorStrategy 'retry'
    tag "${params.collection_name}"
    label 'big_mem'

    input:
    file input_bim from to_miss_bim
    file input_bed from to_miss_bed
    file input_fam from to_miss_fam

    output:
    file 'missingness-excludes-entire' into excludes_miss_entire


shell:
'''
module load IKMB
module load Plink/1.9
MEM=!{task.memory.toMega()-1000}
plink --memory $MEM --bfile "!{new File(input_bim.toString()).getBaseName()}" --missing --out missingness_entire --allow-no-sex
SNPQCI_extract_missingness_entire.py missingness_entire.lmiss !{params.geno_entire_collection} missingness-excludes-entire
'''
}


process determine_missingness_per_batch {
    publishDir params.snpqci_dir ?: '.', mode: 'copy'
    errorStrategy 'retry'
    tag "${params.collection_name}"
    label 'big_mem'

    input:
    file individuals_annotation
    file input_bim from to_miss_batch_bim
    file input_bed from to_miss_batch_bed
    file input_fam from to_miss_batch_fam

    output:
    file 'missingness-excludes-perbatch' into excludes_miss_perbatch


shell:
'''
    module load IKMB
    module load Plink/1.9

MEM=!{task.memory.toMega()-1000}

mawk '{print $1, $2, $7 }' "!{individuals_annotation}" | grep -v "familyID" >cluster_file
plink --memory $MEM --bfile "!{new File(input_bim.toString()).getBaseName()}" --missing \
    --out missingness_perbatch \
    --allow-no-sex \
    --within cluster_file
SNPQCI_extract_missingness_perbatch.py missingness_perbatch.lmiss !{params.geno_batch} "!{individuals_annotation}" missingness-excludes-perbatch
'''
}

process exclude_bad_variants {
    publishDir params.snpqci_dir ?: '.', mode: 'copy'
    errorStrategy 'retry'
    tag "${params.collection_name}"

    input:
    file individuals_annotation
    file input_bim from to_exclude_bim
    file input_bed from to_exclude_bed
    file input_fam from to_exclude_fam
    file excludes_whole from excludes_whole
    file excludes_perbatch from excludes_perbatch
    file excludes_allbatches from excludes_allbatches
    file missingness_excludes_entire from excludes_miss_entire
    file missingness_excludes_perbatch from excludes_miss_perbatch

    output:
    file "${params.collection_name}_QCI.{bim,bed,fam}" into draw_definetti_after
    file "${params.collection_name}_QCI.log"
    file "variant-excludes"


"""
    module load IKMB
    module load Plink/1.9

MEM=${task.memory.toMega()-1000}

TMPDIR=.
NUM_CTRL_BATCHES=\$(tr -s '\\t ' ' ' <${individuals_annotation} | cut -f7,9 -d" " | grep Control | uniq | wc -l)

if [ "\$NUM_CTRL_BATCHES" -gt "4" ]; then
    echo "Found \$NUM_CTRL_BATCHES control batches."
    (tail -n +1 "$excludes_whole" | cut -f1; cat "$excludes_perbatch"; cat "$missingness_excludes_entire"; cat "$missingness_excludes_perbatch") | sort -n | uniq >variant-excludes
else
    echo "Found \$NUM_CTRL_BATCHES control batches. Skipping 'worst batch removed' excludes from HWE testing."
    (tail -n +1 "$excludes_allbatches" | cut -f1; cat "$missingness_excludes_entire"; cat "$missingness_excludes_perbatch") | sort -n | uniq >variant-excludes
#    (cat "$missingness_excludes_entire"; cat "$missingness_excludes_perbatch") | sort -n | uniq >variant-excludes
fi

if [ ${params.skip_snpqc} -eq 0 ]; then
    plink --memory \$MEM --bfile "${new File(input_bim.toString()).getBaseName()}" --exclude variant-excludes --make-bed --out ${params.collection_name}_QCI
else
    BASE="${new File(input_bim.toString()).getBaseName()}"
    NEWBASE="${params.collection_name}_QCI"
    ln -s "\$BASE".bim "\$NEWBASE".bim
    ln -s "\$BASE".bed "\$NEWBASE".bed
    ln -s "\$BASE".fam "\$NEWBASE".fam
    touch "\$NEWBASE".log
fi

plinkinfo.pl "${params.collection_name}_QCI.bim" "${params.collection_name}_QCI.fam" >info.txt

"""
}

process hwe_definetti_qci {
    publishDir params.snpqci_dir ?: '.', mode: 'copy'
    tag "${params.collection_name}"

    def prefix = "${params.collection_name}_QCI_hardy"

    input:
    file new_plink from draw_definetti_after
    file definetti_r
    file indels from postqc_hwe_indels

    output:
    file prefix+".hwe"
    file prefix+"_{cases,controls,cases_controls}_DeFinetti.jpg"
    file prefix+".log"

"""
    module load IKMB
    module load Plink/1.9

MEM=${task.memory.toMega()-1000}
#plink --bfile "${new File(new_plink[0].toString()).getBaseName()}" --exclude ${indels} --make-bed --out no-indels --memory \$MEM
#plink --bfile no-indels --hardy --out ${prefix} --hwe 0.0 --chr 1-22 --allow-no-sex --memory \$MEM


plink --bfile "${new File(new_plink[0].toString()).getBaseName()}" --exclude ${indels} --hardy --hwe 0.0 --chr 1-22 --allow-no-sex --out ${prefix} --memory \$MEM

R --slave --args ${prefix}.hwe ${prefix}_controls_DeFinetti ${prefix}_cases_DeFinetti ${prefix}_cases_controls_DeFinetti <"$definetti_r"
"""
}

