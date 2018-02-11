// Helper closure to check input files
fileExists = { fn ->
              if (fn.exists())
                  return fn;
              else
                  error("File not found: $fn")
}

getCanonicalFileType = { name ->
    def result
    def filename = name.toString()
    switch(filename) {
         // extension is (traditionally) txt, so simple extraction is not possible
        case ~/.*_flag.relatives.txt/:
            result = "relatives"
            break
        case ~/.*annotation.txt/:
            result = "annotation"
            break
        case ~/.*country.eval/:
            result = "country_eval"
            break
        case ~/.*nearly_monomorphic.txt/:
            result = "nearly_monomorphic"
            break
        case ~/.*monomorphic.txt/:
            result = "monomorphic"
            break
        case ~/.*country.pca.evec.*/:
            result = "country_evec"
            break
        case ~/.*PLINKdosage.map.*/:
            result = "PLINKdosage_map"
            break
        case ~/.*PLINKdosage.gz.*/:
            result = "PLINKdosage_gz"
            break
        case ~/.*dat.pca.evec.*/:
            result = "dat_pca_evec"
            break
        // catches bim,bed,fam,log,eval,evec and all the other "simple" types
        default:
            result = filename.drop(filename.lastIndexOf('.')+1) // if the last '.' is on index 10, drop the first 11 characters, leaving the rest
            break
    }
//    println "Filing $filename in $result"
    result
}

// Transform a list of file names (strings) into a map with the file extension as key and the filename as value
mapFileList = { fn -> fn.collectEntries {
        [(getCanonicalFileType(it)): it]
    }
}

ds_stats_input = ["${params.input_stats}.fam", 
                  "${params.input_stats}_annotation.txt"
                  ].collect { fileExists(file(it)) }
ds_input = ["${params.input}.bed", 
            "${params.input}.bim", 
            "${params.input}.fam", 
            "${params.input}.dat.pca.evec",
            "${params.input}_flag.relatives.txt"].collect { fileExists(file(it)) }
ds_imp_input = ["${params.input_imp}.map", 
                "${params.input_imp}.gz", 
                "${params.input_imp}.PLINKdosage.map",
                "${params.input_imp}.PLINKdosage.gz"
                ].collect { fileExists(file(it)) }

process extract_genotyped_variants {
	cpus 4
	
	input:
	file ds_stats_staged from Channel.from(ds_stats_input).collect()
	file ds_staged from Channel.from(ds_input).collect()
	
	output:
	file "${target}.fam" into for_plink_dosage_logistic
	file "${target}.{bim,bed,fam,log}" into for_plink_dosage_logistic_ds, for_qq_manhattan_ds
	
    shell:
    ds_stats = mapFileList(ds_stats_staged)
    ds = mapFileList(ds_staged)
    target = "${params.disease_data_set_prefix_release_statistics}"
    
    // TODO FIXME --keep-fam räumt alles ab "No people remaining after --keep-fam"
'''
module load Plink/1.9

gawk '{ print $1, $2, $6 }' "!{ds.fam}" >"!{ds.fam}.single-id.pheno"
plink --bfile "!{ds.bim.baseName}" \
	  --threads 4 \
	  --remove "!{ds.relatives}" \
 	  --pheno "!{ds.fam}.single-id.pheno" \
	  --make-bed --allow-no-sex \
	  --out "!{target}"

#	  --keep-fam "!{ds_stats.fam}" \
gawk '{ print $1\"_\"$2, $1\"_\"$2, $6 }' "!{target}.fam" >"!{target}.fam.double-id.pheno"
gawk '{ print $1\"_\"$2, $1\"_\"$2, $5 }' "!{target}.fam" >"!{target}.fam.double-id.gender"
gawk '{ print $1\"_\"$2, $1\"_\"$2, $3, $4, $5, $6 }' "!{target}.fam" | sponge "!{target}.fam"
'''
}

// will be overwritten by actual parameter file
params.first_chr = 1
params.last_chr = 22

// Steps 1.1 and 1.2: per-chromosome dosage and logistic calculations,
// results will be merged in merge_log_dos_results
process plink_dosage_logistic {
    cpus 4
    tag "chr$chromosome"
    
    input:
    each chromosome from Channel.from(params.first_chr .. params.last_chr)
    file (fam:"new_fam") from for_plink_dosage_logistic
	file ds_stats_staged from for_plink_dosage_logistic_ds
	file ds_imp_staged from Channel.from(ds_imp_input).collect()
	file ds_release_staged from Channel.from(ds_input).collect()
    
    output:
//    file "${target}_chr${chromosome}.assoc.dosage" into merge_log_dos_results_dosage
    file "${dosage_target}" into merge_log_dos_results_dosage
    file "${log_prefix}.assoc.logistic" into merge_log_dos_results_logistic
    
    shell:
    ds_stats = mapFileList(ds_stats_staged)
    ds_imp = mapFileList(ds_imp_staged)
    ds_release = mapFileList(ds_release_staged)
    covar_name = "PC${params.numofPCs}"
    if(params.numofPCs > 1) {
        covar_name = "PC1-PC${params.numofPCs}"
    }
    target = "${params.disease_data_set_prefix_release_statistics}.imputed"
    dosage_src = "${target}_chr${chromosome}.assoc.dosage"
    dosage_target = "${target}_chr${params.first_chr}-${params.last_chr}.assoc.dosage.${chromosome}"
    log_prefix = "${params.disease_data_set_prefix_release_statistics}.genotyped.chr${chromosome}"
'''
module load IKMB
module load Plink/1.9

plink --threads 4 \
      --fam "!{fam}" \
      --map "!{ds_imp.PLINKdosage_map}" \
      --dosage "!{ds_imp.PLINKdosage_gz}" skip0=2 skip1=0 skip2=1 format=3 case-control-freqs \
      --covar "!{ds_release.dat_pca_evec}" \
      --covar-name "!{covar_name}" \
      --allow-no-sex \
      --ci 0.95 \
      --out "!{target}_chr!{chromosome}"

# Remove NAs
gawk '{ OFS="\t"; print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11 }' < "!{dosage_src}" | grep -v "NA" > "!{dosage_target}"

plink --threads 4 \
      --bfile "!{ds_stats.bim.baseName}" \
      --fam new_fam \
      --logistic hide-covar \
      --covar "!{ds_release.dat_pca_evec}" \
      --covar-name "!{covar_name}" \
      --allow-no-sex \
      --ci 0.95 \
      --chr !{chromosome} \
      --out "!{log_prefix}_tmp1"

plink --threads 4 \
      --bfile "!{ds_stats.bim.baseName}" \
      --fam new_fam \
      --assoc \
      --allow-no-sex \
      --ci 0.95 \
      --chr !{chromosome} \
      --out "!{log_prefix}_tmp2"
      
# Merge --dosage and --assoc results (don't forget to remove all the header lines when merging)
paste "!{log_prefix}_tmp1.assoc.logistic" "!{log_prefix}_tmp2.assoc" >merged

# Header line will also contain "genotyped" where it should say "INFO" instead.
# Merge process will take care of that.
gawk '{ OFS="\t"; print $1, $2, $3, $4, $19, $17, $18, "genotyped", $7, $8, $12 }' <merged > "!{log_prefix}.assoc.logistic"


'''
}

process merge_log_dos_results {
    input:
    file dosage_files from merge_log_dos_results_dosage.collect()
    file log_files from merge_log_dos_results_logistic.collect()
   
    output:
    file "${log_target}" into for_extract_rsq_variants_log
    file "${dos_target}" into for_extract_rsq_variants_dos
    
    shell: 
    log_target = "${params.disease_data_set_prefix_release_statistics}.imputed_chr${params.first_chr}-${params.last_chr}.assoc.logistic"
    dos_target = "${params.disease_data_set_prefix_release_statistics}.imputed_chr${params.first_chr}-${params.last_chr}.assoc.dosage"
'''
for chr in {!{params.first_chr}..!{params.last_chr}}
do
    LOGFILE="!{params.disease_data_set_prefix_release_statistics}.genotyped.chr${chr}.assoc.logistic"
    DOSAGEFILE="!{params.disease_data_set_prefix_release_statistics}.imputed_chr!{params.first_chr}-!{params.last_chr}.assoc.dosage.${chr}"
    
    # Add header line only on the first result file
    if [ $chr == !{params.first_chr} ]; then
        head -n 1 "$LOGFILE" | sed s/genotyped/INFO/ >"!{log_target}"
        head -n 1 "$DOSAGEFILE" >"!{dos_target}"
    fi
    
    # On every other data value, add everything else
    tail -n +2 "$LOGFILE" >>"!{log_target}"
    tail -n +2 "$DOSAGEFILE" >>"!{dos_target}"
done
'''
}


process extract_rsq_variants {
    input:
    file dosage from for_extract_rsq_variants_dos
    file logistic from for_extract_rsq_variants_log
    
    output:
    file "${rsq4}" into for_l95u95_rsq4
    file "${rsq4}.locuszoom" into for_l95u95_rsq4lz
    file "${rsq8}" into for_l95u95_rsqs8
    file "${rsq8}.locuszoom" into for_l95u95_rsq8lz
    
shell:
    rsq4 = "${params.disease_data_set_prefix_release_statistics}.genotyped.imputed_chr${params.first_chr}-${params.last_chr}.assoc.rsq0.4"
    rsq8 = "${params.disease_data_set_prefix_release_statistics}.genotyped.imputed_chr${params.first_chr}-${params.last_chr}.assoc.rsq0.8"
'''
python -c 'from Stats_helpers import *; extract_Rsq_variants(\
    assoc_logistic_input="!{logistic}", \
    assoc_dosage_input="!{dosage}", \
    assoc_merge_output_rsq0_4="!{rsq4}", \
    assoc_merge_output_rsq0_8="!{rsq8}" )'
'''
}

process addL95_U95 {

input:
    file rsq04    from for_l95u95_rsq4
    file rsq04_lz from for_l95u95_rsq4lz
    file rsq08    from for_l95u95_rsqs8
    file rsq08_lz from for_l95u95_rsq8lz

output:
    file rsq04  into for_qq_manhattan_rsq04
    file rsq08  into for_qq_manhattan_rsq08

shell:
pval2qchisq = SCRIPT_DIR + "/pval2qchisq_2.r"
'''
#!/bin/bash

function addL95_U95 {
    local ASSOC=$1
    awk '{ print $9, $11 }' ${ASSOC} >${ASSOC}.P.OR.tmp1
    R --no-save --args ${ASSOC}.P.OR.tmp1 ${ASSOC}.P.OR.tmp2 <"!{pval2qchisq}"
    gawk '{ print $1"\t"$2 }' ${ASSOC}.P.OR.tmp2 >${ASSOC}.P.OR.tmp3
    paste ${ASSOC} ${ASSOC}.P.OR.tmp3 >${ASSOC}.P.OR.tmp4
    mv ${ASSOC}.P.OR.tmp4 ${ASSOC}
}

addL95_U95 "!{rsq04}"
addL95_U95 "!{rsq04_lz}"
addL95_U95 "!{rsq08}"
addL95_U95 "!{rsq08_lz}"

'''
}

process qq_manhattan {
cpus 2

input:
file ds_staged from for_qq_manhattan_ds
file rsq04 from for_qq_manhattan_rsq04
file rsq08 from for_qq_manhattan_rsq08

file qqman from Channel.from(file(PYTHONPATH + "/qqman_dev1.r"))
file qqman2 from Channel.from(file(PYTHONPATH + "/qqman_dev2.r"))



shell:
ds = mapFileList(ds_staged)
logfile = ds.log
manhattan = PYTHONPATH + "/manhattan.r"
manhattan2 = PYTHONPATH + "/manhattan2.r"
qqplotp = PYTHONPATH + "/qqplotP.r"
qqplotp2 = PYTHONPATH + "/qqplotP_v2.r"
//qqman = SCRIPT_DIR + "/qqman_dev1.r"
//qqman2 = SCRIPT_DIR + "/qqman_dev2.r"
pval2chisq = PYTHONPATH + "/qqplot_assoc_pval2CHISQ.r"
rsq04_null = "${rsq04}_null"
rsq08_null = "${rsq08}_null"
rsq04_noxMHC = "${rsq04}_noxMHC"
rsq08_noxMHC = "${rsq08}_noxMHC"
null_variants = ANNOTATION_DIR + "/QQplots/" + params.QQplot_null_variants
snpexclude = BATCH_DIR + "/" + params.QQplot_SNPexcludeList

'''
NUM_CASES=$(grep -P '\\d+ are cases' !{logfile} | cut -d' ' -f4)
NUM_CONTROLS=$(grep -P '\\d+ are controls' !{logfile} | cut -d' ' -f8)

# ------------------- #
# -- null variants -- #
# ------------------- #
# Too few null variants available for QQplots because only high density
# regions and MHC/KIR were extracted here

# qqplotp, qqplotp2, pval2chisq, snpexclude

python -c "from Stats_helpers import *; extractQQplot_null_variants('!{null_variants}', '!{rsq04}', '!{rsq04}_null')"
python -c "from Stats_helpers import *; qqplot_null(file_PLINK_null='!{rsq04}_null', numof_cases=$NUM_CASES, numof_controls=$NUM_CONTROLS, qqplotp='!{qqplotp}', qqplotp2='!{qqplotp2}', pval2chisq='!{pval2chisq}', snpexclude='!{snpexclude}', collection_name='!{params.collection_name}',qqman='!{qqman}',qqman2='!{qqman2}')"
python -c "from Stats_helpers import *; extractQQplot_null_variants('!{null_variants}', '!{rsq08}', '!{rsq08}_null')"
python -c "from Stats_helpers import *; qqplot_null(file_PLINK_null='!{rsq08}_null', numof_cases=$NUM_CASES, numof_controls=$NUM_CONTROLS, qqplotp='!{qqplotp}', qqplotp2='!{qqplotp2}', pval2chisq='!{pval2chisq}', snpexclude='!{snpexclude}', collection_name='!{params.collection_name}',qqman='!{qqman}',qqman2='!{qqman2}')"

# ---------------------------- #
# -- filter for non-xMHC    -- #
# -- xMHC: chr6, [25,34[ Mb -- #
# ---------------------------- #

gawk '{ if (!($1 == 6 && ($3 >= 25000000 && $3 < 34000000))) print }' !{rsq04} >!{rsq04_noxMHC}
gawk '{ if (!($1 == 6 && ($3 >= 25000000 && $3 < 34000000))) print }' !{rsq08} >!{rsq08_noxMHC}
python -c "from Stats_helpers import*; qqplot_xMHC_noxMHC('!{rsq04}', '!{rsq04}_noxMHC', $NUM_CASES, $NUM_CONTROLS, qqplotp='!{qqplotp}', qqplotp2='!{qqplotp2}', pval2chisq='!{pval2chisq}', snpexclude='!{snpexclude}', collection_name='!{params.collection_name}',qqman='!{qqman}',qqman2='!{qqman2}')"
python -c "from Stats_helpers import*; qqplot_xMHC_noxMHC('!{rsq08}', '!{rsq08}_noxMHC', $NUM_CASES, $NUM_CONTROLS, qqplotp='!{qqplotp}', qqplotp2='!{qqplotp2}', pval2chisq='!{pval2chisq}', snpexclude='!{snpexclude}', collection_name='!{params.collection_name}',qqman='!{qqman}',qqman2='!{qqman2}')"

# -------------------------------- #
# -- manhattan plot for chr1-22 -- #
# -------------------------------- #

function manhattan_plot {
    local rsq=$1
    gawk ' $1 !~ /23/ { print $0 } ' ${rsq}      >${rsq}_tmp1
    gawk ' $1 !~ /24/ { print $0 } ' ${rsq}_tmp1 >${rsq}_tmp2
    gawk ' $1 !~ /25/ { print $0 } ' ${rsq}_tmp2 >${rsq}_tmp3
    gawk ' $1 !~ /26/ { print $0 } ' ${rsq}_tmp3 >${rsq}_noChr23
    R --slave --args ${rsq}_noChr23 !{params.collection_name} <!{manhattan}&
    R --slave --args ${rsq}_noChr23 !{params.collection_name} <!{manhattan2}&
    wait
}

manhattan_plot "!{rsq04}"
manhattan_plot "!{rsq08}"
manhattan_plot "!{rsq04}_excludeRegions"
manhattan_plot "!{rsq08}_excludeRegions"
'''
}


process convert_dosages {

input:
each chromosome from Channel.from(params.first_chr .. params.last_chr)

output:
file "$plink_target"

shell:
chrname = ds_imp  // vorsicht: join(Imputation_orig_dir, str(i)+"."+disease_data_set_suffix_release_imputed +".gz") 
flag_relatives_doubleid // disease_data_set_prefix_release + "_flag.relatives.doubleID.txt",
multiallelic_exclude 
target // plink_dosage + rsq0.4.chr + chr + .rs 
plink_target // plink_dosage + rsq0.4.chr + chr 
rs2chrpos = SCRIPT_DIR+"/awk_rs2CHRPOS_bimfiles.awk"
'''
module load Plink/1.9

plink --vcf "!{chrname}" --double-id --remove "!{flag_relatives_doubleid}" --keep "!{ds_stats.ped}" --exclude "!{multialellic_exclude}" --out "!{target}" --allow-no-sex --make-bed --threads 16

NEWTARG="!{target.fam.baseName}.tmp"
cp "!{target.fam}" "${NEWTARG}.fam"
cp "!{target.bed}" "${NEWTARG}.bed"
cp "!{target.log}" "${NEWTARG}.log"

LC_NUMERIC=POSIX gawk -f "!{rs2chrpos}" -- "!{target.bim}" >"${NEWTARG}.bim"

plink --bfile "${NEWTARG}" --exclude "!{multiallelic_exclude}" --out "!{plink_target}" --allow-no-sex --make-bed --threads 16
'''
}

/*
process merge_dosages {
    input:
    file dosages from for_merge_disages.collect()
   
    output:
    
   
    shell: 
    target = // plink dosage name + rsq0.4 + ".chr{first}-{last}
    dosage_basename = // ... 
'''
#!/bin/bash

module load Plink/1.9

truncate merge-list

for ((i=$((${params.first_chr} + 1)); i <= !{params.last_chr}; i++)); do
    echo "!{dosage_basename}.chr${i}.bed !{dosage_basename}.chr${i}.bim !{dosage_basename}.chr${i}.fam" >>merge-list 
done

plink --bfile "!{dosage_basename}.chr!{params.first_chr}" --merge-list merge-list --make-bed --out "!{target}_tmp" --allow-no-sex --threads 16
plink --bfile "!{target}_tmp" --pheno "!{ds_stats.fam_double_id_pheno}" --update-sex "!{ds_stats.fam_double_id_gender}" --make-bed --out "!{target}" --allow-no-sex --threads 16
gawk '{ print $2 }' "!{rsq08}" | tail -n +2 >"!{dosage_basename}.rsq0.8.rs.txt"
plink --bfile "!{rsq04.bim.baseName}" --extract "!{dosage_basename}.rsq0.8.rs.txt" --make-bed --out "!{rsq08.bim.baseName}" --allow-no-sex --threads 16

# Check number of individuals pre and post imputation
COUNT1=$(wc -l "!{dds_stats.fam}" | cut -d" " -f 1)
COUNT2=$(wc -l "!{rsq04.fam}" | cut -d" " -f 1)
if [ "$COUNT1" /= "$COUNT2" ]; then
    (>&2 echo "Abort: different number of individuals in fam file pre versus post imputation.")
    (>&2 echo "    Expected: $COUNT1 individuals in !{dds_stats.fam}")
    (>&2 echo "    Observed: $COUNT2 individuals in !{rsq04.fam}")
    exit 1
fi
'''
}

process cleanup_dataset {
    cpus 2
    
    input:
    rsq4
    rsq8
    dds_stats
    
shell:
'''
module load Plink/1.9

# the tmp version still has duplicates from genotyping (rs numbers) and imputation (chr.:...)
plink --bfile "!{rsq4.bim.baseName}" --bmerge "!{ds_stats.bim.baseName}'" --make-bed --out "!{rsq4.bim.baseName}_tmp" --allow-no-sex &
plink --bfile "!{rsq8.bim.baseName}" --bmerge "!{ds_stats.bim.baseName}'" --make-bed --out "!{rsq8.bim.baseName}_tmp" --allow-no-sex &
wait

# determine the duplicate SNPs from genotyping (rs numbers) and imputation (chr.:...)
# tmp version without duplicate genotyped and imputed SNPs - exclude imputed SNPs which are already genotyped
gawk '{ print $1":"$4 }' "!{rsq4.bim.baseName}_tmp.bim" | sort | uniq -d > "!{rsq4.bim.baseName}.excluded.duplicates.imputed"
plink --bfile "!{rsq4.bim.baseName}_tmp" --exclude "!{rsq4.bim.baseName}.excluded.duplicates.imputed" --make-bed --out "!{rsq4.bim.baseName}_tmp2" --allow-no-sex &
plink --bfile "!{rsq8.bim.baseName}_tmp" --exclude "!{rsq4.bim.baseName}.excluded.duplicates.imputed" --make-bed --out "!{rsq8.bim.baseName}_tmp2" --allow-no-sex &
wait

# determine monomorphic SNPs
plink --bfile "!{rsq4.bim.baseName}_tmp2" --freq --out "!{rsq4.bim.baseName}_tmp2_freq" --allow-no-sex --threads 16
gawk '{ if ($5 == 0) print $2 }' "!{rsq4.bim.baseName}_tmp2_freq.frq" > "!{rsq4.bim.baseName}.excluded.monomorphic"

# final version without duplicate genotyped and imputed SNPs - exclude imputed SNPs which are already genotyped
plink --bfile "!{rsq4.bim.baseName}_tmp2" --exclude "!{rsq4.bim.baseName}.excluded.monomorphic" --make-bed --out "!{rsq4.bim.baseName}.excluded.monomorphic" &
plink --bfile "!{rsq8.bim.baseName}_tmp2" --exclude "!{rsq4.bim.baseName}.excluded.monomorphic" --make-bed --out "!{rsq8.bim.baseName}.excluded.monomorphic" &
wait
plink --bfile "!{rsq4.bim.baseName}_tmp3" --not-chr 0 --make-bed --out "!{rsq04}" --allow-no-sex &
plink --bfile "!{rsq8.bim.baseName}_tmp3" --not-chr 0 --make-bed --out "!{rsq08}" --allow-no-sex &
wait

# another final version without duplicate genotyped and imputed SNPs - 
# exclude imputed SNPs which are already genotyped - but here chr:pos as SNPids for Locuszoom
awk '{ print $1":"$4 }' "!{rsq4.bim}" | sort | uniq -d > "!{rsq4.bim}.duplicates"
awk '{ print $1":"$4 }' "!{rsq8.bim}" | sort | uniq -d > "!{rsq8.bim}.duplicates"
plink --bfile "!{rsq4.bim.baseName}" --make-bed --out "!{rsq4.bim.baseName}.locuszoom.tmp" --allow-no-sex &
plink --bfile "!{rsq8.bim.baseName}" --make-bed --out "!{rsq8.bim.baseName}.locuszoom.tmp" --allow-no-sex &
wait

awk '{ print $1"\t"$1":"$4"\t"$3"\t"$4"\t"$5"\t"$6 }' "!{rsq4.bim}"> "!{rsq4.bim.baseName}.locuszoom.tmp.bim"
awk '{ print $1"\t"$1":"$4"\t"$3"\t"$4"\t"$5"\t"$6 }' "!{rsq8.bim}"> "!{rsq8.bim.baseName}.locuszoom.tmp.bim"
plink --bfile "!{rsq4.bim.baseName}.locuszoom.tmp" --exclude "!{rsq4.bim}.duplicates" --make-bed --out "!{rsq4.bim.baseName}.locuszoom" &
plink --bfile "!{rsq8.bim.baseName}.locuszoom.tmp" --exclude "!{rsq8.bim}.duplicates" --make-bed --out "!{rsq8.bim.baseName}.locuszoom" &
wait
'''
}


process definetti_plots {
    input:

    
shell:
'''
python -c 'from Stats_Immunochip import *; DeFinetti_plots(file_PLINK_geno_imp="!{file_PLINK_geno_imp}")'
'''
}

process plink_clumping {
input:
    rsq4_lz
    rsq8_lz
shell:
'''
python -c 'from Stats_Immunochip import *; PLINK_clumping("!{imp.bim.baseName}", "!{rsq4_lz.bim.baseName}", "!{rsq8_lz.bim.baseName}")'
'''
}

// TODO: This could be parallelized over SNPs (split/merge approach) 
process locuszoom {
    
shell:
'''
python -c 'from Stats_Immunochip import *; \
    locuszoom_run(snp_colnr=2, \
    file_suffix="_clump.clumped_all.noXMHC", \
    file_PLINK_assoc="!{rsq4_assoc}", \
    file_PLINK="!{rsq4}")'
'''
}


*/



