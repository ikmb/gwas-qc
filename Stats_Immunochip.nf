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
        case ~/.*_flag.relatives.doubleID.txt/:
            result = "relatives_doubleID"
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

// Originally from exclude_relatives.py but we don't have relatives anymore.
// TODO: Check if obsolete.
//ds_stats_input = ["${params.input_stats}.fam",
//              "${params.input_stats}.ped",
//              "${params.input_stats}_annotation.txt"].collect { fileExists(file(it)) }

// From "QCed" folder, output of QC pipeline
ds_input = ["${params.input}.bed",
         "${params.input}.bim",
         "${params.input}.fam",
         ].collect { fileExists(file(it)) }

// From Imputation server (1.vcf.gz, 2.vcf.gz, ...)
ds_imp_input = Channel.fromPath("${params.input_imp}/{[1-9],[1-2][0-9]}.vcf?gz", checkIfExists: true)


params.min_info_score = 0.3
params.first_chr = 1
params.last_chr = 22
params.disease_data_set_prefix_release_statistics = "${params.collection_name}_Release"

process preprocess_infofilter_dosage {
    tag "${params.collection_name}.${chrom}"
    cpus 2
    time { 48.h * task.attempt }
    errorStrategy { task.exitStatus == 128 ? 'retry' : 'terminate' }

    input:
    file vcfgz from ds_imp_input

    output:
    file "${chrom}.vcf.*" into vcfmaps_preproc
    file "${chrom}.INFO${params.min_info_score}.vcf.*" into infomaps_preproc
    set val(chrom), file("${chrom}.vcf.gz"), file("${chrom}.vcf.map"), file("${chrom}.INFO${params.min_info_score}.vcf.PLINKdosage.gz"), file("${chrom}.INFO${params.min_info_score}.vcf.PLINKdosage.map") into plinkdosage_gz_map, for_dosage_convert_gz_map

    shell:
    m = vcfgz =~ /(\d+).vcf.gz$/
    if(!m.getCount()) {
    println "No input files from '" + vcfgz + "': " + m
    } else {
    chrom = m[0][1]
    }
'''
module load Plink/1.9

AWK=mawk

# for each chromosome:
# vcf.gz -> vcf.map, uncompressed tmp.vcf

# Here we would need to uncompress the large vcfgz two times. This takes ages
# on seriously large files. A way to circumvent:
# 1. Create FIFO (uncompressed data will be put here later)
# 2. Connect some process that will consume uncompressed data
# 3. Start creating uncompressed data and duplicate the stream into the FIFO

# 1.
rm -f vcfgz_fifo
mkfifo vcfgz_fifo

# 2.
# vcf.gz -> vcf.map
(cat vcfgz_fifo | $AWK '$0!~/^#/ {NF=8; print $0}') >!{chrom}.vcf.map &

# 3.
# vcf.gz -> INFO.vcf.gz, uncompressed INFO.vcf.gz
bgzip -c -d !{vcfgz} | tee vcfgz_fifo \
    | $AWK '
    $0~/^#/ { 
        print $0
        next
    }
    {
        split($8, fields, /;|=/)
        for(i=1; i in fields; i++) {
            if(fields[i]=="INFO" || fields[i]=="DR2") {
                val=fields[i+1]
                break
            }
        }
        
        if(val >= !{params.min_info_score})
            {print $0}
        }' \
    | bgzip >!{chrom}.INFO!{params.min_info_score}.vcf.gz 

wait

tabix -p vcf !{chrom}.INFO!{params.min_info_score}.vcf.gz &

rm -f vcfgz_map_fifo
mkfifo vcfgz_map_fifo
# INFO.vcf.gz -> INFO.vcf.map

(cat vcfgz_map_fifo | $AWK '$0!~/^#/{NF=8; print $0}') >!{chrom}.INFO!{params.min_info_score}.vcf.map &

# AwkProg v2, more robust and faster than v1
cat <<'AwkProg2' >awkprog2.awk
BEGIN{
    printf "CHR BP SNP A1 A2 INFO"
}
{
    # skip header lines
    if (($1 ~ /^##/)) {
    } else {
        # Copy header line
        if (($1 ~ /^#CHROM/)) {
            for (i=10; i<=NF; i++)
                printf " "$i" "$i
            printf "\\n"
        } else {
            # Extract r^2 score from INFO field
            {
                split($8, fields, /(;|=)/)
                for(i=1; i in fields; i++) {
                    if(fields[i] == "INFO" || fields[i] == "DR2") {
                            info = fields[i+1]
                            break
                    }
                }
            }


            printf $1" "$2" "$1":"$2" "$4" "$5" "info
            for (i=10; i<=NF; i++) {
                split($i,array,":")
                # Missing fields are specified as ".", not splittable by ","
                if(array[4] == ".") {
                    printf " . . ."
                } else {
                    split(array[4],array_GP,",")
                    printf " "array_GP[1]" "array_GP[2]" "array_GP[3]
                }
            }
            printf "\\n"
        }
    }
}
AwkProg2

cat <<'AwkProg' >awkprog.awk
BEGIN {
printf "CHR BP SNP A1 A2 INFO"
}
{
        if (($1 ~ /^##/))  { # if comment line starting with hash character
        } else {
            if (($1 ~ /^#CHROM/)) { # if header line
                for (i=10; i<=NF; i++)
                    printf " "$i" "$i
                printf "\\n"
            } else {
                info = substr($8, match($8, /INFO=/)+5)
                printf $1" "$2" "$1":"$2" "$4" "$5" "info
                for (i=10; i<=NF; i++) {
                    split($i,array,":")
                    # Missing fields are specified as ".", not splittable by ","
                    if(array[4] == ".") {
                        printf " . . ."
                    } else {
                        split(array[4],array_GP,",")
                        printf " "array_GP[1]" "array_GP[2]" "array_GP[3]
                    }
                }
                printf "\\n"
            }
        }
}
AwkProg

rm -f plinkdosage
mkfifo plinkdosage

(cat plinkdosage    | tail -n +2 \
    | $AWK '{print $1, $3, 0, $2}' \
    | uniq) >!{chrom}.INFO!{params.min_info_score}.vcf.PLINKdosage.map &

bgzip -c -d !{chrom}.INFO!{params.min_info_score}.vcf.gz | tee vcfgz_map_fifo \
    | $AWK -f awkprog2.awk \
    | tee plinkdosage \
    | gzip >!{chrom}.INFO!{params.min_info_score}.vcf.PLINKdosage.gz

wait

'''
}



process preprocess_hrc_vcf_map {
publishDir params.output ?: '.', mode: 'copy', overwrite: true
    time { 24.h * task.attempt }
    tag "${params.collection_name}"
    input:
    file(mapfiles) from vcfmaps_preproc.collect()
    file(infomapfiles) from infomaps_preproc.collect()

    output:
    file "1-22.vcf.map"
    file "1-22.INFO${params.min_info_score}.vcf.map" into for_gen_sumstats_info

    shell:
    '''
    cp 1.vcf.map 1-22.vcf.map
    for ((i=2; i<=22; i++))
    do
        cat ${i}.vcf.map >> 1-22.vcf.map
    done

    cat 1.INFO!{params.min_info_score}.vcf.map >1-22.INFO0.3.vcf.map
    for ((i=2; i<=22; i++))
    do
        cat ${i}.INFO!{params.min_info_score}.vcf.map >> 1-22.INFO!{params.min_info_score}.vcf.map
    done
    '''
}

process extract_genotyped_variants {
	cpus 4
    tag "${params.collection_name}"
    errorStrategy { task.exitStatus == 128 ? 'retry' : 'terminate' }

	input:
//	file ds_stats_staged from Channel.from(ds_stats_input).collect()
	file ds_staged from Channel.from(ds_input).collect()

	output:
	file "${target}_fid_iid.fam" into for_plink_dosage_logistic_fam, for_cleanup_dataset_fam, for_merge_dosages_fam, for_clumping_fam
	file "${target}.{bim,bed,fam,log,dat.pca.evec}" into for_plink_dosage_logistic_ds, for_qq_manhattan_ds, for_convert_dosages_stats, for_merge_dosages_stats, for_cleanup_dataset_stats, for_clumping_ds
	file "${target}.fam.double-id.pheno" into for_merge_dosages_pheno
	file "${target}.fam.double-id.gender" into for_merge_dosages_gender

    shell:
//    ds_stats = mapFileList(ds_stats_staged)
    ds = mapFileList(ds_staged)
    target = "dummy${params.disease_data_set_prefix_release_statistics}"

'''
module load Plink/1.9
AWK=mawk

$AWK '{ print $1, $2, $6 }' "!{ds.fam}" >"!{ds.fam}.single-id.pheno"
plink --bfile "!{ds.bim.baseName}" \
	--threads 4 \
	--pheno "!{ds.fam}.single-id.pheno" \
	--make-bed --allow-no-sex \
	--out "!{target}" || true


#	--remove "${ds.relatives}" \
#	  --keep-fam "${ds_stats.fam}"

$AWK '{ print $1\"_\"$2, $1\"_\"$2, $6 }' "!{target}.fam" >"!{target}.fam.double-id.pheno"
$AWK '{ print $1\"_\"$2, $1\"_\"$2, $5 }' "!{target}.fam" >"!{target}.fam.double-id.gender"
$AWK '{ print $1\"_\"$2, $1\"_\"$2, $3, $4, $5, $6 }' "!{target}.fam" | TMPDIR="." sponge "!{target}_fid_iid.fam"

cp !{params.covar} covar
$AWK 'FILENAME~/fam$/{samples[$2]=1}(FILENAME~/covar$/ && FNR==1) {print $0} (FILENAME~/covar$/ && ($2 in samples || FNR!=1)){$1=$1"_"$2;$2=$1;print $0}' !{target}.fam covar >"!{target}.dat.pca.evec"
'''
}

// Steps 1.1 and 1.2: per-chromosome dosage and logistic calculations,
// results will be merged in merge_log_dos_results
process plink_dosage_logistic {
    validExitStatus 0,128
    cpus 4
    tag "${params.collection_name}.$chromosome"

    input:
    set val(chromosome), file(vcfgz), file(vcfmap), file(dosage_gz), file(dosage_map) from plinkdosage_gz_map
    file (fam:"new_fam") from for_plink_dosage_logistic_fam
	file ds_stats_staged from for_plink_dosage_logistic_ds
//	file ds_imp_staged from Channel.from(ds_imp_input).collect()
	file ds_release_staged from Channel.from(ds_input).collect()
//	file ds_stats_orig_staged from Channel.from(ds_stats_input).collect()

    output:
//    file "${target}_chr${chromosome}.assoc.dosage" into merge_log_dos_results_dosage
    file "${dosage_target}" into merge_log_dos_results_dosage
    file "${log_prefix}.assoc.logistic" into merge_log_dos_results_logistic

    shell:
    ds_stats = mapFileList(ds_stats_staged)
    //ds_imp = mapFileList(ds_imp_staged)
    ds_release = mapFileList(ds_release_staged)
//    ds_stats_orig = mapFileList(ds_stats_orig_staged)

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

echo "DS_STATS: !{ds_stats}"


plink --threads 4 \
      --fam "!{fam}" \
      --map "!{dosage_map}" \
      --dosage "!{dosage_gz}" skip0=2 skip1=0 skip2=1 format=3 case-control-freqs \
      --covar "!{ds_stats.dat_pca_evec}" \
      --covar-name "!{covar_name}" \
      --allow-no-sex \
      --ci 0.95 \
      --out "!{target}_chr!{chromosome}"

# Remove NAs
gawk '{ OFS="\t"; print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11 }' < "!{dosage_src}" | grep -v "NA" > "!{dosage_target}"

#      --fam new_fam
plink --threads 4 \
      --fam "!{fam}" \
      --bfile "!{ds_stats.bim.baseName}" \
      --logistic hide-covar \
      --covar "!{ds_stats.dat_pca_evec}" \
      --covar-name "!{covar_name}" \
      --allow-no-sex \
      --ci 0.95 \
      --chr !{chromosome} \
      --out "!{log_prefix}_tmp1"

# --fam new_fam
plink --threads 4 \
      --fam "!{fam}" \
      --bfile "!{ds_stats.bim.baseName}" \
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
publishDir params.output ?: '.', mode: 'copy', overwrite: true
    validExitStatus 0,128
    tag "${params.collection_name}"
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
    tag "${params.collection_name}"
    input:
    file dosage from for_extract_rsq_variants_dos
    file logistic from for_extract_rsq_variants_log

    output:
    file "${rsq3}" into for_l95u95_rsq3, for_gen_sumstats_rsq3
    file "${rsq3}.locuszoom" into for_l95u95_rsq3lz
    file "${rsq8}" into for_l95u95_rsq8, for_gen_sumstats_rsq8
    file "${rsq8}.locuszoom" into for_l95u95_rsq8lz

shell:
    rsq3 = "${params.disease_data_set_prefix_release_statistics}.genotyped.imputed_chr${params.first_chr}-${params.last_chr}.assoc.rsq0.3"
    rsq8 = "${params.disease_data_set_prefix_release_statistics}.genotyped.imputed_chr${params.first_chr}-${params.last_chr}.assoc.rsq0.8"
'''
python -c 'from Stats_helpers import *; extract_Rsq_variants(\
    assoc_logistic_input="!{logistic}", \
    assoc_dosage_input="!{dosage}", \
    assoc_merge_output_rsq0_3="!{rsq3}", \
    assoc_merge_output_rsq0_8="!{rsq8}" )'
'''
}

process addL95_U95 {
    tag "${params.collection_name}"
input:
    file rsq03    from for_l95u95_rsq3
    file rsq03_lz from for_l95u95_rsq3lz
    file rsq08    from for_l95u95_rsq8
    file rsq08_lz from for_l95u95_rsq8lz

output:
    file rsq03  into for_qq_manhattan_rsq03
    file rsq08  into for_qq_manhattan_rsq08, for_merge_dosages_rsq08
    file rsq03_lz into for_clumping_rsq03, for_lz_lz03
    file rsq08_lz into for_clumping_rsq08

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

addL95_U95 "!{rsq03}"
addL95_U95 "!{rsq03_lz}"
addL95_U95 "!{rsq08}"
addL95_U95 "!{rsq08_lz}"

'''
}

process qq_manhattan {
tag "${params.collection_name}"
publishDir params.output ?: '.', mode: 'copy', overwrite: true
validExitStatus 0,1
time 10.h
cpus 2

input:
file ds_staged from for_qq_manhattan_ds
file rsq03 from for_qq_manhattan_rsq03
file rsq08 from for_qq_manhattan_rsq08

file qqman from Channel.from(file(PYTHONPATH + "/qqman_dev1.r"))
file qqman2 from Channel.from(file(PYTHONPATH + "/qqman_dev2.r"))

output:
file "*.jpg"
file "*.pdf"

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
rsq03_null = "${rsq03}_null"
rsq08_null = "${rsq08}_null"
rsq03_noxMHC = "${rsq03}_noxMHC"
rsq08_noxMHC = "${rsq08}_noxMHC"
null_variants = ANNOTATION_DIR + "/QQplots/" + params.QQplot_null_variants
snpexclude = BATCH_DIR + "/" + params.QQplot_SNPexcludeList

if(!file(null_variants).exists()) {
    null_variants=file("/dev/null")
}

if(!file(snpexclude).exists()) {
    snpexclude=file("/dev/null")
}

'''
NUM_CASES=$(grep -P '\\d+ are cases' !{logfile} | cut -d' ' -f4)
NUM_CONTROLS=$(grep -P '\\d+ are controls' !{logfile} | cut -d' ' -f8)

# ------------------- #
# -- null variants -- #
# ------------------- #
# Too few null variants available for QQplots because only high density
# regions and MHC/KIR were extracted here

# qqplotp, qqplotp2, pval2chisq, snpexclude

python -c "from Stats_helpers import *; extractQQplot_null_variants('!{null_variants}', '!{rsq03}', '!{rsq03}_null')"
python -c "from Stats_helpers import *; qqplot_null(file_PLINK_null='!{rsq03}_null', numof_cases=$NUM_CASES, numof_controls=$NUM_CONTROLS, qqplotp='!{qqplotp}', qqplotp2='!{qqplotp2}', pval2chisq='!{pval2chisq}', snpexclude='!{snpexclude}', collection_name='!{params.collection_name}',qqman='!{qqman}',qqman2='!{qqman2}')"&
python -c "from Stats_helpers import *; extractQQplot_null_variants('!{null_variants}', '!{rsq08}', '!{rsq08}_null')"
python -c "from Stats_helpers import *; qqplot_null(file_PLINK_null='!{rsq08}_null', numof_cases=$NUM_CASES, numof_controls=$NUM_CONTROLS, qqplotp='!{qqplotp}', qqplotp2='!{qqplotp2}', pval2chisq='!{pval2chisq}', snpexclude='!{snpexclude}', collection_name='!{params.collection_name}',qqman='!{qqman}',qqman2='!{qqman2}')"&
wait

# ---------------------------- #
# -- filter for non-xMHC    -- #
# -- xMHC: chr6, [25,34[ Mb -- #
# ---------------------------- #

gawk '{ if (!($1 == 6 && ($3 >= 25000000 && $3 < 34000000))) print }' !{rsq03} >!{rsq03_noxMHC}
gawk '{ if (!($1 == 6 && ($3 >= 25000000 && $3 < 34000000))) print }' !{rsq08} >!{rsq08_noxMHC}
python -c "from Stats_helpers import*; qqplot_xMHC_noxMHC('!{rsq03}', '!{rsq03}_noxMHC', $NUM_CASES, $NUM_CONTROLS, qqplotp='!{qqplotp}', qqplotp2='!{qqplotp2}', pval2chisq='!{pval2chisq}', snpexclude='!{snpexclude}', collection_name='!{params.collection_name}',qqman='!{qqman}',qqman2='!{qqman2}')"&
python -c "from Stats_helpers import*; qqplot_xMHC_noxMHC('!{rsq08}', '!{rsq08}_noxMHC', $NUM_CASES, $NUM_CONTROLS, qqplotp='!{qqplotp}', qqplotp2='!{qqplotp2}', pval2chisq='!{pval2chisq}', snpexclude='!{snpexclude}', collection_name='!{params.collection_name}',qqman='!{qqman}',qqman2='!{qqman2}')"&
wait

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

manhattan_plot "!{rsq03}"
manhattan_plot "!{rsq08}"
manhattan_plot "!{rsq03}_excludeRegions"
manhattan_plot "!{rsq08}_excludeRegions"
'''
}


process convert_dosages {
validExitStatus 0,128
tag "${params.collection_name}.$chromosome"

input:
//each chromosome from Channel.from(params.first_chr .. params.last_chr)
set val(chromosome), file(vcfgz), file(vcfmap), file(dosage_gz), file(dosage_map) from for_dosage_convert_gz_map
//file ds_imp_staged from Channel.from(ds_imp_input).collect()
file ds_release_staged from Channel.from(ds_input).collect()
file ds_stats_staged from for_convert_dosages_stats
//file ds_stats_orig_staged from Channel.from(ds_stats_input).collect()

output:
// file "$plink_target.{bim,bed,fam,log}" into for_merge_dosages_rsq03//,for_merge_dosages
file "$target.{bim,bed,fam,log}" into for_merge_dosages_rsq03


shell:
//ds_stats_orig = mapFileList(ds_stats_orig_staged)
ds_release = mapFileList(ds_release_staged)
ds_stats = mapFileList(ds_stats_staged)
flag_relatives_doubleid = ds_release.relatives_doubleID
multiallelic_exclude = params.Multiallelic_SNPexcludeList

target = "${params.disease_data_set_prefix_release_statistics}.genotyped.imputed.rsq${params.min_info_score}.chr${chromosome}"
plink_target = "${params.disease_data_set_prefix_release_statistics}.genotyped.imputed.rsq${params.min_info_score}.chr${chromosome}.chrpos"

rs2chrpos = SCRIPT_DIR+"/awk_rs2CHRPOS_bimfiles.awk"
'''
module load Plink/1.9
TMPDIR=.

### Convert VCF to Plink which will include multiallelics and same-id-different-position SNPs
plink --vcf "!{vcfgz}" --double-id --make-bed --out "!{target}_with_multiallelics" --allow-no-sex

### Find multi-position rs ids
cut -f2 "!{target}_with_multiallelics.bim" | sort | uniq -d >multi-position.rs
awk 'NR==FNR{rsids[$0];next} {if($2 in rsids){$2=$1":"$4} print}' multi-position.rs "!{target}_with_multiallelics.bim" >new.bim
cp new.bim "!{target}_with_multiallelics.bim"
rm new.bim

### Convert Rs to chr:pos so we can identify different Rs IDs on the same position, i.e. multiallelics
LC_NUMERIC=POSIX gawk -f "!{rs2chrpos}" -- "!{target}_with_multiallelics.bim" >"!{target}_chrpos_with_multiallelics.bim"
<"!{target}_chrpos_with_multiallelics.bim" awk '{print $2}' | TMPDIR=. sort | uniq -d >"!{target}.multiallelics"

# plink --bim "!{target}_chrpos_with_multiallelics.bim" --bed "!{target}_with_multiallelics.bed" --fam "!{target}_with_multiallelics.fam" --exclude "!{target}.multiallelics" --make-bed --out "!{plink_target}"

### Find the variant names with the listed chr:pos multiallelics, so we can remove them
gawk 'NR==FNR{ids[$0];next} {if( ($1 ":" $4) in ids) {print $2}}' "!{target}.multiallelics" "!{target}_with_multiallelics.bim" >"!{target}.multiallelics.rs"

# Remove multiallelics
if [ -e "!{multiallelic_exclude}" ]; then
    cat "!{multiallelic_exclude}" >>"!{target}.multiallelics.rs"
fi

plink --bfile "!{target}_with_multiallelics" --exclude "!{target}.multiallelics.rs" --make-bed --out "!{target}"
# Replace '.'-Variants with chr:pos names
mv "!{target}.bim" "!{target}.bim_with_dots"
gawk '{if($2==".") {$2=$1 ":" $4} print $0}' <"!{target}.bim_with_dots" >"!{target}.bim"


#NEWTARG="!{target}.tmp"
#ln -s "!{target}.fam" "${NEWTARG}.fam"
#ln -s "!{target}.bed" "${NEWTARG}.bed"
#ln -s "!{target}.log" "${NEWTARG}.log"

# LC_NUMERIC=POSIX gawk -f "!{rs2chrpos}" -- "!{target}.bim" >"${NEWTARG}.bim"

#if [ -e "!{multiallelic_exclude}" ]; then
#    plink --bfile "${NEWTARG}" --exclude "!{multiallelic_exclude}" --out "!{plink_target}" --allow-no-sex --make-bed --threads 16
#else
#    plink --bfile "${NEWTARG}" --out "!{plink_target}" --allow-no-sex --make-bed --threads 16
#fi

'''
}


process merge_dosages {
    errorStrategy 'retry'
    memory { 20.GB * task.attempt }
    time { 8.h * task.attempt }
    tag "${params.collection_name}"
    input:
    //file dosages from for_merge_dosages.collect()
    file rsq03_staged from for_merge_dosages_rsq03.collect()
    file rsq08 from for_merge_dosages_rsq08
    file ds_stats_staged from for_merge_dosages_stats
    file pheno from for_merge_dosages_pheno
    file gender from for_merge_dosages_gender
    file (fam: "${params.disease_data_set_prefix_release_statistics}.fam") from for_merge_dosages_fam

    output:
    file "$target.{bim,bed,fam,log}" into for_cleanup_dataset_rsq03
    file "$rsq08_target_name.{bim,bed,fam,log}" into for_cleanup_dataset_rsq08

    shell:
    ds_stats = mapFileList(ds_stats_staged)
    rsq03 = mapFileList(rsq03_staged)
//    rsq08 = "${params.disease_data_set_prefix_release_statistics}.genotyped.imputed.rsq0.8.chr${params.first_chr}-${params.last_chr}"

    target = "${params.disease_data_set_prefix_release_statistics}.genotyped.imputed.rsq${params.min_info_score}.chr${params.first_chr}-${params.last_chr}"// plink dosage name + rsq0.4 + ".chr{first}-{last}
    dosage_basename = "${params.disease_data_set_prefix_release_statistics}.genotyped.imputed.rsq${params.min_info_score}"
    rsq08_target_name = "${params.disease_data_set_prefix_release_statistics}.genotyped.imputed.rsq0.8.chr${params.first_chr}-${params.last_chr}"
//    #
//
'''
#!/bin/bash

module load Plink/1.9

rm -f merge-list

for ((i=$((!{params.first_chr})); i <= !{params.last_chr}; i++)); do
    echo "!{dosage_basename}.chr${i}.bed !{dosage_basename}.chr${i}.bim !{dosage_basename}.chr${i}.fam" >>merge-list
done

plink --bfile "!{dosage_basename}.chr!{params.first_chr}" \
    --merge-list merge-list \
    --memory 16000\
    --make-bed --out "!{target}_tmp" --allow-no-sex

#TRIES=0
#RETRY=1
#while [ "$RETRY" -eq 1 ]
#do
#    touch "multiallelic-excludes"
#    plink --bfile "!{dosage_basename}.chr!{params.first_chr}" \
#        --merge-list merge-list \
#        --exclude "multiallelic-excludes" \
#        --make-bed --out "!{target}_tmp" --allow-no-sex \
#        || true
#
#    (( TRIES++ ))
#    if [ "$TRIES" -le 10 ] && [ -s "!{target}_tmp-merge.missnp" ]; then
#        cat "!{target}_tmp-merge.missnp" >>multiallelic-excludes
#        rm "!{target}_tmp-merge.missnp"
#    elif ! [ -e "!{target}_tmp-merge.missnp" ]; then
#        RETRY=0
#        echo "All $(wc -l multiallelic-excludes) multiallelics removed."
#    else
#        RETRY=0
#        echo "Error: could not remove all multiallelics - check input data" >>/dev/stderr
#        exit 1
#    fi
#done

plink --bfile "!{target}_tmp" \
      --pheno "!{pheno}" \
      --update-sex "!{gender}" \
      --memory 16000 \
      --make-bed --out "!{target}" --allow-no-sex 
gawk '{ print $2 }' "!{rsq08}" | tail -n +2 >"!{dosage_basename}.rsq0.8.rs.txt"
plink --bfile "!{target}" \
      --extract "!{dosage_basename}.rsq0.8.rs.txt" \
      --memory 16000 \
      --make-bed --out "!{rsq08_target_name}" --allow-no-sex

# Check number of individuals pre and post imputation
COUNT1=$(wc -l "!{fam}" | cut -d" " -f 1)
COUNT2=$(wc -l "!{rsq03.fam}" | cut -d" " -f 1)
if [ "$COUNT1" -ne "$COUNT2" ]; then
    (>&2 echo "Abort: different number of individuals in fam file pre versus post imputation.")
    (>&2 echo "    Expected: $COUNT1 individuals in !{ds_stats.fam}")
    (>&2 echo "    Observed: $COUNT2 individuals in !{rsq03.fam}")
#    exit 1
fi
'''
}

/*for_cleanup_dataset_rsq03.subscribe { println "Channel rsq3: $it\n" }*/
/*for_cleanup_dataset_rsq08.subscribe { println "Channel rsq8: $it\n" }*/
/*for_cleanup_dataset_fam.subscribe { println "Channel fam: $it\n" }*/
/*for_cleanup_dataset_stats.subscribe { println "Channel stats: $it" }*/

/*for_definetti_rsq03 = Channel.create()*/
/*for_definetti_rsq08 = Channel.create()*/
/*for_clump_rsq03_ds = Channel.create()*/
/*for_clump_rsq08_ds = Channel.create()*/

/*
process cleanup_dataset {
    input:
    file rsq03_staged from for_cleanup_dataset_rsq03
    file rsq08_staged from for_cleanup_dataset_rsq08
    file (fam: "${params.disease_data_set_prefix_release_statistics}.fam") from for_cleanup_dataset_fam
    file ds_stats_staged from for_cleanup_dataset_stats

    output:
    file 'dummy' into for_definetti_rsq03, for_definetti_rsq08, for_clump_rsq03_ds, for_clump_rsq08_ds

    shell:
    rsq3 = mapFileList(rsq03_staged)
    rsq8 = mapFileList(rsq08_staged)
    ds_stats = mapFileList(ds_stats_staged)
    '''
    echo "DS_STATS: !{ds_stats}"
    echo "RSQ3: !{rsq3}"
    echo "RSQ8: !{rsq8}"
    '''
}
*/

process cleanup_dataset {
memory 35.GB
time 12.h
publishDir params.output ?: '.', mode: 'copy', overwrite: true
tag "${params.collection_name}"
    cpus 2
    stageInMode 'copy' // some files need to be overwritten, cannot stage them as links

    input:
    file rsq03_staged from for_cleanup_dataset_rsq03
    file rsq08_staged from for_cleanup_dataset_rsq08
    file (fam: "${params.disease_data_set_prefix_release_statistics}.fam") from for_cleanup_dataset_fam
    file ds_stats_staged from for_cleanup_dataset_stats

    output:
//    file "${rsq3base}.{bim,bed,fam,log}" into for_definetti_rsq03
    set file("${rsq3base}.bim"), file("${rsq3base}.bed"), file ("${rsq3base}.fam"), file ("${rsq3base}.log") into for_definetti_rsq03, for_gen_sumstats
    set file("${rsq8base}.bim"), file("${rsq8base}.bed"), file ("${rsq8base}.fam"), file ("${rsq8base}.log") into for_definetti_rsq08
//    file "${rsq8base}.{bim,bed,fam,log}" into for_definetti_rsq08
    set file("${rsq3.bim.baseName}.locuszoom.bim"),file("${rsq3.bim.baseName}.locuszoom.bed"),file("${rsq3.bim.baseName}.locuszoom.fam"),file("${rsq3.bim.baseName}.locuszoom.log") into for_clump_rsq03_ds, for_lz_rsq03
    set file("${rsq8.bim.baseName}.locuszoom.bim"),file("${rsq8.bim.baseName}.locuszoom.bed"),file("${rsq8.bim.baseName}.locuszoom.fam"),file("${rsq8.bim.baseName}.locuszoom.log") into for_clump_rsq08_ds

    shell:
    
    ds_stats = mapFileList(ds_stats_staged)
    rsq3 = mapFileList(rsq03_staged)
    rsq8 = mapFileList(rsq08_staged)
    rsq3base = rsq3.bim.baseName
    rsq8base = rsq8.bim.baseName
    
    
'''
module load Plink/1.9
echo "DS_STATS: !{ds_stats}"
echo "rsq3: !{rsq3}"
echo "RSQ8: !{rsq8}"
echo "DS_STATS_staged: !{ds_stats_staged}"
echo "rsq3_staged: !{rsq03_staged}"
echo "rsq8_staged: !{rsq08_staged}"

# if things don't work out, remove snps with more than 2 alleles
plink --memory 16000 --bfile "!{rsq3.bim.baseName}" --bmerge "!{ds_stats.bim.baseName}" --make-bed --out "!{rsq3.bim.baseName}_tmp" --allow-no-sex || true
if [ -e "!{rsq3.bim.baseName}_tmp-merge.missnp" ]; then
    plink --memory 16000 --bfile "!{rsq3.bim.baseName}" --exclude !{rsq3.bim.baseName}_tmp-merge.missnp --make-bed --allow-no-sex --out rsq3 &
    plink --memory 16000 --bfile "!{ds_stats.bim.baseName}" --exclude !{rsq3.bim.baseName}_tmp-merge.missnp --make-bed --allow-no-sex --out ds &
    wait
    plink --memory 16000 --bfile "!{rsq8.bim.baseName}" --exclude !{rsq3.bim.baseName}_tmp-merge.missnp --make-bed --allow-no-sex --out rsq8
    plink --memory 16000 --bfile rsq3 --bmerge ds --make-bed --out "!{rsq3.bim.baseName}_tmp" --allow-no-sex &
    plink --memory 16000 --bfile rsq8 --bmerge ds --make-bed --out "!{rsq8.bim.baseName}_tmp" --allow-no-sex &
    wait
else
    plink --memory 16000 --bfile "!{rsq8.bim.baseName}" --bmerge "!{ds_stats.bim.baseName}" --make-bed --out "!{rsq8.bim.baseName}_tmp" --allow-no-sex
fi

# determine the duplicate SNPs from genotyping (rs numbers) and imputation (chr.:...)
# tmp version without duplicate genotyped and imputed SNPs - exclude imputed SNPs which are already genotyped
gawk '{ print $1":"$4 }' "!{rsq3.bim.baseName}_tmp.bim" | sort | uniq -d > "!{rsq3.bim.baseName}.excluded.duplicates.imputed"
plink --memory 16000 --bfile "!{rsq3.bim.baseName}_tmp" --exclude "!{rsq3.bim.baseName}.excluded.duplicates.imputed" --make-bed --out "!{rsq3.bim.baseName}_tmp2" --allow-no-sex &
plink --memory 16000 --bfile "!{rsq8.bim.baseName}_tmp" --exclude "!{rsq3.bim.baseName}.excluded.duplicates.imputed" --make-bed --out "!{rsq8.bim.baseName}_tmp2" --allow-no-sex &
wait

# determine monomorphic SNPs
plink --memory 16000 --bfile "!{rsq3.bim.baseName}_tmp2" --freq --out "!{rsq3.bim.baseName}_tmp2_freq" --allow-no-sex --threads 16
gawk '{ if ($5 == 0) print $2 }' "!{rsq3.bim.baseName}_tmp2_freq.frq" > "!{rsq3.bim.baseName}.excluded.monomorphic"

# final version without duplicate genotyped and imputed SNPs - exclude imputed SNPs which are already genotyped
plink --memory 16000 --bfile "!{rsq3.bim.baseName}_tmp2" --exclude "!{rsq3.bim.baseName}.excluded.monomorphic" --make-bed --out "!{rsq3.bim.baseName}_tmp3" &
plink --memory 16000 --bfile "!{rsq8.bim.baseName}_tmp2" --exclude "!{rsq3.bim.baseName}.excluded.monomorphic" --make-bed --out "!{rsq8.bim.baseName}_tmp3" &
wait
plink --memory 16000 --bfile "!{rsq3.bim.baseName}_tmp3" --not-chr 0 --make-bed --out "!{rsq3.bim.baseName}" --allow-no-sex &
plink --memory 16000 --bfile "!{rsq8.bim.baseName}_tmp3" --not-chr 0 --make-bed --out "!{rsq8.bim.baseName}" --allow-no-sex &
wait

# another final version without duplicate genotyped and imputed SNPs -
# exclude imputed SNPs which are already genotyped - but here chr:pos as SNPids for Locuszoom
awk '{ print $1":"$4 }' "!{rsq3.bim}" | sort | uniq -d > "!{rsq3.bim}.duplicates"
awk '{ print $1":"$4 }' "!{rsq8.bim}" | sort | uniq -d > "!{rsq8.bim}.duplicates"
plink --memory 16000 --bfile "!{rsq3.bim.baseName}" --make-bed --out "!{rsq3.bim.baseName}.locuszoom.tmp" --allow-no-sex &
plink --memory 16000 --bfile "!{rsq8.bim.baseName}" --make-bed --out "!{rsq8.bim.baseName}.locuszoom.tmp" --allow-no-sex &
wait

awk '{ print $1"\t"$1":"$4"\t"$3"\t"$4"\t"$5"\t"$6 }' "!{rsq3.bim}"> "!{rsq3.bim.baseName}.locuszoom.tmp.bim"
awk '{ print $1"\t"$1":"$4"\t"$3"\t"$4"\t"$5"\t"$6 }' "!{rsq8.bim}"> "!{rsq8.bim.baseName}.locuszoom.tmp.bim"
plink --memory 16000 --bfile "!{rsq3.bim.baseName}.locuszoom.tmp" --exclude "!{rsq3.bim}.duplicates" --make-bed --out "!{rsq3.bim.baseName}.locuszoom" &
plink --memory 16000 --bfile "!{rsq8.bim.baseName}.locuszoom.tmp" --exclude "!{rsq8.bim}.duplicates" --make-bed --out "!{rsq8.bim.baseName}.locuszoom" &
wait

'''
}

process generate_sumstats {
publishDir params.output ?: '.', mode: 'copy', overwrite: true
tag "${params.collection_name}"
cpus 2
time {4.h * task.attempt }
    input:
    set file(bim), file(bed), file(fam), file(log) from for_gen_sumstats
    file(assoc03) from for_gen_sumstats_rsq3
    file(assoc08) from for_gen_sumstats_rsq8
    file(info) from for_gen_sumstats_info

    output:
    file "${params.collection_name}.INFO0.3.sumstats.gz"
    file "${params.collection_name}.INFO0.8.sumstats.gz"
shell:
'''

plink2sumstats.py "!{log}" "!{assoc03}" >"!{params.collection_name}.INFO0.3.sumstats.tmp" &
plink2sumstats.py "!{log}" "!{assoc08}" >"!{params.collection_name}.INFO0.8.sumstats.tmp" &
wait

(vcf2sumstats.pl "!{info}" "!{params.collection_name}.INFO0.3.sumstats.tmp" | gzip >"!{params.collection_name}.INFO0.3.sumstats.gz") &
(vcf2sumstats.pl "!{info}" "!{params.collection_name}.INFO0.8.sumstats.tmp" | gzip >"!{params.collection_name}.INFO0.8.sumstats.gz") &
wait
'''

}


process definetti_plots {
tag "${params.collection_name}"
publishDir params.output ?: '.', mode: 'copy', overwrite: true
    input:
    file rsq03_staged from for_definetti_rsq03
    file rsq08_staged from for_definetti_rsq08

    output:
    file "*.jpg"

shell:
rsq03base = mapFileList(rsq03_staged).bim.baseName
rsq08base = mapFileList(rsq08_staged).bim.baseName
definetti = SCRIPT_DIR + "/DeFinetti_hardy.r"
'''
module load Plink/1.9

plink --bfile "!{rsq03base}" --hardy --out "!{rsq03base}_hardy" --hwe 0.0 --threads 4
R --slave --args "!{rsq03base}_hardy.hwe" "!{rsq03base}_controls_DeFinetti" "!{rsq03base}_cases_DeFinetti" "!{rsq03base}_cases_controls_DeFinetti" <"!{definetti}"

plink --bfile "!{rsq08base}" --hardy --out "!{rsq08base}_hardy" --hwe 0.0 --threads 4
R --slave --args "!{rsq08base}_hardy.hwe" "!{rsq08base}_controls_DeFinetti" "!{rsq08base}_cases_DeFinetti" "!{rsq08base}_cases_controls_DeFinetti" <"!{definetti}"
'''
}


process plink_clumping {
tag "${params.collection_name}"
memory 24.GB
input:
    file ds_imp_staged from for_clumping_ds
//    file ds_imp_fam from Channel.from(ds_stats_input).collect()
    file rsq3_lz from for_clumping_rsq03
    file rsq8_lz from for_clumping_rsq08
    file rsq3_staged from for_clump_rsq03_ds
    file rsq8_staged from for_clump_rsq08_ds

    // python -c 'from Stats_Immunochip import *; PLINK_clumping("!{imp.fam.baseName}", "!{rsq3_lz}", "!{rsq8_lz}")'

output:
    file "${rsq3_lz}_clump.clumped*" into for_lz_clump03


shell:
imp = mapFileList(ds_imp_staged)
rsq3 = mapFileList(rsq3_staged)
rsq8 = mapFileList(rsq8_staged)
'''
module load Plink/1.9


# replace runs of tabs and spaces with a single space, take column 2 (sample names) and list duplicates
echo "Detecting duplicates..."
<"!{rsq3.bim}" tr -s '\t ' ' ' | cut -f2 -d' ' | uniq -d >rsq3-duplicates
<"!{rsq8.bim}" tr -s '\t ' ' ' | cut -f2 -d' ' | uniq -d >rsq8-duplicates

# output those records that are not in the 'duplicates' list
echo "Checking for and removing duplicates from locuszoom tables..."
if [ ! -s rsq3-duplicates ]; then
    awk 'NR==FNR{duplicates[$0];next} {f=!($2 in duplicates)} f' rsq3-duplicates "!{rsq3_lz}" >rsq3_lz
    plink --bfile "!{rsq3.bim.baseName}" --memory 10000 --threads 1 --exclude rsq3-duplicates --make-bed --out rsq3 --allow-no-sex
fi

if [ ! -s rsq8-duplicates ]; then
    awk 'NR==FNR{duplicates[$0];next} {f=!($2 in duplicates)} f' rsq8-duplicates "!{rsq8_lz}" >rsq8_lz
    plink --bfile "!{rsq8.bim.baseName}" --memory 10000 --threads 1 --exclude rsq8-duplicates --make-bed --out rsq8 --allow-no-sex
fi

echo "Calculate clumping..."
python -c 'from Stats_helpers import *; \
    PLINK_clumping("rsq3", \
        "rsq8", \
        "!{rsq3_lz}", \
        "!{rsq8_lz}", \
        !{params.clumpr2}, \
        !{params.clumpp1}, \
        !{params.clumpp2}, \
        !{params.clumpkb})'

# rename rsq3_lz "!{rsq3.bim.baseName}" rsq3_lz_clump.clumped*
#mv rsq3_lz_clumped.clumped "!{rsq3_lz}_clump.clumped"
#mv rsq3_lz_clumped.clumped_all "!{rsq3_lz}_clump.clumped_all"
#mv rsq3_lz_clumped.clumped_all.noXMHC "!{rsq3_lz}_clump.clumped_all.noXMHC"
#mv rsq3_lz_clumped.clumped_groups "!{rsq3_lz}_clump.clumped_groups"
#mv rsq3_lz_clumped.clumped_groups.noXMHC "!{rsq3_lz}_clump.clumped_groups.noXMHC"

'''
}

