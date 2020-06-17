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

// From "QCed" folder, output of QC pipeline
ds_input = ["${params.input}.bed",
         "${params.input}.bim",
         "${params.input}.fam",
         ].collect { fileExists(file(it)) }

covars = fileExists(file(params.covar))

params.min_info_score = 0.3
params.first_chr = 1
params.last_chr = 23
params.disease_data_set_prefix_release_statistics = "dummy"

params.additional_covars = ""


// Converts covariates and phenotype annotations to SAIGE compatible format.
// As case/control status is lost during imputation, its re-added to the Plink dataset.
process preprocess {
    tag "${params.collection_name}"
    time {12.h * task.attempt}

    input:
//    set file(bim), file(bed), file(fam)
    file qced_staged from Channel.from(ds_input).collect()
    file covars
//    file additional_covars

    output:
    //set file(bim), file(bed), file(fam_doubleid) into saige_
    tuple file("qced.doubleid.bim"), file("qced.doubleid.bed"), file("qced.doubleid.fam") into for_prune_qced,for_filter_qced,for_prepare_spa_genotyped
    file pheno into for_fit_pheno
    file covar_cols into for_fit_covar_list

    shell:
    qced = mapFileList(qced_staged)
    fam_doubleid = "qced.doubleid.fam"
    doubleid_basename = "qced.doubleid"
    pheno = "${params.collection_name}.pheno"
    covar_cols = "covar_cols.txt"
'''
echo -n "Generating double ID: !{qced.fam} -> !{fam_doubleid}... "
<"!{qced.fam}" mawk '{$1=$1 "_" $2; $2=$1; print $0}' >"!{fam_doubleid}"
ln -s "!{qced.bed}" "!{doubleid_basename}.bed"
ln -s "!{qced.bim}" "!{doubleid_basename}.bim"
echo "OK"


echo -n "Generating pheno file !{pheno}..."

# Take evec file as a whole, SAIGE does not care about additional columns.
# Filter evec file according to samples contained in FAM (if not already done)

mawk 'FNR==NR{samples[$2];next} {if($2 in samples) {$1=$1 "_" $2; $2=$1; print $0}}' !{qced.fam} "!{covars}"  >"evec.double-id"

EVEC_LINES=$(wc -l <evec.double-id)
FAM_LINES=$(wc -l <!{qced.fam})
if [ $EVEC_LINES -ne $FAM_LINES ]; then
    echo I could not find covariates in !{covars} for every sample specified in !{qced.fam}. Please check.
    echo Aborting.
    exit 1
fi

echo "FID IID PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 batch" >evec.double-id.withheader
cat evec.double-id >>evec.double-id.withheader

# Take phenotype info from FAM, translate to SAIGE-encoding
echo "Pheno" >pheno-column
<"!{fam_doubleid}"  mawk '{if($6=="2") {$6="1";} else if($6=="1") {$6="0";} print $6}' >>pheno-column

# Merge both
paste -d" " evec.double-id.withheader pheno-column >"!{pheno}"

perl -e "print join(',', map {'PC'.\\$_}(1 .. !{params.numofPCs}) )" >covar_cols.txt
'''
}

process prune {
    tag "${params.collection_name}"
    time {24.h * task.attempt}
    input:
    tuple file(bim), file(bed), file(fam) from for_prune_qced

    output:
    tuple file("qced_pruned.bim"), file("qced_pruned.bed"), file("qced_pruned.fam") into for_fit_qced_pruned

    shell:
    '''
    module load IKMB
    module load Plink/1.9
plink --bim "!{bim}" --bed "!{bed}" --fam "!{fam}" --indep-pairwise 50 5 0.2 --out prunedata
plink --bfile "!{bim.baseName}" --extract prunedata.prune.in --make-bed --out qced_pruned
    '''
}

process saige_fit_model {
    tag "${params.collection_name}"
    container "docker://wzhou88/saige:0.36.6"
    time {48.h * task.attempt}
    input:
    tuple file(bim), file(bed), file(fam) from for_fit_qced_pruned
    file covar_list from for_fit_covar_list
    file pheno from for_fit_pheno

    output:
    tuple file("saige_model.rda"), file("saige_model.varianceRatio.txt") into for_saige_spa_imputed_model, for_saige_spa_genotyped_model

shell:
    '''
# for resuming:
rm -f saige_model.varianceRatio.txt
rm -f saige_model.rda

step1_fitNULLGLMM.R \
    --plinkFile=!{bim.baseName} \
    --phenoFile=!{pheno} \
    --phenoCol=Pheno \
    --covarColList=$(cat !{covar_list}) \
    --sampleIDColinphenoFile=IID \
    --traitType=binary        \
    --outputPrefix=saige_model \
    --nThreads=15 \
    --LOCO=FALSE
'''
}

process stage_plink_intermediates {
    tag "${params.collection_name}"
    input:
    file(trace) from Channel.from(file(params.plink_trace))

    output:
    file '*.vcf.gz' into for_filter_vcfgz
    file "*.genotyped.imputed.rsq0.3.chr1-${params.last_chr}.{bed,bim,fam,log}" into for_generate_sumstats_ds, for_clumping_ds
    file "*.genotyped.imputed.rsq0.3.chr1-${params.last_chr}.locuszoom.{bed,bim,fam}" into for_clumping_lz
    file "1-${params.last_chr}.INFO*.vcf.map" into for_generate_sumstats_infomap
    shell:
'''
### Stage INFO-filtered vcf.gz
# extract list of dir hashes from trace
mawk 'FS=";" {if($1=="preprocess_infofilter_dosage") print $4}' <"!{trace}" >hashlist
# complete dir hashes to full directory names
while IFS= read -r line; do compgen -d "!{workflow.workDir}/$line"; done <hashlist >dirlist
# find INFO
xargs -a dirlist -l -I {} -- find {} -name '*INFO*.vcf.gz' >filelist
xargs -a filelist -l -- ln -s

### Stage cleaned BIM/BED/FAM from cleanup_dataset
mawk 'FS=";" {if($1=="cleanup_dataset") print $4}' <"!{trace}" >hashlist
while IFS= read -r line; do compgen -d "!{workflow.workDir}/$line"; done <hashlist >dirlist
# find bim/bed/fam
xargs -a dirlist -l -I {} -- find {} -name '*genotyped.imputed.rsq0.3.chr1-!{params.last_chr}.*' >filelist
xargs -a filelist -l -- ln -s

### Stage VCF info0.3 map for sumstats generation
mawk 'FS=";" {if($1=="preprocess_hrc_vcf_map") print $4}' <"!{trace}" >hashlist
while IFS= read -r line; do compgen -d "!{workflow.workDir}/$line"; done <hashlist >dirlist
# find 1-2.INFO*.vcf.map
xargs -a dirlist -l -I {} -- find {} -name '1-!{params.last_chr}.INFO*.vcf.map' >filelist
xargs -a filelist -l -- ln -s
'''
}

// Extract chromosome number from VCF file name
// 12.vcf.gz -> 12
def get_chrom(vcf) {
    def m = vcf =~ /\/(\d+).INFO.*vcf.gz$/
    return m[0][1]
}

process filter_samples {
    tag "${params.collection_name}.${chrom}"
	time {12.h * task.attempt}
    input:
    tuple val(chrom), file(vcfgz) from for_filter_vcfgz.flatten().map{[ get_chrom(it), it ] }
    tuple file(bim), file(bed), file(fam) from for_filter_qced

    output:
    tuple val(chrom), file(samplefile), file("filtered.vcf.gz"), file("filtered.vcf.gz.tbi") into for_saige_spa_imputed_vcf

    shell:
    targetvcf = "filtered.vcf.gz"
    samplefile="samplefile"
    '''
    mawk '{print $2}' "!{fam}" >keeplist
    vcftools --gzvcf "!{vcfgz}" --keep keeplist --recode --stdout | bgzip >!{targetvcf}
    tabix -p vcf -f !{targetvcf} &

    zcat !{targetvcf} | head -n200 | grep '#CHR' | cut -f 10- | tr -s '[:blank:]' "\\n" >samplefile

    wait
    '''
}

process prepare_spa_genotyped {
    tag "${params.collection_name}.$chrom"
    time {12.h * task.attempt}

    input:
    tuple file(bim), file(bed), file(fam) from for_prepare_spa_genotyped
    val chrom from Channel.of(1..params.last_chr).flatten()

    output:
    tuple val(chrom), file("samplefile"), file("${chrom}.genotyped.vcf.gz"), file("${chrom}.genotyped.vcf.gz.tbi") into for_saige_spa_genotyped_vcf

    shell:
    '''
    module load Plink/1.9
    plink --bim "!{bim}" --bed "!{bed}" --fam "!{fam}" --chr !{chrom} \
        --recode vcf-iid bgz \
        --out !{chrom}.genotyped
    tabix -f !{chrom}.genotyped.vcf.gz
    <!{fam} mawk '{print $2}' >samplefile
'''
}

process saige_spa_imputed {
    tag "${params.collection_name}.$chrom"
    container "docker://wzhou88/saige:0.36.3.2"
    errorStrategy 'retry'
    maxRetries 5
    time {24.h * task.attempt}
    input:
    tuple val(chrom), file(samplefile), file(filteredvcf), file(filteredtbi) from for_saige_spa_imputed_vcf
    tuple file(model), file(varianceRatio) from for_saige_spa_imputed_model

    output:
    file "${params.collection_name}.imputed.SAIGE.chr${chrom}.txt" into for_merge_results_imp

shell:
    '''
step2_SPAtests.R \
    --vcfFile="!{filteredvcf}" \
    --vcfFileIndex="!{filteredtbi}" \
    --vcfField=GT \
    --chrom="!{chrom}" \
    --minMAF=0.001 \
    --minMAC=1 \
    --sampleFile="!{samplefile}" \
    --GMMATmodelFile="!{model}" \
    --varianceRatioFile="!{varianceRatio}" \
    --SAIGEOutputFile=temp.stats \
    --numLinesOutput=2 \
    --IsOutputAFinCaseCtrl=TRUE

# add odds ratio
<temp.stats mawk 'NR==1{print $0 " OR";next} {print $0 " " exp($10)}' \
    >"!{params.collection_name}.imputed.SAIGE.chr!{chrom}.txt"
    '''
}

process saige_spa_genotyped {
    tag "${params.collection_name}.$chrom"
    container "docker://wzhou88/saige:0.36.3.2"
    errorStrategy 'retry'
    maxRetries 5
    time {24.h * task.attempt}
    input:
    tuple val(chrom), file(samplefile), file(vcf), file(tbi) from for_saige_spa_genotyped_vcf
    tuple file(model), file(varianceRatio) from for_saige_spa_genotyped_model

    output:
    file "${params.collection_name}.genotyped.SAIGE.chr${chrom}.txt" into for_merge_results_gen
    shell:
    '''
step2_SPAtests.R \
    --vcfFile="!{vcf}" \
    --vcfFileIndex="!{tbi}" \
    --vcfField=GT \
    --chrom="!{chrom}" \
    --minMAF=0.001 \
    --minMAC=1 \
    --sampleFile="!{samplefile}" \
    --GMMATmodelFile="!{model}" \
    --varianceRatioFile="!{varianceRatio}" \
    --SAIGEOutputFile=temp.stats \
    --numLinesOutput=2 \
    --IsOutputAFinCaseCtrl=TRUE \
    --IsOutputHetHomCountsinCaseCtrl=TRUE \
    --IsOutputNinCaseCtrl=TRUE

# add odds ratio
<temp.stats mawk 'NR==1{print $0 " OR";next} {print $0 " " exp($10)}' \
    >"!{params.collection_name}.genotyped.SAIGE.chr!{chrom}.txt"
    '''
}

process merge_results {
    tag "${params.collection_name}"
    publishDir params.output ?: '.', mode: 'copy', overwrite: true
    time {24.h * task.attempt}

    input:
    file (imp_results) from for_merge_results_imp.collect()
    file (gen_results) from for_merge_results_gen.collect()

    output:
    file "${params.collection_name}.genotyped.SAIGE.chr1-${params.last_chr}.txt"
    file "${params.collection_name}.imputed.SAIGE.chr1-${params.last_chr}.txt"
    file "${params.collection_name}.SAIGE.txt" into for_generate_sumstats_stats

    shell:
    
'''
# Cat everything together keeping only one header line
GENO="!{params.collection_name}.genotyped.SAIGE.chr1-!{params.last_chr}.txt"
cp "!{params.collection_name}.genotyped.SAIGE.chr1.txt" "$GENO"
for ((i=2; i<= !{params.last_chr}; i++)); do
    tail -n +2 "!{params.collection_name}.genotyped.SAIGE.chr${i}.txt" >>"$GENO"
done

IMP="!{params.collection_name}.imputed.SAIGE.chr1-!{params.last_chr}.txt"
cp "!{params.collection_name}.imputed.SAIGE.chr1.txt" "$IMP"
for ((i=2; i<= !{params.last_chr}; i++)); do
    tail -n +2 "!{params.collection_name}.imputed.SAIGE.chr${i}.txt" >>"$IMP"
done

# Now merge. Keep genotyped, add imputed, and sort by chr:pos
mawk 'NR==1{print $0;next} FNR==1{next} FNR==NR{samples[$1":"$2];print $0} {if(!($1":"$2 in samples)) print $0}' "$GENO" "$IMP" \
    | sort -n -k1,2 >"!{params.collection_name}.SAIGE.txt"

'''
}

process generate_sumstats {
publishDir params.output ?: '.', mode: 'copy', overwrite: true
tag "${params.collection_name}"
cpus 2

input:
    file ds_staged from for_generate_sumstats_ds
    file stats from for_generate_sumstats_stats
    file infomap from for_generate_sumstats_infomap

output:
    file "${params.collection_name}.SAIGE.sumstats.gz"
    file "stats-filtered" into for_addL95_U95

shell:
    ds = mapFileList(ds_staged)
'''
# Keep only variants that we have in the cleaned dataset
mawk 'NR==FNR{variants[$1":"$4];next} FNR==1{print $0}{ if($1":"$2 in variants) print $0 }'\
    "!{ds.bim}" "!{stats}" \
    >stats-filtered


saige2sumstats.py "!{ds.log}" stats-filtered >sumstats.tmp
vcf2sumstats.pl "!{infomap}" sumstats.tmp | gzip >"!{params.collection_name}.SAIGE.sumstats.gz"
'''
}

/*
process addL95_U95 {
    tag "${params.collection_name}"
input:
    file stats from for_addL95_U95

output:
    file "${stats}.l95u95" into for_clumping_l95u95, for_clumping_assoc_lz

shell:
pval2qchisq = SCRIPT_DIR + "/pval2qchisq_2.r"
'''

# Extract odds ratio and p value
<"!{stats}" mawk 'NR==1{print "OR P";next}{ print $20, $13 }' >p.or.tmp1

# Make L95 and U95 columns
R --no-save --args p.or.tmp1 p.or.tmp2 <"!{pval2qchisq}"

# Convert space to tab
mawk '{ print $1"\t"$2 }' <p.or.tmp2 >p.or.tmp3

# Attach columns to input file
paste "!{stats}" p.or.tmp3 >"!{stats}.l95u95"

'''
}


process clumping {
tag "${params.collection_name}"
publishDir params.output ?: '.', mode: 'copy', overwrite: true
memory 24.GB
input:
    file rsq3_lz from for_clumping_assoc_lz
    file rsq3_staged from for_clumping_lz

    // python -c 'from Stats_Immunochip import *; PLINK_clumping("!{imp.fam.baseName}", "!{rsq3_lz}", "!{rsq8_lz}")'

output:
    file "${rsq3_lz}_clump.clumped*"


shell:
rsq3 = mapFileList(rsq3_staged)
'''
module load Plink/1.9


# replace runs of tabs and spaces with a single space, take column 2 (sample names) and list duplicates
echo "Detecting duplicates..."
<"!{rsq3.bim}" tr -s '\t ' ' ' | cut -f2 -d' ' | uniq -d >rsq3-duplicates

# output those records that are not in the 'duplicates' list
echo "Checking for and removing duplicates from locuszoom tables..."
if [ ! -s rsq3-duplicates ]; then
    awk 'NR==FNR{duplicates[$0];next} {f=!($2 in duplicates)} f' rsq3-duplicates "!{rsq3_lz}" >rsq3_lz
    plink --bfile "!{rsq3.bim.baseName}" --memory 10000 --threads 1 --exclude rsq3-duplicates --make-bed --out rsq3 --allow-no-sex
fi

echo "Calculate clumping..."
python -c 'from Stats_helpers import *; \
    PLINK_clumping("rsq3", \
        "", \
        "!{rsq3_lz}", \
        "", \
        !{params.clumpr2}, \
        !{params.clumpp1}, \
        !{params.clumpp2}, \
        !{params.clumpkb})'


'''
}
*/
