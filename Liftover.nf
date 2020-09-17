
process hg38_liftover {
publishDir "${params.hg38_dir}", mode: 'copy'
tag "${dataset}"
input:
    tuple file(bed), file(bim), file(fam),val(dataset) from Channel.of([file(params.bed), file(params.bim), file(params.fam), params.collection_name])
    
output:
    tuple file("${bed.baseName}_b38.bed"), file("${bim.baseName}_b38.bim"), file("${fam.baseName}_b38.fam"),val(dataset) into split_vcf
    tuple file("${bed.baseName}_b38_noATCG.bed"), file("${bim.baseName}_b38_noATCG.bim"), file("${fam.baseName}_b38_noATCG.fam"),file("${fam.baseName}_b38_noATCG.log")
    file "indels" into split_vcf_indels
    file "atcg" into split_vcf_atcg
    file "unmapped.bed"
shell:
'''
module load Plink/1.9

LIFTOVER=!{params.ucsc_liftover}/liftOver
BASENAME=!{bed.getBaseName()}_b38
CHAIN=!{params.ucsc_liftover}/hg19ToHg38.over.chain.gz

## Translate to chr:pos-pos and chr23 -> chrX
#mawk '{print "chr"$1":"$4"-"$4}' !{bim} \\
#    | sed 's/chr23/chrX/g' \\
#    | sed 's/chr24/chrY/g' \\
#    | sed 's/chr25/chrX/g' \\
#    | sed 's/chr26/chrMT/g' \\
#    > prelift.pos

# Create UCSC BED file from BIM
# Columns: chrX pos-1 pos rsID
<!{bim} mawk '{ print "chr"$1, $4-1, $4, $4, $2 }' \\
    | sed 's/^chr23/chrX/' \\
    | sed 's/^chr24/chrY/' \\
    | sed 's/^chr25/chrX/' \\
    | sed 's/^chr26/chrMT/' >prelift.bed

# lift-over
$LIFTOVER prelift.bed $CHAIN postlift.bed unmapped.bed

# Generate exclude list for unmapped variants
<unmapped.bed mawk '$0~/^#/{next} {print $5}' >unmapped-variants

# Generate update-pos list for Plink
<postlift.bed mawk '{print $5,$3}' >new-pos

/opt/plink2 --bed !{bed} --bim !{bim} --fam !{fam} --exclude unmapped-variants --update-map new-pos --make-pgen --sort-vars  --out ${BASENAME}
/opt/plink2 --pfile ${BASENAME} --make-bed --out ${BASENAME}.tmp

plink --bfile ${BASENAME}.tmp --merge-x no-fail --make-bed --out ${BASENAME}


# collect indels for later removal
mawk 'BEGIN{c["A"]="T";c["C"]="G";c["G"]="C";c["T"]="A"} { if($5==c[$6]) print $2 }' ${BASENAME}.bim >atcg
/opt/plink2 --bfile ${BASENAME} --exclude atcg --make-bed --out ${BASENAME}_noATCG
grep D ${BASENAME}.bim | mawk '{print $2}' | sort | uniq >indels

# determine monomorphic variants
plink --bfile ${BASENAME} --freq --out "base_freq" --allow-no-sex
python -c 'from SNPQC_helpers import *; frq =  Frq(frq_file="base_freq.frq", write_monomorphic_file="monomorphic"); frq.write_monomorphic_variants_file(); del frq'
cat monomorphic >>indels
'''
}

process split_vcf {

//publishDir "${params.hg38_dir}", mode: 'copy'
tag "${dataset}_${chrom}${atcg}"
    input:
    tuple file(bed), file(bim), file(fam), val(dataset) from split_vcf
    file(indels) from split_vcf_indels
    file(atcg_vars) from split_vcf_atcg
    each chrom from Channel.of(1..24)
    each atcg from Channel.of('', '.noATCG')

    output:
    file "*.vcf.gz"
    file "*.vcf.gz.tbi"

shell:
'''
INFIX="!{atcg}"
# ANNOTATION=/work_ifs/ikmb_repository/references/genomes/homo_sapiens/UCSC/hg38/hg38.fa
ANNOTATION=!{params.ucsc_liftover}/hg38.fa

rm -f remove-vars
if [ "$INFIX" = ".noATCG" ]; then
    cat !{indels} !{atcg_vars} | sort | uniq >remove-vars
else
    ln -s !{indels} remove-vars
fi

/opt/plink2 --bed !{bed} --bim !{bim} --fam !{fam} \\
        --exclude remove-vars \\
        --chr !{chrom} --output-chr chrM --export vcf-4.2 --out !{chrom}$INFIX || true

if [ -e "!{chrom}$INFIX.vcf" ]; then
    bgzip <!{chrom}$INFIX.vcf >!{chrom}_tmp.vcf.gz
    tabix !{chrom}_tmp.vcf.gz

    bcftools norm -m -both -N --check-ref s -f $ANNOTATION !{chrom}_tmp.vcf.gz | bgzip >nochmal.vcf.gz
    bcftools norm -m -both -N --check-ref s -f $ANNOTATION nochmal.vcf.gz | bgzip >chr!{chrom}$INFIX.vcf.gz
    tabix chr!{chrom}$INFIX.vcf.gz


    rm -f !{chrom}$INFIX.vcf !{chrom}_tmp.vcf.gz !{chrom}_tmp.vcf.gz.tbi atcg indels
else
    touch NA_!{chrom}.vcf.gz
    touch NA_!{chrom}.vcf.gz.tbi
fi

cp *.vcf.gz !{params.hg38_dir}/
cp *.tbi !{params.hg38_dir}/
'''
}


