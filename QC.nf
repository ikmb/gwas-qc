// vim: syntax=nextflow

/* Verify config parameter*/
params.dataset_prefixes = [:]
if(params.dataset_prefixes.count({true}) == 0) {
    System.err.println("You didn't specify a configuration or the configuration does not contain dataset definitions!")
    return
}

/* Verify dataset file list */
def assertFileExists(f) {
    if( !f.exists() ) {
        System.err.println("Could not find file " + f)
        System.exit(2)
    }
    return true
}

input_datasets = Channel.create()
params.dataset_prefixes.each { key, value ->
    value.each { b -> 
        assertFileExists(file(b + ".bim")) 
        assertFileExists(file(b + ".bed")) 
        assertFileExists(file(b + ".fam")) 
        assertFileExists(file(b + "_individuals_annotation.txt"))
        input_datasets << [key, b, file(b+".bim").getBaseName()]
    }
}
input_datasets.close()

/* Additional external files */
hapmap2_samples = file("/ifs/data/nfs_share/sukmb388/regeneron_annotations/hapmap2-annotations.txt")

Rs_script = file("Rs.nf")
SNPQCI_script = file("SNPQCI.nf")
SampleQC_script = file("SampleQCI.nf")
SNPQCII_script = file("SNPQCII.nf")
FinalAnalysis_script = file("FinalAnalysis.nf")

process PrepareFiles {
    tag "$dataset/$batch"

    input:
    set dataset,filebase,batch from input_datasets

    output:
    set dataset,filebase,batch into Rs_datasets

shell:
'''
echo Dummy
'''
}

process Rs {
    tag "$dataset/$batch"
input:
    set dataset,filebase,batch from Rs_datasets

output:
    file "Rs-${dataset}-${batch}.trace.txt" into trace_reports_rs
    set val(dataset), val(filebase), file("${batch}_Rs.bed"), file("${batch}_Rs.bim"), file("${batch}_Rs.fam"), file("${batch}.indels") into SNPQCI_batches
    //set dataset,batch into SNPQCI_info
    file "Rs-${dataset}-${batch}.trace.txt" into Rs_traces
shell:
    dsconfig = params.dataset_config[dataset]
'''
echo !{task}
echo !{workflow}
MYPWD=$(pwd)
mkdir -p !{workflow.workDir}/!{dataset}-!{batch}-Rs
cd !{workflow.workDir}/!{dataset}-!{batch}-Rs

NXF_PARAMS="!{Rs_script} -c !{params.qc_config} \\
    -c !{dsconfig} \\
    --filebase=!{filebase} \\
    --batch_name=!{batch} \\
    --ds_name=!{dataset} \\
    --chip_defs=!{workflow.projectDir}/config/ChipDefinitions.groovy \\
    --rs_dir=!{params.output}/!{dataset}/Rs \\
    -with-report "!{params.output}/!{dataset}/Rs/!{batch}_execution-report.html" \\
    -resume"

nextflow run $NXF_PARAMS
mv trace.txt Rs-!{dataset}-!{batch}.trace.txt
cp Rs-!{dataset}-!{batch}.trace.txt $MYPWD

cd $MYPWD
    ln -fs !{params.output}/!{dataset}/Rs/!{batch}_Rs.bed
    ln -fs !{params.output}/!{dataset}/Rs/!{batch}_Rs.bim
    ln -fs !{params.output}/!{dataset}/Rs/!{batch}_Rs.fam
    ln -fs !{params.output}/!{dataset}/Rs/!{batch}.indels

'''
}

process SNPQCI {
    tag "${dataset}"
input:
    set val(dataset), val(filebase), file(bed), file(bim), file(fam), file(indels) from SNPQCI_batches.groupTuple()
output:
    set val(dataset), val(filebase), file("${prefix}.bed"), file("${prefix}.bim"), file("${prefix}.fam") into SampleQC_ds
    file "SNPQCI-${dataset}.trace.txt" into SNPQCI_trace
    
shell:
    individuals_annotation = file(filebase[0]).getParent().toString() + "/" + dataset + "_individuals_annotation.txt"
    prefix = "${dataset}_QCI"
'''
echo Checking !{dataset}:
echo BEDs: "!{bed}"
echo BIMs: "!{bim}"
echo FAMs: "!{fam}"
echo INDELS: "!{indels}"
echo Annotation: "!{filebase}_individuals_annotation.txt"

set +x

MYPWD=$(pwd)
mkdir -p !{workflow.workDir}/!{dataset}-SNPQCI
cd !{workflow.workDir}/!{dataset}-SNPQCI

mkdir -p !{params.output}/!{dataset}/SNPQCI

FIRSTLINE=1
BASES=(!{filebase.join(" ")})
for base in "${BASES[@]}"
do
    echo "Processing ${base}_individuals_annotation.txt into !{dataset}_individuals_annotation.txt"
    if [ "$FIRSTLINE" -gt 0 ]; then
        cp ${base}_individuals_annotation.txt !{dataset}_individuals_annotation.txt
        FIRSTLINE=0
    else
        tail -n +2 ${base}_individuals_annotation.txt >>!{dataset}_individuals_annotation.txt
    fi
done

cp !{dataset}_individuals_annotation.txt !{params.output}/!{dataset}/SNPQCI/

nextflow run !{SNPQCI_script} -c !{params.qc_config} \\
    --rs_dir="$MYPWD" \\
    --snpqci_dir="!{params.output}/!{dataset}/SNPQCI" \\
    --collection_name="!{dataset}" \\
    --chip_defs=!{workflow.projectDir}/config/ChipDefinitions.groovy \\
    --rs_bims="!{bim}" \\
    --rs_beds="!{bed}" \\
    --rs_fams="!{fam}" \\
    --rs_indels="!{indels}" \\
    --individuals_annotation="!{params.output}/!{dataset}/SNPQCI/!{dataset}_individuals_annotation.txt" \\
    -with-report "!{params.output}/!{dataset}/SNPQCI/execution-report.html" \\
    -resume

mv trace.txt SNPQCI-!{dataset}.trace.txt
cp SNPQCI-!{dataset}.trace.txt $MYPWD
cd $MYPWD

    ln -fs !{params.output}/!{dataset}/SNPQCI/!{prefix}.bed
    ln -fs !{params.output}/!{dataset}/SNPQCI/!{prefix}.bim
    ln -fs !{params.output}/!{dataset}/SNPQCI/!{prefix}.fam
'''
}

process SampleQC {
    tag "${dataset}"
input:
    set val(dataset), val(filebase), file(bed), file(bim), file(fam) from SampleQC_ds
output:
    set val(dataset), val(filebase), file("${prefix}.bed"), file("${prefix}.bim"), file("${prefix}.fam") into SNPQCII_ds
    file("${prefix}.pca.evec") into FinalAnalysis_evec
    file "SampleQC-${dataset}.trace.txt" into SampleQC_trace
shell:
//    individuals_annotation = file(filebase[0]).getParent().toString() + "/" + dataset + "_individuals_annotation.txt"
    prefix = "${dataset}_SampleQCI_final_withoutRelatives"
'''
MYPWD=$(pwd)
mkdir -p !{workflow.workDir}/!{dataset}-SampleQC
cd !{workflow.workDir}/!{dataset}-SampleQC
ANNO=!{params.output}/!{dataset}/SNPQCI/!{dataset}_individuals_annotation.txt
cat $ANNO !{hapmap2_samples} >!{dataset}_individuals_annotation_hapmap2.txt

nextflow run !{SampleQC_script} -c !{params.qc_config} -c !{params.dataset_config[dataset]}\\
    --snpqci_dir="!{params.output}/!{dataset}/SNPQCI" \\
    --sampleqci_dir="!{params.output}/!{dataset}/SampleQCI" \\
    --collection_name="!{dataset}" \\
    --individuals_annotation="$ANNO" \\
    --individuals_annotation_hapmap2="$(pwd)/!{dataset}_individuals_annotation_hapmap2.txt" \\
    -with-report "!{params.output}/!{dataset}/SampleQCI/execution-report.html" \\
    -resume


mv trace.txt SampleQC-!{dataset}.trace.txt
cp SampleQC-!{dataset}.trace.txt $MYPWD

cd $MYPWD
ln -fs !{params.output}/!{dataset}/SampleQCI/!{prefix}.bed
ln -fs !{params.output}/!{dataset}/SampleQCI/!{prefix}.bim
ln -fs !{params.output}/!{dataset}/SampleQCI/!{prefix}.fam
ln -fs !{params.output}/!{dataset}/SampleQCI/!{prefix}.pca.evec
'''
}

process SNPQCII {
    tag "${dataset}"
input:
    set val(dataset), val(filebase), file(bed), file(bim), file(fam) from SNPQCII_ds
output:
    set val(dataset), val(filebase), file("${prefix}.bed"), file("${prefix}.bim"), file("${prefix}.fam"), file("${prefix}_annotation.txt") into FinalAnalysis_ds
    file "SNPQCII-${dataset}.trace.txt" into SNPQCII_trace
shell:
    prefix = "${dataset}_QCed"
//    individuals_annotation = file(filebase[0]).getParent().toString() + "/" + dataset + "_individuals_annotation.txt"
'''
MYPWD=$(pwd)
mkdir -p !{workflow.workDir}/!{dataset}-SNPQCII
cd !{workflow.workDir}/!{dataset}-SNPQCII
ANNO=!{params.output}/!{dataset}/SNPQCI/!{dataset}_individuals_annotation.txt

nextflow run !{SNPQCII_script} -c !{params.qc_config} -c !{params.dataset_config[dataset]} \\
    --sampleqci_dir="!{params.output}/!{dataset}/SampleQCI" \\
    --qc_dir="!{params.output}/!{dataset}/SNPQCII" \\
    --collection_name="!{dataset}" \\
    --individuals_annotation="$ANNO" \\
    -with-report "!{params.output}/!{dataset}/SNPQCII/execution-report.html" \\
    -resume

mv trace.txt SNPQCII-!{dataset}.trace.txt
cp SNPQCII-!{dataset}.trace.txt $MYPWD

cd $MYPWD
ln -fs !{params.output}/!{dataset}/SNPQCII/!{prefix}.bed
ln -fs !{params.output}/!{dataset}/SNPQCII/!{prefix}.bim
ln -fs !{params.output}/!{dataset}/SNPQCII/!{prefix}.fam
ln -fs !{params.output}/!{dataset}/SNPQCII/!{prefix}_annotation.txt
'''

}

process FinalAnalysis {
tag "${dataset}"
input:
    set val(dataset), val(filebase), file(bed), file(bim), file(fam), file(anno) from FinalAnalysis_ds
    file(evec) from FinalAnalysis_evec
output:
    file "FinalAnalysis-${dataset}.trace.txt" into FinalAnalysis_trace
    set val(dataset), val(filebase) into ReportData
shell:
'''
MYPWD=$(pwd)
mkdir -p !{workflow.workDir}/!{dataset}-FinalAnalysis
cd !{workflow.workDir}/!{dataset}-FinalAnalysis

nextflow run !{FinalAnalysis_script} -c !{params.qc_config} -c !{params.dataset_config[dataset]} \\
    --collection_name="!{dataset}" \\
    --individuals_annotation="$MYPWD/!{anno}" \\
    --bed="$MYPWD/!{bed}" \\
    --bim="$MYPWD/!{bim}" \\
    --fam="$MYPWD/!{fam}" \\
    --evec="$MYPWD/!{evec}" \\
    --qc_dir="!{params.output}/!{dataset}/QCed" \\
    -with-report "!{params.output}/!{dataset}/QCed/execution-report.html" \\
    -resume

mv trace.txt FinalAnalysis-!{dataset}.trace.txt
cp FinalAnalysis-!{dataset}.trace.txt $MYPWD
ln -fs !{params.output}/!{dataset}/SNPQCII/!{dataset}_QCed.bed !{params.output}/!{dataset}/QCed/
ln -fs !{params.output}/!{dataset}/SNPQCII/!{dataset}_QCed.bim !{params.output}/!{dataset}/QCed/
ln -fs !{params.output}/!{dataset}/SNPQCII/!{dataset}_QCed.fam !{params.output}/!{dataset}/QCed/
ln -fs !{params.output}/!{dataset}/SNPQCII/!{dataset}_QCed.log !{params.output}/!{dataset}/QCed/
'''
}

process Report {
tag "${dataset}"
echo true
input:
    set val(dataset), val(filebase) from ReportData
    file "Rs-trace" from Rs_traces
    file "SNPQCI-trace" from SNPQCI_trace
    file "SampleQC-trace" from SampleQC_trace
    file "SNPQCII-trace" from SNPQCII_trace
    file "Final-trace" from FinalAnalysis_trace
shell:
'''
echo "$(ls)"
'''
}
