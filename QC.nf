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
        assertFileExists(file(new File(b).getParent()+"/"+key+"_individuals_annotation.txt"))
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
    -resume"

echo "Would run 'nextflow $NXF_PARAMS'"
nextflow run $NXF_PARAMS
mv trace.txt $MYPWD/Rs-!{dataset}-!{batch}.trace.txt
cd $MYPWD
if [ ! -e !{batch}_Rs.bed ]; then
    ln -s !{params.output}/!{dataset}/Rs/!{batch}_Rs.bed
    ln -s !{params.output}/!{dataset}/Rs/!{batch}_Rs.bim
    ln -s !{params.output}/!{dataset}/Rs/!{batch}_Rs.fam
    ln -s !{params.output}/!{dataset}/Rs/!{batch}.indels
fi

'''
}

process SNPQCI {
    publishDir "!{params.output}/!{dataset}/SNPQCI", mode: 'copy'
    tag "${dataset}"
input:
    set val(dataset), val(filebase), file(bed), file(bim), file(fam), file(indels) from SNPQCI_batches.groupTuple()
output:
    set val(dataset), val(filebase), file("${prefix}.bed"), file("${prefix}.bim"), file("${prefix}.fam") into SampleQC_ds
    
shell:
    individuals_annotation = file(filebase[0]).getParent().toString() + "/" + dataset + "_individuals_annotation.txt"
    prefix = "${dataset}_QCI"
'''
echo Checking !{dataset}:
echo BEDs: "!{bed}"
echo BIMs: "!{bim}"
echo FAMs: "!{fam}"
echo INDELS: "!{indels}"

MYPWD=$(pwd)
mkdir -p !{workflow.workDir}/!{dataset}-SNPQCI
cd !{workflow.workDir}/!{dataset}-SNPQCI

nextflow run !{SNPQCI_script} -c !{params.qc_config} \\
    --rs_dir="$MYPWD" \\
    --snpqci_dir="!{params.output}/!{dataset}/SNPQCI" \\
    --collection_name="!{dataset}" \\
    --chip_defs=!{workflow.projectDir}/config/ChipDefinitions.groovy \\
    --rs_bims="!{bim}" \\
    --rs_beds="!{bed}" \\
    --rs_fams="!{fam}" \\
    --rs_indels="!{indels}" \\
    --individuals_annotation="!{individuals_annotation}" \\
    -resume

mv trace.txt $MYPWD/SNPQCI-!{dataset}.trace.txt
cd $MYPWD

if [ ! -e !{prefix}.bed ]; then
    ln -s !{params.output}/SNPQCI/!{prefix}.bed
    ln -s !{params.output}/SNPQCI/!{prefix}.bim
    ln -s !{params.output}/SNPQCI/!{prefix}.fam
fi
'''
}

process SampleQC {
    tag "${dataset}"
input:
    set val(dataset), val(filebase), file(bed), file(bim), file(fam) from SampleQC_ds
shell:
    individuals_annotation = file(filebase[0]).getParent().toString() + "/" + dataset + "_individuals_annotation.txt"
'''
MYPWD=$(pwd)
mkdir -p !{workflow.workDir}/!{dataset}-SampleQC
cd !{workflow.workDir}/!{dataset}-SampleQC

cat !{individuals_annotation} !{hapmap2_samples} >!{dataset}_individuals_annotation_hapmap2.txt

nextflow run !{SampleQC_script} -c !{params.qc_config} -c !{params.dataset_config[dataset]}\\
    --snpqci_dir="!{params.output}/SNPQCI" \\
    --collection_name="!{dataset}" \\
    --individuals_annotation="!{individuals_annotation}" \\
    --individuals_annotation_hapmap2="$(pwd)/!{dataset}_individuals_annotation_hapmap2.txt" \\
    -resume

mv trace.txt $MYPWD/SampleQC-!{dataset}.trace.txt
'''
}

process SNPQCII {
//echo true
shell:
'''
echo !{task}
'''

}

process FinalAnalysis {
//echo true
shell:
'''
echo !{task}
'''
}

process Report {
//echo true
shell:
'''
echo !{task}
'''
}
