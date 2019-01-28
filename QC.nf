
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
    file hapmap2_samples

    output:
    set dataset,filebase,batch into Rs_datasets
    set file("${batch}_individuals_annotation_hapmap2.txt"), file("${batch}_individuals_annotation_cases_controls_hapmap2.txt") into SampleQC_hapmap


shell:
'''
cat "!{filebase}_individuals_annotation.txt" "!{hapmap2_samples}" >"!{batch}_individuals_annotation_hapmap2.txt"
cat "!{filebase}_individuals_annotation.txt" "!{hapmap2_samples}" >"!{batch}_individuals_annotation_cases_controls_hapmap2.txt"
'''
}

process Rs {
    tag "$dataset/$batch"
echo true
input:
    set dataset,filebase,batch from Rs_datasets

shell:
    dsconfig = params.dataset_config[dataset]
'''
NXF_PARAMS="!{Rs_script} -c !{params.qc_config} -c !{dsconfig}"
echo "Would run 'nextflow $NXF_PARAMS'"
nextflow $NXF_PARAMS
'''
}

process SNPQCI {
echo true
shell:
'''
echo !{task}
'''
}

process SampleQC {
echo true
    input:
    // bla
    set file(hapmap), file(hapmap_cc) from SampleQC_hapmap
shell:
'''
echo "Would run 'nextflow run !{SampleQC_script} --individuals_annotation_hapmap2 '!{hapmap}''"
'''
}

process SNPQCII {
echo true
shell:
'''
echo !{task}
'''

}

process FinalAnalysis {
echo true
shell:
'''
echo !{task}
'''
}

process Report {
echo true
shell:
'''
echo !{task}
'''
}
