// vim: syntax=nextflow

params.nxfdir = "."

/* Verify config parameter*/
params.imputed_dataset = [:]
if(params.imputed_dataset.count({true}) == 0) {
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

Plink_script = file("${params.nxfdir}/Stats_Immunochip.nf")
Saige_script = file("${params.nxfdir}/Stats_SAIGE.nf")

input_datasets_plink = Channel.create()
input_datasets_saige = Channel.create()
params.imputed_dataset.each { key, value ->
        input_datasets_plink << [key, value, params.qced_dataset[key], params.covars[key], params.dataset_config[key]]
        input_datasets_saige << [key, value, params.qced_dataset[key], params.covars[key], params.dataset_config[key]]
}
//input_datasets.dump(tag: 'debug')

process PlinkAssoc {
	publishDir "${params.output}/${ds_name}/Plink_Assoc", mode: 'copy'
    tag "$ds_name"
    input:
    set ds_name, impdir, qced, covars, config from input_datasets_plink

    output:
    file "trace.txt" into saige_trace
    file "*.sumstats.gz"
    file "*.jpg"
    file "*.bim"
    file "*.bed"
    file "*.fam"
    file "dag.html"
    
shell:
'''
NXF_PARAMS="-c !{config} !{Plink_script} --input=!{qced} --output="." --input_imp=!{impdir} --covar=!{covars} -resume -ansi-log false -with-dag dag.html --collection_name=!{ds_name}" 
echo params: $NXF_PARAMS
nextflow run $NXF_PARAMS
'''
}

process SAIGEAssoc {
    publishDir "${params.output}/${ds_name}/SAIGE_Assoc", mode: 'copy'
    tag "$ds_name"

    input:
    set ds_name, impdir, qced, covars, config from input_datasets_saige
    file "plink_trace.txt" from saige_trace

    output:
    stdout into for_fertig
    file "${ds_name}.imputed.SAIGE.chr1-22.txt"
    file "${ds_name}.genotyped.SAIGE.chr1-22.txt"
    file "${ds_name}.SAIGE.txt"
    file "trace.txt"
    file "dag.html"
    file "${ds_name}.SAIGE.sumstats.gz"
    file "*clump*"
shell:
'''
NXF_PARAMS="-c !{config} !{Saige_script} --input=!{qced} --output="." --plink_trace=plink_trace.txt --covar=!{covars} -resume -ansi-log false -with-dag dag.html --collection_name=!{ds_name}"
echo params: $NXF_PARAMS
nextflow run $NXF_PARAMS
'''
}

process Fertig {
input:
val dummy from for_fertig
shell:
'''

'''
}
