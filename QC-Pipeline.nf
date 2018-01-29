// -*- mode:java -*-

/* 
   Processes the configuration and feeds the batches into the first 
   pipeline stage
*/

// Target folders for different stages
params.Rs_out = "./out-Rs"

Rs_out = file(params.Rs_out)

process run_rs {
    echo true
    executor "local"
    tag "$batch_prefix"
    
    input:
    val batch_prefix from Channel.from(params.disease_names.split(","))

    shell:
    input_set = BATCH_DIR + "/$batch_prefix/orig_files/$batch_prefix"
'''
    mkdir -p "!{params.Rs_out}"
    nextflow -c "!{NXF_DIR}/Rs.nf" -c "!{batch_prefix}/!{batch_prefix}.config" \
             run -profile=subpipeline "!{NXF_DIR}/Rs.nf" \
             --input="!{input_set}" \
             --output="!{Rs_out}" \
             --resume
'''
}
/*
process run_SNPQCI {
    input:
    file rs from for_SNPQCI

    output:
    file "*${params.collection_name}_QCI.{bed,bim,fam}" into for_SampleQCI

    module "IKMB"
    module "Nextflow"

"""
mkdir -p "${snpqci_out}"
nextflow -c ${snpqci_config} -c ${set_config} run ${snpqci_nxf} \
         --input=${params.input} --output=${snpqci_out} \
         -resume
"""
}

process run_SampleQCI {
    input:
    file dataset from for_SampleQCI

    output:
    file "*${params.collection_name}_SampleQCI.{bed,bim,fam}" for_SNPQCII
    file "*${params.collection_name}_SampleQCI_withoutRelatives_pruned.{bed,bim,fam}" into for_SNPQCII_withoutRelatives
    file "${params.collection_name}_SampleQCI_withoutRelatives.pca.evec" into for_SNPQCII_withoutRelatives_pc

    module "IKMB"
    module "Nextflow"

"""
mkdir -p "${sampleqci_out}"
nextflow -c ${sampleqci_config} -c ${set_config} run ${sampleqci_nxf} \
         --input=${dataset[0].baseName} --output=${sampleqci_out} \
         -resume
"""
}

process run_SNPQCII {
    input:
    file dataset from for_SNPQCII
    file dataset_withoutRelatives from for_SNPQCII_withoutRelatives
    file pcfile from for_SNPQCII_withoutRelatives_pc

    output:
    // ?

    module "IKMB"
    module "Nextflow"

"""
mkdir -p "${snpqcii_out}"
nextflow -c "${snpqcii_config}" -c ${set_config} run ${snpqcii_nxf} \
         --input=${dataset[0].baseName} --input_wr=${dataset_withoutRelatives[0].baseName} \
         --output=${snpqcii_out} \
         -resume
"""
}
*/
