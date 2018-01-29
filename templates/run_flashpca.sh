#!/bin/sh

# Nextflow template.
# Expected variables: prefix dataset[bed] annotation params.preQCIMDS_1kG_sample

module load "IKMB"
module load "Plink/1.7"
module load "FlashPCA"

flashpca2 -d 2000 --bfile "!{dataset.bed.baseName}" \
    --outval !{prefix}_eigenvalues_flashpca2 \
    --outvec !{prefix}_eigenvectors_flashpca2 \
    --outpc  !{prefix}_pcs_flashpca2 \
    --numthreads 15 \
    --outload !{prefix}_loadings_flashpca2 \
    --outmeansd !{prefix}_meansfd_flashpca2 \
    --memory 64000

echo Adding batch info | ts
python -c 'from SampleQCI_helpers import *; addbatchinfo_10PCs("!{prefix}_pcs_flashpca2", "!{prefix}_eigenvalues_flashpca2", "!{prefix}.pca.evec", "!{prefix}.eval", "!{annotation}", "!{params.preQCIMDS_1kG_sample}")'

echo Drawing FLASHPCA2 eigenvectors for set with batch info | ts
R --slave --args "!{dataset.bed.baseName}" "!{params.preQCIMDS_1kG_sample}" <"!{draw_evec_FLASHPCA2}"

echo Adding country info | ts
python -c 'from SampleQCI_helpers import *; addcountryinfo_10PCs("!{prefix}_pcs_flashpca2", "!{prefix}_eigenvalues_flashpca2", "!{prefix}.country.pca.evec", "!{prefix}.country.eval", "!{annotation}", "!{params.preQCIMDS_1kG_sample}")'

echo Drawing FLASHPCA2 eigenvectors for set with country info | ts
R --slave --args "!{dataset.bed.baseName}.country" "!{params.preQCIMDS_1kG_sample}" <"!{draw_evec_FLASHPCA2}"
