#!/bin/sh

# Nextflow template.
# Expected variables: dataset[log,bed] pca[eval] params.twtable

NUM_CASES=$(grep -P 'After.*\d+.cases' !{dataset.log} | cut -d' ' -f3)
NUM_CONTROLS=$(grep -P 'After .*\d+.controls' !{dataset.log} | cut -d' ' -f5)
echo Cases: $NUM_CASES, controls: $NUM_CONTROLS
NUM_SAMPLES=$(($NUM_CASES + $NUM_CONTROLS))
NUM_SNPS=$(grep -P 'After.*\d+.SNPs' !{dataset.log} | cut -d' ' -f8)
echo Samples: $NUM_SAMPLES, markers: $NUM_SNPS

head -n 2000 !{pca.eval} >eval.tw

if [ "!{params.projection_on_populations_CON_only}" == "True" ]; then
    if [ "$NUM_CONTROLS" -gt "$NUM_SNPS" ]; then
        let NUM_CONTROLS--
        twstats -t !{params.twtable} -i eval.tw -o !{dataset.bed.baseName}.eval.tracy_widom -n "$NUM_CONTROLS"
    else
        twstats -t !{params.twtable} -i eval.tw -o !{dataset.bed.baseName}.eval.tracy_widom
    fi
else
    if [ "$NUM_SAMPLES" -gt "$NUM_SNPS" ]; then
        let NUM_SAMPLES--
        twstats -t !{params.twtable} -i eval.tw -o !{dataset.bed.baseName}.eval.tracy_widom -n "$NUM_SAMPLES"
    else
        twstats -t !{params.twtable} -i eval.tw -o !{dataset.bed.baseName}.eval.tracy_widom
    fi
fi
