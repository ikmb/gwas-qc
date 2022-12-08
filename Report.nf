// params.traces - filename with list of traces, in order for report script
// params.system - system info file

in_system = file(params.system)
in_traces = file(params.traces)
in_dataset = params.dataset


process Report {
label "smallmem"
label "shortrunning"

tag "${in_dataset}"
publishDir ".", mode: 'copy'
input:
    file in_system
    file in_traces
output:
    file "${in_dataset}-report.pdf"
shell:
    report_dir = "${workflow.projectDir}/report"
'''

export NXF_WORK=!{workflow.workDir}
export WARN_SNPQC=!{params.skip_snpqc}
export WARN_SAMPLEQC=!{params.skip_sampleqc}
export WARN_RELATED=!{params.keep_related}
export PERL5LIB=!{report_dir}
export HAPMAP_PCA_ENABLED=!{params.activate_hapmap_pca}
echo "Plink 1.9;$(plink --version)" >>container.txt
echo "Plink 2;$(plink2 --version)" >>container.txt
echo "FlashPCA;$(flashpca2 --version 2>&1 | head -n2 | tail -n1)" >>container.txt
echo "bcftools;$(bcftools version | head -n1)" >>container.txt
echo "Eigensoft;6.1.4" >>container.txt
echo "R;$(R --version | head -n1)" >>container.txt
echo "perl;$(perl -e 'print $^V')" >>container.txt
echo "python;$(python -V 2>&1)" >>container.txt

cat !{in_system} container.txt >lastpage.txt

LINE=$(tr '\\n' ' ' <!{in_traces})

perl !{report_dir}/report.pl \\
    !{workflow.workDir} \\
    !{report_dir}/preamble.tex \\
    $LINE \\
    lastpage.txt

latexmk -lualatex report
mv report.pdf "!{in_dataset}-report.pdf"

'''
}

