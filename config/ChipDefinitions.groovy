class ChipDefinitions {

    static String Producer(String prod) {
        def chip_producer_allowed = ["Illumina" : "Illumina", "Affymetrix" : "Affymetrix"]
        return chip_producer_allowed[prod];
    }

    static String SNPAnnotations(String ver) {
        def chip_versions_allowed = [
            "Illu300v3" : "unknown",
            "Illu550"   : "unknown",
            "Immunochip"    : "ichip.hg18.hg19.dbsnpID.chr1-26.txt",
            "Exomechipv1"     : "v1_exomearray.hg19.dbsnpID.chr1-26.txt",
            "Exomechipv1-1"   : "v1-1_exomearray.hg19.dbsnpID.chr1-26.txt",
            "HumanCoreExome24v1" : "HumanCoreExome24_v1.0.hg19.dbsnpID.chr1-26.txt",
            "AgiFranceCustom" : "Agifrance_custom.annot_nofilter.txt",
            "Affy6" : "GenomeWideSNP_6.na24.annot_nofilter.txt",
            "Affy5" : "GenomeWideSNP_5.na24.annot_nofilter.txt",
            "Affy500kSet" : "Mapping250K.na24.annot_nofilter.txt",
            "GSAarrayv1" : "GSAarrayv1.hg19.chr1-26.txt"
        ]

        return chip_versions_allowed[ver];
    }

    static String StrandInfo(String ver) {
        def chip_strand_info_allowed = [
            "Immunochip_orig_annotation":"ichip.orig_annotation.hg18.hg19.dbsnpID.chr1-26.minusStrandOnly.rs.hg19.txt",
            "Immunochip_TOP_annotation":"ichip.TOP_annotation.hg18.hg19.dbsnpID.chr1-26.FlipToPlusStrandOnly.rs.hg19.txt",
            "Exomechipv1_orig_annotation":"v1_exomearray.hg19.dbsnpID.chr1-26.minusStrandOnly.rs.txt",
            "Exomechipv1-1_orig_annotation":"v1-1_exomearray.hg19.dbsnpID.chr1-26.minusStrandOnly.rs.txt",
            "HumanCoreExome24v1_orig_annotation":"HumanCoreExome24_v1.0.hg19.dbsnpID.chr1-26.minusStrandOnly.rs.txt",
            "GSAarrayv1_orig_annotation" : "GSAarrayv1.hg19.chr1-26.minusStrandOnly.rs.txt"
        ]

        return chip_strand_info_allowed[ver];
    }

    static String RsExclude(String ver) {
        def chip_rs_exclude = [
            "Immunochip":"ichip.hg18.hg19.dbsnpID.chr1-26.exclude.HailiangHuang.chr25.chr26.txt",
            "HumanCoreExome24v1":"HumanCoreExome24_v1.0.hg19.dbsnpID.chr1-26.duplicates.txt",
            "GSAarrayv1" : "GSAarrayv1.hg19.chr1-26.exclude.txt"
        ]

        return chip_rs_exclude[ver];
    }


    static String RsAutosomes(String ver) {
        def chip_rs_autosomes  = [
            "Immunochip" : "ichip.hg19.dbsnpID.chr1-22.rs_id.txt",
            "HumanCoreExome24v1" : "HumanCoreExome24_v1.0.hg19.dbsnpID.chr1-22.rs_id.txt",
	    "GSAarrayv1" : "GSAarrayv1.hg19.snp147Common.dbsnpID.chr1-22.rs_id.txt"
        ]
        return chip_rs_autosomes[ver];
    }
}

