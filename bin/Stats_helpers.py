#!/usr/bin/env python

import sys
import re
import os

from os.path import *
#import string

# import gzip
# import math
# import decimal
# import datetime
# from os import listdir
# import subprocess

# may also need some of these:

# import Ingos lib
#sys.path.append(os.path.join(os.path.dirname[0], "../../all_scripts"))
sys.path.append(os.environ['PYLIB_DIR'] + "/all_scripts")
sys.path.append(os.environ['PYLIB_DIR'] + "/lib")
from all_common import Command

# import my lib
# sys.path.append(join(sys.path[0], "../lib"))
# sys.path.append(os.environ['PYLIB_DIR'] + "/lib")

#from plink_classes import PackedPed
from eigenstrat_classes import PackedPed


def extractQQplot_null_variants(assoc_input, assoc_output):
    """ extract null variants from assoc file """

    try:
        fh_r = file(QQplot_null_variants, "r")
    except IOError, e:
      print e
      sys.exit(1)
    
    null = {}
    line = fh_r.readline().rstrip('\n')
    while line:
        list = re.split("\s+",line)
        rs   = list[0]
        null[rs] = True
        line = fh_r.readline().rstrip('\n')
    
    fh_r.close()
    
    try:
        fh_r = file(assoc_input, "r")
        fh_w = file(assoc_output, "w")
    except IOError, e:
      print e
      sys.exit(1)
    
    line = fh_r.readline().rstrip('\n')
    fh_w.writelines(line + "\n")
    line = fh_r.readline().rstrip('\n')
    while line:
        list = re.split("\s+",line)
        if list[0] == "":
            del list[0]
        rs = list[1]
        if null.has_key(rs):
            fh_w.writelines(line + "\n")
        line = fh_r.readline().rstrip('\n')
    
    fh_r.close()
    fh_w.close()


def extract_Rsq_variants(assoc_logistic_input, assoc_dosage_input, assoc_merge_output_rsq0_4, assoc_merge_output_rsq0_8):
    """ extract imputed variants based on Rsq value """

    try:
        fh_w1 = file(assoc_merge_output_rsq0_8, "w")
        fh_w2 = file(assoc_merge_output_rsq0_4, "w")
    except IOError, e:
        print e
        sys.exit(1)

    store_lines_rsq0_4 = []
    store_lines_rsq0_8 = []
    
    snp2line = {}

    # -- assoc_logistic_input - all genotyped marker -- #
    try:
        fh_r = file(assoc_logistic_input, "r")
    except IOError, e:
        print e
        sys.exit(1)
    # header
    line = fh_r.readline().rstrip('\n')
    line = fh_r.readline().rstrip('\n')
    while line:
        list = re.split("\s+",line)
        if list[0] == "":
            del list[0]
        if list[-1] == "":
            del list[-1]

        # rs numbers instead of chr:pos in column SNP
        snp2line[list[0]+":"+list[2]] = line
        chr  = decimal.Decimal(list[0])
        pos  = decimal.Decimal(list[2])
        store_lines_rsq0_4.append( (chr, pos, line) )
        store_lines_rsq0_8.append( (chr, pos, line) )
        line = fh_r.readline().rstrip('\n')
    fh_r.close()

    # -- assoc_dosage_input - all imputed marker -- #
    try:
        fh_r = file(assoc_dosage_input, "r")
    except IOError, e:
        print e
        sys.exit(1)

    # header
    header = fh_r.readline().rstrip('\n')
    fh_w1.writelines(header+"\n")
    fh_w2.writelines(header+"\n")
    line = fh_r.readline().rstrip('\n')
    while line:
        list = re.split("\s+",line)
        if list[0] == "":
            del list[0]
        if list[-1] == "":
            del list[-1]
        # if not already genotyped --> add
        if not snp2line.has_key(list[0]+":"+list[2]):
            if float(list[7]) >= 0.4:
                chr  = decimal.Decimal(list[0])
                pos  = decimal.Decimal(list[2])
                store_lines_rsq0_4.append( (chr, pos, line) )
                if float(list[7]) >= 0.8:
                    store_lines_rsq0_8.append( (chr, pos, line) )
        line = fh_r.readline().rstrip('\n')
    fh_r.close()

    # sort merged files by chr, pos
    s_rsq0_4 = sorted(store_lines_rsq0_4, key=lambda tupel: tupel[1])
    t_rsq0_4 = sorted(s_rsq0_4, key=lambda tupel: tupel[0], reverse=False)
    s_rsq0_8 = sorted(store_lines_rsq0_8, key=lambda tupel: tupel[1])
    t_rsq0_8 = sorted(s_rsq0_8, key=lambda tupel: tupel[0], reverse=False)

    for tupel in t_rsq0_4:
        line = tupel[2]
        fh_w2.writelines(line+"\n")

    for tupel in t_rsq0_8:
        line = tupel[2]
        fh_w1.writelines(line+"\n")

    fh_w1.close()
    fh_w2.close()

    # save another version with chr:pos as SNPids for Locuszoom
    try:
        fh_w1 = file(assoc_merge_output_rsq0_8+".locuszoom", "w")
        fh_w2 = file(assoc_merge_output_rsq0_4+".locuszoom", "w")
    except IOError, e:
        print e
        sys.exit(1)
    
    fh_w1.writelines(header+"\n")
    fh_w2.writelines(header+"\n")
    
    check_duplicates_rsq0_4 = {} 
    check_duplicates_rsq0_8 = {} 

    for tupel in t_rsq0_4:
        line = tupel[2]
        list = re.split("\s+",line)
        if not check_duplicates_rsq0_4.has_key(list[0]+":"+list[2]):
            fh_w2.writelines(list[0]+"\t"+\
                         list[0]+":"+list[2] +"\t"+\
                         list[2]+"\t"+\
                         list[3]+"\t"+\
                         list[4]+"\t"+\
                         list[5]+"\t"+\
                         list[6]+"\t"+\
                         list[7]+"\t"+\
                         list[8]+"\t"+\
                         list[9]+"\t"+\
                         list[10]+"\n")
            check_duplicates_rsq0_4[list[0]+":"+list[2]] = True

    for tupel in t_rsq0_8:
        line = tupel[2]
        list = re.split("\s+",line)
        if not check_duplicates_rsq0_8.has_key(list[0]+":"+list[2]):
            fh_w1.writelines(list[0]+"\t"+\
                         list[0]+":"+list[2] +"\t"+\
                         list[2]+"\t"+\
                         list[3]+"\t"+\
                         list[4]+"\t"+\
                         list[5]+"\t"+\
                         list[6]+"\t"+\
                         list[7]+"\t"+\
                         list[8]+"\t"+\
                         list[9]+"\t"+\
                         list[10]+"\n")
            check_duplicates_rsq0_8[list[0]+":"+list[2]] = True

    fh_w1.close()
    fh_w2.close()


def excludeQQplot_variants(assoc_input, assoc_output):
    """ exclude specific regions from assoc null file """

    try:
        fh_r = file(QQplot_SNPexcludeList, "r")
    except IOError, e:
      print e
      sys.exit(1)
    
    exclude_regions = []
    line = fh_r.readline().rstrip('\n')
    while line:
        list = re.split("\s+",line)
        chr    = int(list[0])
        left   = int(list[1])
        right  = int(list[2])
        exclude_regions.append((chr, left, right))
        line = fh_r.readline().rstrip('\n')
    
    fh_r.close()
    
    try:
        fh_r = file(assoc_input, "r")
        fh_w = file(assoc_output, "w")
    except IOError, e:
      print e
      sys.exit(1)
    
    line = fh_r.readline().rstrip('\n')
    fh_w.writelines(line + "\n")
    line = fh_r.readline().rstrip('\n')
    while line:
        list = re.split("\s+",line)
        if list[0] == "":
            del list[0]
        chr = int(list[0])
        pos = int(list[2])
        exclude_this_region = False
        for region in exclude_regions:
            if chr == region[0] and pos >= region[1] and pos <= region[2]:
                exclude_this_region = True
                break
        if exclude_this_region:
            pass
        else:
            fh_w.writelines(line + "\n")
        line = fh_r.readline().rstrip('\n')
    
    fh_r.close()
    fh_w.close()


def qqplot_xMHC_noxMHC(file_PLINK, file_PLINK_noxMHC, numof_cases, numof_controls, qqplotp, qqplotp2, pval2chisq, snpexclude):
    """ extract variants based on Rsq value """

    # ---------------------------- #
    # -- filtered for non-xMHC  -- #
    # -- xMHC: chr6, [25,34[ Mb -- #
    # ---------------------------- #

    print "\n    qq-plot from association results ...\n\n"
    os.system("R --slave --args %s %s %s %s < %s" \
        %(file_PLINK_noxMHC,\
        str(numof_cases),\
        str(numof_controls),\
        collection_name,\
        pval2chisq))
    
    print "\n    another qq-plotpval from association results ...\n\n"
    os.system("R --slave --args %s %s %s %s < %s" \
        %(file_PLINK_noxMHC,\
        str(numof_cases),\
        str(numof_controls),\
        collection_name,\
        qqplotp))
    
    print "\n    another qq-plotpval from association results ...\n\n"
    os.system("R --slave --args %s %s %s %s < %s" \
        %(file_PLINK_noxMHC,\
        str(numof_cases),\
        str(numof_controls),\
        collection_name,\
        qqplotp2))
    
    # -------------------------------------- #
    # -- including xMHC: chr6, [25,34[ Mb -- #
    # -------------------------------------- #
     
    print "\n    qq-plot from association results ...\n\n"
    os.system("R --slave --args %s %s %s %s < %s" \
        %(file_PLINK,\
        str(numof_cases),\
        str(numof_controls),\
        collection_name,\
        pval2chisq))
    
    print "\n    another qq-plotpval from association results ...\n\n"
    os.system("R --slave --args %s %s %s %s < %s" \
        %(file_PLINK,\
        str(numof_cases),\
        str(numof_controls),\
        collection_name,\
        qqplotp))
    
    print "\n    another qq-plotpval from association results ...\n\n"
    os.system("R --slave --args %s %s %s %s < %s" \
        %(file_PLINK,\
        str(numof_cases),\
        str(numof_controls),\
        collection_name,\
        qqplotp2))

    if QQplot_SNPexcludeList_exists:
    
        excludeQQplot_variants(assoc_input=file_PLINK, \
                               assoc_output=file_PLINK + "_excludeRegions")
    
        excludeQQplot_variants(assoc_input=file_PLINK_noxMHC, \
                               assoc_output=file_PLINK_noxMHC + "_excludeRegions")
    
        # ---------------------------- #
        # -- filtered for non-xMHC  -- #
        # -- xMHC: chr6, [25,34[ Mb -- #
        # ---------------------------- #

        print "\n    qq-plot from association results ...\n\n"
        os.system("R --slave --args %s %s %s %s < %s" \
            %(file_PLINK_noxMHC + "_excludeRegions",\
            str(numof_cases),\
            str(numof_controls),\
            collection_name,\
            pval2chisq))
        
        print "\n    another qq-plotpval from association results ...\n\n"
        os.system("R --slave --args %s %s %s %s < %s" \
            %(file_PLINK_noxMHC + "_excludeRegions",\
            str(numof_cases),\
            str(numof_controls),\
            collection_name,\
            qqplotp))
        
        print "\n    another qq-plotpval from association results ...\n\n"
        os.system("R --slave --args %s %s %s %s < %s" \
            %(file_PLINK_noxMHC + "_excludeRegions",\
            str(numof_cases),\
            str(numof_controls),\
            collection_name,\
            qqplotp2))
        
        os.system("rm -f %s" %(file_PLINK_noxMHC + "_excludeRegions"))
        
        # -------------------------------------- #
        # -- including xMHC: chr6, [25,34[ Mb -- #
        # -------------------------------------- #
         
        print "\n    qq-plot from association results ...\n\n"
        os.system("R --slave --args %s %s %s %s < %s" \
            %(file_PLINK + "_excludeRegions",\
            str(numof_cases),\
            str(numof_controls),\
            collection_name,\
            pval2chisq))
        
        print "\n    another qq-plotpval from association results ...\n\n"
        os.system("R --slave --args %s %s %s %s < %s" \
            %(file_PLINK + "_excludeRegions",\
            str(numof_cases),\
            str(numof_controls),\
            collection_name,\
            qqplotp))
        
        print "\n    another qq-plotpval from association results ...\n\n"
        os.system("R --slave --args %s %s %s %s < %s" \
            %(file_PLINK + "_excludeRegions",\
            str(numof_cases),\
            str(numof_controls),\
            collection_name,\
            qqplotp2))



def qqplot_null(file_PLINK_null, numof_cases, numof_controls, qqplotp, qqplotp2, pval2chisq, snpexclude):
    """ extract null variants """

    # use normal qqplot scripts instead of null qqplot scripts because usage of
    # ld pruned variants as null SNPs
    print "\n    qq-plotpval for null SNVs from association results ...\n\n"
    os.system("R --no-save --args %s %s %s %s < %s" \
        %(file_PLINK_null,\
        str(numof_cases),\
        str(numof_controls),\
        collection_name,\
        qqplotp))

    print "\n    qq-plotpval for null SNVs from association results ...\n\n"
    os.system("R --no-save --args %s %s %s %s < %s" \
        %(file_PLINK_null,\
        str(numof_cases),\
        str(numof_controls),\
        collection_name,\
        qqplotp2))
    
    print "\n    qq-plot for null SNVs from association results ...\n\n"
    os.system("R --no-save --args %s %s %s %s < %s" \
        %(file_PLINK_null,\
        str(numof_cases),\
        str(numof_controls),\
        collection_name,\
        pval2chisq))
    
    if snpexclude:
    
        excludeQQplot_variants(assoc_input=file_PLINK_null, \
                               assoc_output=file_PLINK_null + "_excludeRegions")
    
        print "\n    qq-plotpval for null SNVs with specific regions excluded from association results ...\n\n"
        os.system("R --no-save --args %s %s %s %s < %s" \
            %(file_PLINK_null + "_excludeRegions",\
            str(numof_cases),\
            str(numof_controls),\
            collection_name,\
            qqplotp))

        os.system("R --no-save --args %s %s %s %s < %s" \
            %(file_PLINK_null + "_excludeRegions",\
            str(numof_cases),\
            str(numof_controls),\
            collection_name,\
            qqplot2))


def generateHitSpecfile(chromosome, pos, snp, GWAS_regions, file_HitSpec):
    """ """ 
    try:
        fh_w = open(file_HitSpec, "w")
    except IOError, e:
        print e
        sys.exit(1)
    
    sep = "\t"
    fh_w.writelines("snp" +sep+\
          "chr" +sep+\
          "start" +sep+\
          "end" +sep+\
          "flank" +sep+\
          "run" +sep+\
          "m2zargs\n")
    
    flank        = "500kb"
    run          = "yes"
    start        = "NA" 
    end          = "NA" 
    #phenos       = list[5]
    
    m2zargs      = "title=\""+ basename(disease_data_set_prefix_release_statistics) +" "+ snp +"\""
    
    # -- highlight known GWAS region if in plotting region -- #
    if QQplot_SNPexcludeList_exists:
        left  = int(pos) - 500000
        right = int(pos) + 500000
        for region in GWAS_regions:
            if chromosome == str(region[0]):
                # whole plotting region is known locu
                if int(region[1]) <= left and right <= int(region[2]):
                    hiStart      = str(left)
                    hiEnd        = str(right) 
                    m2zargs      += " hiStart=\""+ hiStart + "\"" + " hiEnd=\""+ hiEnd + "\" hiColor=blue"
                # known locus intersects left plotting border
                elif int(region[1]) <= left and (left < int(region[2]) and int(region[2]) <= right):
                    hiStart      = str(left)
                    hiEnd        = str(region[2]) 
                    m2zargs      += " hiStart=\""+ hiStart + "\"" + " hiEnd=\""+ hiEnd + "\" hiColor=blue"
                # known locus intersects right plotting border
                elif (left <= int(region[1]) and int(region[1]) < right) and right <= int(region[2]):
                    hiStart      = str(region[1])
                    hiEnd        = str(right) 
                    m2zargs      += " hiStart=\""+ hiStart + "\"" + " hiEnd=\""+ hiEnd + "\" hiColor=blue"
                # known locus within left and right plotting border
                elif (left <= int(region[1]) and int(region[1]) < right) and (left < int(region[2]) and int(region[2]) <= right):
                    hiStart      = str(region[1])
                    hiEnd        = str(region[2]) 
                    m2zargs      += " hiStart=\""+ hiStart + "\"" + " hiEnd=\""+ hiEnd + "\" hiColor=blue"
    
    m2zargs      += " showAnnot=T showRefsnpAnnot=T annotPch=\"21,24,24,25,22,22,8,7\""
    m2zargs      += " rfrows=4" 
    m2zargs      += " refsnpTextSize=0.8" 
    m2zargs      += " legendSize=0.4" 
    #m2zargs      += " snpset=\"HapMap\"" 
    #m2zargs      += " metalRug=\"Immunochip\"" 
    #m2zargs      += " theme=\"publication\"" 

    if (snp == chromosome+":"+pos):
        fh_w.writelines("chr" +snp +sep+\
              chromosome +sep+\
              start +sep+\
              end +sep+\
              flank +sep+\
              run +sep+\
              m2zargs +"\n")
    elif (snp[0:3] == "imm" or snp[0:3] == "seq"):
        fh_w.writelines("chr" +chromosome+":"+pos +sep+\
              chromosome +sep+\
              start +sep+\
              end +sep+\
              flank +sep+\
              run +sep+\
              m2zargs +"\n")
    else:
        fh_w.writelines(snp +sep+\
              chromosome +sep+\
              start +sep+\
              end +sep+\
              flank +sep+\
              run +sep+\
              m2zargs +"\n")

    fh_w.close()




def generateCredibleSets(chromosome, pos, snp, pval, file_CredibleSets):
    """ """ 
    
    sep = "\t"
    #counter2color = ["black", "red", "green3", "blue", "cyan", "magenta", "yellow", "gray"]
    counter2color = ["black", "red", "green3", "blue", "cyan", "magenta", "yellow", "gray", "antiquewhite4", "aquamarine4", "azure4", "bisque4", "blueviolet", "brown4", "cadetblue4", "chartreuse", "chocolate4", "cornflowerblue", "cornsilk4", "cyan4", "darkgoldenrod4", "darkmagenta", "darkolivegreen4", "darkorange1", "darkorange4", "darkorchid", "darkorchid4", "darkred", "blue4", "aquamarine1", "brown2", "cornsilk2"]
    snp_counter = 0

    try:
        fh_w = open(file_CredibleSets, "w")
    except IOError, e:
        print e
        sys.exit(1)
    
    fh_w.writelines("snp" +sep+\
          "chr" +sep+\
          "pos" +sep+\
          "pp" +sep+\
          "p-value" +sep+\
          "group" +sep+\
          "color\n")

    snp_counter = 0
    pp  = "NA"

    if (snp == chromosome+":"+pos):
        fh_w.writelines("chr" +snp +sep+\
              chromosome +sep+\
              pos +sep+\
              pp +sep+\
              pval +sep+\
              snp +sep+\
              counter2color[snp_counter] + "\n")
    elif (snp[0:3] == "imm" or snp[0:3] == "seq"):
        fh_w.writelines("chr" +chromosome+":"+pos +sep+\
              chromosome +sep+\
              pos +sep+\
              pp +sep+\
              pval +sep+\
              snp +sep+\
              counter2color[snp_counter] + "\n")
    else:
        fh_w.writelines(snp +sep+\
              chromosome +sep+\
              pos +sep+\
              pp +sep+\
              pval +sep+\
              snp +sep+\
              counter2color[snp_counter] + "\n")
    
    ####LDvariants = sp2.replace("_typed(1)","").split(",")
    ####for variant in LDvariants:
    ####    fh_w.writelines(variant +sep+\
    ####          chromosome +sep+\
    ####          pos +sep+\
    ####          pp +sep+\
    ####          pval +sep+\
    ####          snp +" and LD variants"+ +sep+\
    ####          counter2color[snp_counter] + "\n")
        
    ####snp_counter += 1
    fh_w.close()




def generateDenoteMarker(chromosome, pos, snp, file_DenoteMarker):
    """ """ 

    sep = "\t"
    
    try:
        fh_w = open(file_DenoteMarker, "w")
    except IOError, e:
        print e
        sys.exit(1)
    
    fh_w.writelines("snp" +sep+\
          "chr" +sep+\
          "pos" +sep+\
          "string" +sep+\
          "color\n")
    
    if (snp == chromosome+":"+pos):
        fh_w.writelines("chr" +snp +sep+\
              chromosome +sep+\
              pos +sep+\
              "" +sep+\
              "black\n")
    elif (snp[0:3] == "imm" or snp[0:3] == "seq"):
        fh_w.writelines("chr" +chromosome+":"+pos +sep+\
              chromosome +sep+\
              pos +sep+\
              "" +sep+\
              "black\n")
    else:
        fh_w.writelines(snp +sep+\
              chromosome +sep+\
              pos +sep+\
              "" +sep+\
              "black\n")
    
    fh_w.close()




def run_locuszoom(chromosome, file_PLINK_assoc, file_HitSpec, file_CredibleSets, file_DenoteMarker, prefix_out, file_PLINK_ld):
    """ """

    # workaround: Because Locuszoom tries to display all loci in file_CredibleSets, I need a separate input for each locus
    try:
        fh_r1 = open(file_CredibleSets, "r")
        fh_r2 = open(file_HitSpec, "r")
        fh_r3 = open(file_DenoteMarker, "r")
        # write new tmp files
        fh_w1_tmp = open(file_CredibleSets +".tmp", "w")
        fh_w2_tmp = open(file_HitSpec +".tmp", "w")
        fh_w3_tmp = open(file_DenoteMarker +".tmp", "w")
    except IOError, e:
        print e
        sys.exit(1)
    
    # read header
    header_1 = fh_r1.readline().rstrip('\n')
    header_2 = fh_r2.readline().rstrip('\n')
    header_3 = fh_r3.readline().rstrip('\n')
    header = True
    line_1 = header_1
    while line_1:
    
        if header: 
            # print header
            fh_w1_tmp.writelines(header_1 +"\n")
            fh_w2_tmp.writelines(header_2 +"\n")
            fh_w3_tmp.writelines(header_3 +"\n")
    
            # print first line
            line_1   = fh_r1.readline().rstrip('\n')
            line_2   = fh_r2.readline().rstrip('\n')
            line_3   = fh_r3.readline().rstrip('\n')
            fh_w1_tmp.writelines(line_1 +"\n")
            fh_w2_tmp.writelines(line_2 +"\n")
            fh_w3_tmp.writelines(line_3 +"\n")
            list = re.split("\s+",line_1)
            locus_nr_prev = list[5]
            header = False
        else:
            list = re.split("\s+",line_1)
            locus_nr = list[5]
            if locus_nr == locus_nr_prev:
                fh_w1_tmp.writelines(line_1 +"\n")
                fh_w2_tmp.writelines(line_2 +"\n")
                fh_w3_tmp.writelines(line_3 +"\n")
                locus_nr_prev = locus_nr
            else:
                # close current file and plot with Locuszoom
                fh_w1_tmp.close()
                fh_w2_tmp.close()
                fh_w3_tmp.close()
                os.system("locuszoom --metal=%s --markercol SNP --pvalcol P --ld %s --build hg19 --hitspec=%s fineMap=\"%s\" --denote-markers-file %s --gwas-cat whole-cat_significant-only --plotonly --no-date --prefix %s" %(file_PLINK_assoc, file_PLINK_ld, file_HitSpec +".tmp", file_CredibleSets +".tmp", file_DenoteMarker +".tmp", prefix_out))
                #os.system("locuszoom --metal=%s --markercol SNP --pvalcol P --ld-vcf %s --build hg19 --hitspec=%s fineMap=\"%s\" --denote-markers-file %s --gwas-cat whole-cat_significant-only --plotonly --no-date --prefix %s" %(file_PLINK_assoc, join(Imputation_orig_dir, str(chromosome)+"."+disease_data_set_suffix_release_imputed +".gz"), file_HitSpec +".tmp", file_CredibleSets +".tmp", file_DenoteMarker +".tmp", prefix_out))
                # vcf file tabixed from filtered PLINK data does not work, don't know why, but LD cannot be calculated by Locuszoom
                #os.system("locuszoom --metal=%s --markercol SNP --pvalcol P --ld-vcf %s --build hg19 --hitspec=%s fineMap=\"%s\" --denote-markers-file %s --gwas-cat whole-cat_significant-only --plotonly --no-date --prefix %s" %(file_PLINK_assoc, file_PLINK_rsq0_4, file_HitSpec +".tmp", file_CredibleSets +".tmp", file_DenoteMarker +".tmp", prefix_out))
    
                # write new tmp files
                try:
                    fh_w1_tmp = open(file_CredibleSets +".tmp", "w")
                    fh_w2_tmp = open(file_HitSpec +".tmp", "w")
                    fh_w3_tmp = open(file_DenoteMarker +".tmp", "w")
                except IOError, e:
                    print e
                    sys.exit(1)
    
                # print header
                fh_w1_tmp.writelines(header_1 +"\n")
                fh_w2_tmp.writelines(header_2 +"\n")
                fh_w3_tmp.writelines(header_3 +"\n")
                
                fh_w1_tmp.writelines(line_1 +"\n")
                fh_w2_tmp.writelines(line_2 +"\n")
                fh_w3_tmp.writelines(line_3 +"\n")
    
                list = re.split("\s+",line_1)
                locus_nr_prev = list[5]
    
        line_1   = fh_r1.readline().rstrip('\n')
        line_2   = fh_r2.readline().rstrip('\n')
        line_3   = fh_r3.readline().rstrip('\n')
    
    # close current file and plot with Locuszoom
    fh_w1_tmp.close()
    fh_w2_tmp.close()
    fh_w3_tmp.close()
    os.system("locuszoom --metal=%s --markercol SNP --pvalcol P --ld %s --build hg19 --hitspec=%s fineMap=\"%s\" --denote-markers-file %s --gwas-cat whole-cat_significant-only --plotonly --no-date --prefix %s" %(file_PLINK_assoc, file_PLINK_ld, file_HitSpec +".tmp", file_CredibleSets +".tmp", file_DenoteMarker +".tmp", prefix_out))
    #os.system("locuszoom --metal=%s --markercol SNP --pvalcol P --ld-vcf %s --build hg19 --hitspec=%s fineMap=\"%s\" --denote-markers-file %s --gwas-cat whole-cat_significant-only --plotonly --no-date --prefix %s" %(file_PLINK_assoc, join(Imputation_orig_dir, str(chromosome)+"."+disease_data_set_suffix_release_imputed +".gz"), file_HitSpec +".tmp", file_CredibleSets +".tmp", file_DenoteMarker +".tmp", prefix_out))
    # vcf file tabixed from filtered PLINK data does not work, don't know why, but LD cannot be calculated by Locuszoom
    #os.system("locuszoom --metal=%s --markercol SNP --pvalcol P --ld-vcf %s --build hg19 --hitspec=%s fineMap=\"%s\" --denote-markers-file %s --gwas-cat whole-cat_significant-only --plotonly --no-date --prefix %s" %(file_PLINK_assoc, file_PLINK_rsq0_4, file_HitSpec +".tmp", file_CredibleSets +".tmp", file_DenoteMarker +".tmp", prefix_out))
    os.system("rm -f %s" %(file_CredibleSets +".tmp"))
    os.system("rm -f %s" %(file_HitSpec +".tmp"))
    os.system("rm -f %s" %(file_DenoteMarker +".tmp"))
    fh_r1.close()
    fh_r2.close()
    fh_r3.close()



