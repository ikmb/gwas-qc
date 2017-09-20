#!/usr/bin/env python

import sys
import re
import os

#from os.path import *
#import string

# import gzip
# import math
# import decimal
# import datetime
# from os import listdir
# import subprocess

# may also need some of these:

# import Ingos lib
sys.path.append(join(sys.path[0], "../../all_scripts"))
sys.path.append(os.environ['PYLIB_DIR'] + "/all_scripts")
from all_common import Command

# import my lib
# sys.path.append(join(sys.path[0], "../lib"))
# sys.path.append(os.environ['PYLIB_DIR'] + "/lib")

from plink_classes import PackedPed
# from eigenstrat_classes import *

def pca_convert(plink, eigenstrat_parameter_file, annotation_file):
    """ convert PLINK file data set to eigenstrat format """

    # ----------------------------- #
    # - generate parameter file m - #
    # ----------------------------- #

    packedped = PackedPed(write_file=eigenstrat_parameter_file)

    packedped.set_input_PLINK_binary(
        bed=plink + ".bed",
        bim=plink + ".bim",
        fam=plink + ".fam")

    packedped.write_par_file()
    del packedped

    # ------------------------ #
    # - run convertf program - #
    # ------------------------ #

    cmd = Command("convertf -p %s" % (eigenstrat_parameter_file))
    cmd.run()
    del cmd

    os.system("mv %s.ind %s.ind.orig" % (plink, plink))

    # read individualIDs and HAPMAP info from from hapmap2 fam file
    try:
        fh_anno = file(annotation_file, "r")
    except IOError, e:
        print e
        sys.exit(1)

    individuals2batch_id = {}
    # skip header
    line = fh_anno.readline().replace("\n", "")
    line = fh_anno.readline().replace("\n", "")
    while line:

        list = re.split("\s+", line)
        IID = list[1]
        batch_id = list[6]
        individuals2batch_id[IID] = batch_id

        line = fh_anno.readline().replace("\n", "")

    fh_anno.close()

    # re-write ind file with info on HapMap samples and batch_info
    try:
        fh_ind     = file(plink + ".ind.orig", "r")
        fh_ind_new = file(plink + ".ind", "w")
    except IOError, e:
        print e
        sys.exit(1)

    batches = []
    batches_dict = {}

    # no header line
    line = fh_ind.readline().replace("\n", "")
    while line:

        list = re.split("\s+", line)
        if list[0] == "":
            del list[0]

        # change info last column from "Case/Control" to batch_id
        if list[0] in individuals2batch_id:

            batch_id = individuals2batch_id[list[0]]
            if batch_id in batches_dict.has_key:
                batches.append(batch_id)
                batches_dict[batch_id] = True
            if list[-1] == "Case":
                line = line.replace("Case", batch_id)
            elif list[-1] == "Control":
                line = line.replace("Control", batch_id)
            # nothing to replace
            else:
                print >> sys.stderr, "\n    warning: could not replace case/control status for sample " + list[0] + " by batch_id in file pca.evec file " + plink_pca + ".pca.evec ...\n\n"
            fh_ind_new.writelines(line + "\n")

        # nothing to replace
        else:
            print >> sys.stderr, "\n    warning: could not found sample " + list[0] + " in annotation file " + individuals_annotation_cases_controls_hapmap2 + " ...\n\n"
            fh_ind_new.writelines(line + "\n")

        line = fh_ind.readline().replace("\n", "")

    fh_ind.close()
    fh_ind_new.close()
    del batches_dict


pca_main_program = "smartpca.perl.DE"

def pca_run(plink, sigmathreshold, projection_on_populations, numof_pc, numof_threads=1):
    """ run eigenstrat program """

    # ------------------------ #
    # - run eigenstrat program - #
    # ------------------------ #

    plink_pca = plink + "_" + numof_pc + "PC"

    teststring = "%s -i %s.eigenstratgeno -a %s.snp -b %s.ind -k %s -o %s.pca -p %s.plot -e %s.eval -l %s.log -m 5 -t %s -s %s -w %s -f %s -g %s.snpweights" \
                 % (pca_main_program,
                    plink,
                    plink,
                    plink,
                    numof_pc,
                    plink_pca,
                    plink_pca,
                    plink_pca,
                    plink_pca,
                    numof_pc,
                    sigmathreshold,
                    projection_on_populations,
                    numof_threads,
                    plink_pca)
    print >> sys.stderr, teststring
    cmd = Command("%s -i %s.eigenstratgeno -a %s.snp -b %s.ind -k %s -o %s.pca -p %s.plot -e %s.eval -l %s.log -m 5 -t %s -s %s -w %s -f %s -g %s.snpweights"
                  % (pca_main_program,
                     plink,
                     plink,
                     plink,
                     numof_pc,
                     plink_pca,
                     plink_pca,
                     plink_pca,
                     plink_pca,
                     numof_pc,
                     sigmathreshold,
                     projection_on_populations,
                     numof_threads,
                     plink_pca))
    cmd.run()
    del cmd

    # draw first two PCs
    os.system("R --slave --args %s < %s" % (plink_pca, "raw_evec_EIGENSTRAT.r"))

    # read which batches (HapMap) were used for projection
    projection_batches = {}
    try:
        fh_proj     = file(projection_on_populations, "r")
    except IOError, e:
        print e
        sys.exit(1)

    line = fh_proj.readline().rstrip('\n')
    while line:

        list = re.split("\s+", line)
        if list[0] == "":
            del list[0]

        projection_batches[list[0]] = list[0]
        line = fh_proj.readline().rstrip('\n')

    fh_proj.close()

    # re-write pca.evec file without projection samples (HapMap samples)
    try:
        fh_pcaevec     = file(plink_pca + ".pca.evec", "r")
        fh_pcaevec_new = file(plink_pca + ".withoutProjection.pca.evec", "w")
    except IOError, e:
        print e
        sys.exit(1)

    # skip header line
    line = fh_pcaevec.readline().rstrip('\n')

    line = fh_pcaevec.readline().rstrip('\n')
    while line:

        list = re.split("\s+", line)
        if list[0] == "":
            del list[0]
        if list[-1] == "":
            del list[-1]

        if not list[-1] in projection_batches:
            for i in xrange(len(list)):
                if i == 0:
                    fh_pcaevec_new.writelines(list[i])
                else:
                    fh_pcaevec_new.writelines("\t" + list[i])
            fh_pcaevec_new.writelines("\n")

        line = fh_pcaevec.readline().rstrip('\n')

    fh_pcaevec.close()
    fh_pcaevec_new.close()

    # draw first two PCs without HapMap samples
    os.system("R --slave --args %s %s < %s" % (plink_pca, plink_pca + ".withoutProjection", "draw_evec_withoutProjection.r"))


def write_snps_autosomes_noLDRegions_noATandGC_noIndels(bim, outfile):
    """  write only autosomal snps, remove SNPs from high LD regions (also MHC),
    remove A/T and C/G SNPs, remove Indels """

    print "\n        remove SNPs from high LD regions ..\n\n"
    print "\n        remove A/T and C/G SNPs ...\n\n"
    print "\n        remove insertions/deletions ...\n\n"

    try:
        bim = file(bim, "r")
        out = file(outfile, "w")
    except IOError, e:
        print e
        sys.exit(1)

    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'D': 'D', 'I': 'I'}
    indels = {'D': 'D', 'I': 'I'}

    line = bim.readline().replace("\n", "")
    while line:

        list = re.split("\s+", line)
        chr  = int(list[0])
        pos  = int(list[3])
        a1   = list[4]
        a2   = list[5]
        # exclude non-autosomes
        if 0 < chr and chr < 23:
            # exclude xMHC SNPs AND exclude A/T and C/G SNPs AND exclude D/I SNPs
            if (not ((1 == chr and (48000000 <= pos and pos < 52000000)) or
                     (2 == chr and (86000000 <= pos and pos < 100500000)) or
                     (2 == chr and (134500000 <= pos and pos < 138000000)) or
                     (2 == chr and (183000000 <= pos and pos < 183000000)) or
                     (3 == chr and (47500000 <= pos and pos < 50000000)) or
                     (3 == chr and (83500000 <= pos and pos < 87000000)) or
                     (3 == chr and (89000000 <= pos and pos < 97500000)) or
                     (5 == chr and (44500000 <= pos and pos < 50500000)) or
                     (5 == chr and (98000000 <= pos and pos < 100500000)) or
                     (5 == chr and (129000000 <= pos and pos < 132000000)) or
                     (5 == chr and (135500000 <= pos and pos < 138500000)) or
                     (6 == chr and (25500000 <= pos and pos < 33500000)) or
                     (6 == chr and (57000000 <= pos and pos < 64000000)) or
                     (6 == chr and (140000000 <= pos and pos < 142500000)) or
                     (7 == chr and (55000000 <= pos and pos < 66000000)) or
                     (8 == chr and (8000000 <= pos and pos < 12000000)) or
                     (8 == chr and (43000000 <= pos and pos < 50000000)) or
                     (8 == chr and (112000000 <= pos and pos < 115000000)) or
                     (10 == chr and (37000000 <= pos and pos < 43000000)) or
                     (11 == chr and (46000000 <= pos and pos < 57000000)) or
                     (11 == chr and (87500000 <= pos and pos < 90500000)) or
                     (12 == chr and (33000000 <= pos and pos < 40000000)) or
                     (12 == chr and (109500000 <= pos and pos < 112000000)) or
                     (20 == chr and (32000000 <= pos and pos < 34500000))))\
                    and (a1 != complement[a2]) and (not (a1 in indels or a2 in indels)):

                # write variants for inclusion
                out.writelines("%s\n" % (list[1]))

        line = bim.readline().replace("\n", "")

    bim.close()
    out.close()
