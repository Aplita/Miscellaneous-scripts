#!/bin/bash

################################################################################
# Run MHC-2 epitope prediction     #
####################################
# This script simply runs locally IEDB epitope prediction based on a set of
# predefined alleles and a .fasta document with protein sequences.
#
# Written by Ana Paula Vargas R.
# Bioinformatics lab - UPCH. Lima, Peru.
# Last updated: 05-jan-20.
###############################################################################

# Set vars:
ALLS='HLA-DRB1*01:01,HLA-DRB1*01:03,HLA-DRB1*03:01,HLA-DRB1*04:01,HLA-DRB1*04:02,HLA-DRB1*04:03,HLA-DRB1*04:04,HLA-DRB1*04:05,HLA-DRB1*07:01,HLA-DRB1*08:01,HLA-DRB1*08:02,HLA-DRB1*09:01,HLA-DRB1*10:01,HLA-DRB1*11:01,HLA-DRB1*12:01,HLA-DRB1*13:01,HLA-DRB1*13:02,HLA-DRB1*15:01,HLA-DRB1*16:02,HLA-DRB3*01:01,HLA-DRB3*02:02,HLA-DRB3*03:01,HLA-DRB4*01:01,HLA-DRB4*01:03,HLA-DRB5*01:01,HLA-DPA1*01/DPB1*04:01,HLA-DPA1*01:03/DPB1*02:01,HLA-DPA1*02:01/DPB1*05:01,HLA-DPA1*02:01/DPB1*01:01,HLA-DPA1*03:01/DPB1*04:02,HLA-DQA1*01:01/DQB1*05:01,HLA-DQA1*01:02/DQB1*06:02,HLA-DQA1*03:01/DQB1*03:02,HLA-DQA1*04:01/DQB1*04:02,HLA-DQA1*05:01/DQB1*02:01,HLA-DQA1*05:01/DQB1*03:01'
METHOD="IEDB_recommended"


###############
# Execute:
###############

# MSP7C:
python ~/Downloads/mhc_ii/./mhc_II_binding.py $METHOD $ALLS MSP7C_all_Prot.fasta 12-18 > IEDB_12-18mer_MSP7C.txt

# MSP7H:
python ~/Downloads/mhc_ii/./mhc_II_binding.py $METHOD $ALLS MSP7H_all_Prot.fasta 12-18 > IEDB_12-18mer_MSP7H.txt

# MSP8:
python ~/Downloads/mhc_ii/./mhc_II_binding.py $METHOD $ALLS MSP8_all_Prot.fasta 12-18 > IEDB_12-18mer_MSP8.txt
