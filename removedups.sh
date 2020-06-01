#!/bin/bash

################################################################################
# Remove extension duplicates     #
###################################
# This script reads all files in a directory and removes those with duplicate
# extension (PE.PE) from a bad analysis.
#
# Written by Ana Paula Vargas.
# Bioinformatics lab - UPCH. Lima, Peru.
# Last updated: 09-mar-20.
###############################################################################



# Remove duplicate files with PE.PE

cd ~/Desktop/21_2000GenomesProc/results_all/
DirPEPE=`ls *.vcf | grep "PE.PE"`
DirPE=`ls *.vcf | grep -v "PE.PE"`


for f in $DirPEPE; do
  bnamePEPE=`basename $f`
  jnamePEPE=${bnamePEPE::-40}


if [[ "${DirPE[@]}" =~ "${jnamePEPE}" ]]; then
  echo "FOUND: $jnamePEPE"
  rm $f
else
  echo "NOT FOUND: $jnamePEPE"
  mv $f "$jnamePEPE.PE.bwa.MTb.passed.realn.flt.anot.vcf"
fi


done
