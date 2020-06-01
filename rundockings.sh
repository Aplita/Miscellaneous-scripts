#!/bin/bash

################################################################################
# Run dockings            #
###########################
# Script to run a series of dockings.
# All structures had previously been prepared.
# 3 different grid files prepared previously.
# Total of 3 x 7 = 21 dockings.
#
# Written by Ana Paula Vargas.
# Bioinformatics lab - UPCH. Lima, Peru.
# Last updated: 31-may-20.
###############################################################################

# Get filenames
Dir=/home/apa/Desktop/17_CLPC1/06_MutDockings/02_PreparedStructures/*.pdbqt

declare -a REC
for file in $Dir; do
  bname=`basename $file`
  jname=${bname::-6}
  REC=("${REC[@]}" "$jname")
done

# Remove ligand from that list
unset REC[7]

# Iterate over files to dock
DN=0

while [ $DN -le 6 ]; do
  N=$(( $DN + 1 ))
  ID="d$N"

# Generate .gpf file from template, according to docking
# Templates are 3 to modify: a, b, c
echo "Generating GPF from templates. Docking N $ID"
sed "s/3wdb_prot_min/${REC[DN]}/g" gridA_lit.gpf > "$ID"_A.gpf
sed "s/3wdb_prot_min/${REC[DN]}/g" gridB_meta1.gpf > "$ID"_B.gpf
sed "s/3wdb_prot_min/${REC[DN]}/g" gridC_meta2.gpf > "$ID"_C.gpf
echo "Finished generating GPF"

# Run autogrid4
echo "Running autogrid - A. Docking ID: $ID"
autogrid4 -p "$ID"_A.gpf -l "$ID"_A.glg
echo "Running autogrid - B. Docking ID: $ID"
autogrid4 -p "$ID"_B.gpf -l "$ID"_B.glg
echo "Running autogrid - C. Docking ID: $ID"
autogrid4 -p "$ID"_C.gpf -l "$ID"_C.glg
echo "Finished autogrid"

# Generate .dpf file from template
echo "Generating DPF from templates. Docking ID: $ID"
sed "s/3wdb_prot_min/${REC[DN]}/g" template.dpf > "$ID"_A.dpf
sed "s/3wdb_prot_min/${REC[DN]}/g" template.dpf > "$ID"_B.dpf
sed "s/3wdb_prot_min/${REC[DN]}/g" template.dpf > "$ID"_C.dpf
echo "Finished generating DPF"

# Run autodock4
echo "Running autodock - A. Docking ID: $ID"
autodock4 -p "$ID"_A.dpf -l "$ID"_A.glg
echo "Running autodock - B. Docking ID: $ID"
autodock4 -p "$ID"_B.dpf -l "$ID"_B.glg
echo "Running autodock - C. Docking ID: $ID"
autodock4 -p "$ID"_C.dpf -l "$ID"_C.glg
echo "Finished autodock"


DN=$(( $DN + 1 ))
done
