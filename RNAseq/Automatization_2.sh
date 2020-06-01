#!/bin/bash

################################################################################
# Automatization 2      #
#########################
# This script automatizes the de novo assembly of test MTB sequences
# for RNA-seq analysis.
#
# Written by Ana Paula Vargas.
# Bioinformatics lab - UPCH. Lima, Peru.
# Last updated: 10-ago-18.
###############################################################################


# Piloto/
# $1: FileName
# $2: RootName + path

PATH=$PATH:~/Downloads/oases/
export PATH
source ~/.bashrc
PATH=$PATH:~/Downloads/velvet/
export PATH
source ~/.bashrc
PATH=$PATH:~/Downloads/bowtie2-2.3.4.1-linux-x86_64/
export PATH
source ~/.bashrc
PATH=$PATH:~/Downloads/cufflinks-2.2.1.Linux_x86_64/
export PATH
source ~/.bashrc
PATH=$PATH:~/Downloads/tophat-2.1.1.Linux_x86_64/
export PATH
source ~/.bashrc

Analysis() {
	# Create support variables.
	TrimName="$1.trim.fastq"
	FTrimName="$1.trim.fasta"
	BowName="Piloto-$1"

	# Create and prepare directories.
	mkdir $1
	# cp $2 $1
	# Not practical to copy file. Instead I'll use full path to file.
	cd $1

	# Programs.
	echo "Trimming"
	TrimmomaticSE -phred33 $2 $TrimName ILLUMINACLIP:TruSeq2-SE.fa:2:30:10 LEADING:19 TRAILING:19 SLIDINGWINDOW:4:15 MINLEN:34
	echo "Finished trimming"
	echo ""
	echo "Converting files to fasta"
        fastq_to_fasta -Q33 -n -i $TrimName -o $FTrimName
	echo "Finished converting"
	echo ""
	echo "Excecuting Oases"
        python2 ~/workshop/oases/scripts/oases_pipeline.py -m 21 -M 23 -d $FTrimName
	echo "Oases ready"
	echo ""
	echo "Excecuting Bowtie"
	mkdir BowResults
	cd oasesPipelineMerged
        bowtie2-build transcripts.fa /home/vargas/Desktop/RNA-seq/Piloto/$1/BowResults/$BowName
	echo "Finished Bowtie"
	echo ""
	cd ..
        mkdir THResults
	echo "Excecuting TopHat"
        tophat2 -o THResults/ -g 1 --library-type fr-firststrand -x 1 -p 5 /home/vargas/Desktop/RNA-seq/Piloto/$1/BowResults/$BowName $TrimName &> THResults/run_tophat.log
	echo "Finished TopHat"
	echo ""
        cd THResults
	echo "Excecuting SamTools"
        samtools view -h -o Hits.sam accepted_hits.bam
	echo "Finished SamTools"
	echo ""
        cd ..
        mkdir CuffOut
	echo "Excecuting Cufflinks"
        cufflinks -p 5 --library-type fr-firststrand -o CuffOut/ THResults/Hits.sam &> CuffOut/CuffHits.log
	echo "Finished Cufflinks"
	echo ""
        cd ..

}



shopt -s nullglob; # Make sure doesn't expand.
for File in *.fq; do
	echo "PROCESSING FILE $File"
	FileName="$(basename $File .fq)"; # Remove .fq extension. Now file name is ready and stored in FileName.
	RName="$File"; # Name, extension included.
	RPName="$PWD/$RName"; # Name + path + extension.
	Analysis $FileName $RPName; # Pass FileName and RootName to function.
	echo "FINISHED PROCESSING FILE $File"
	echo ""
	echo ""
done
