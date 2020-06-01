################################################################################
# Automatization 1      #
#########################
# This script automatizes the de novo assembly of test MTB sequences
# for RNA-seq analysis.
#
# Written by Ana Paula Vargas.
# Bioinformatics lab - UPCH. Lima, Peru.
# Last updated: 14-jun-18.
###############################################################################

proc gen_analysis {FileName RootName} {

# Create support vars
        set TrimName "$FileName.trim.fastq"
        set FTrimName "$FileName.trim.fasta"
        set BowName "Piloto-$FileName"
        puts $TrimName

# Create and prepare a directory
        file mkdir $FileName
        file copy $RootName "$FileName/"
        cd "$FileName/"

# Execute programs
        TrimmomaticSE -phred33 $FileName $TrimName ILLUMINACLIP:adapters/adapters/TruSeq2-SE.fa:2:30:10 LEADING:19 TRAILING:19 SLIDINGWINDOW:4:15 MINLEN:34 CROP:123 HEADCROP:30
        fastq_to_fasta -Q33 -n -i $TrimName -o $FTrimName
        python oases_0.2.8/scripts/oases_pipeline.py -m 21 -M 23 -d $FTrimName
        bowtie2-build oasesPipelineMerged/transcripts.fa $BowName
        file mkdir THResults
        tophat2 -o THResults/ -g 1 --library-type fr-firststrand -x 1 -p 5 bowtie2transc/"Piloto-$FileName" "$FileName.trim.fasta" &> THResults/run_tophat.log
        cd THResults/
        samtools view -h -o accepted_hits.bam Hits.sam
        cd ..
        file mkdir CuffOut
        cufflinks -p 5 --library-type fr-firststrand -o CuffOut/ THResults/Hits.sam &> CuffOut/CuffHits.log
        cd ..

}


foreach file [glob *.fq] {
        # Prepare names
        set Name [file rootname $file] ; # Remove extension and save in Name
        set RName "$file"; # Name + extension
        gen_analysis $Name $RName; # Execute procedure for each file
}
