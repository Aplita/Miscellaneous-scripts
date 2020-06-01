#!/bin/bash

################################################################################
# Move PC               #
#########################
# Just a collection of programs I use on Ubuntu.
#
# Written by Ana Paula Vargas.
# Bioinformatics lab - UPCH. Lima, Peru.
# Last updated: 13-mar-20.
###############################################################################


sudo apt-get update
sudo apt-get upgrade

# Install restricted extras (media codecs)
sudo apt -y install ubuntu-restricted-extras

# Install tweak tools, toggle minimize on click
# Also don't forget to activate multiple workspaces (manual)
sudo apt-get -y install gnome-tweak-tool
gsettings set org.gnome.shell.extensions.dash-to-dock click-action 'minimize'

# Activate preferences on gedit
gsettings set org.gnome.gedit.preferences.editor display-line-numbers true
# Also list them if needed:
# gsettings list-recursively | grep -i gedit
# Links to download themes:
# https://github.com/mig/gedit-themes/
# http://scribes.sourceforge.net/themegenerator.php
# Last used theme: DESERT (from GitHub)

# Essentials
sudo apt -y install build-essential zlib1g-dev libncurses5-dev libgdbm-dev libnss3-dev libssl-dev libreadline-dev libffi-dev libgtk-3-dev wget gcc

# Anaconda (download first)
sudo apt -y install python3.7
cd Downloads/
chmod +x Anaconda3-2019.10-Linux-x86_64.sh
bash Anaconda3-2019.10-Linux-x86_64.sh
source ~/.bashrc
conda update conda
conda update anaconda

# R and RStudio (download first)
# First check if new version, this is for 3.6
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
# Libraries for bioconductor
sudo apt-get install curl
sudo apt-get install libssl-dev
sudo apt-get install libcurl4-openssl-dev
sudo apt-get install libxml2-dev

sudo apt -y install r-base
sudo dpkg -i rstudio-1.2.5033-amd64.deb
sudo apt-get install -f


# BLAST
sudo apt-get -y install ncbi-blast+

# Foxit reader (download first)
tar -xvzf FoxitReader.enu.setup.2.4.4.0911.x64.run.tar.gz
./FoxitReader.enu.setup.2.4.4.0911\(r057d814\).x64.run

# ClustalOmega & ClustalX
sudo apt-get -y install clustalo
sudo apt-get -y install clustalx

# Git
sudo apt -y install git
git config --global user.name "Ana Paula Vargas"
git config --global user.email "ana.vargas.r@upch.pe"

# Htop
sudo apt install htop

# Docking
sudo apt-get -y install autodock autodocktools autogrid

# VIM
sudo apt-get -y install vim

# VMD & NAMD
tar -xvzf vmd-1.9.4a12.bin.LINUXAMD64-CUDA9-OptiX411-OSPRay140.opengl.tar.gz
cd vmd-1.9.4a12/
./configure
cd src/
sudo make install
cd ..
tar -xvzf NAMD_Git-2019-06-02_Linux-x86_64-multicore-CUDA.tar.gz
cd NAMD_Git-2019-06-02_Linux-x86_64-multicore-CUDA/
sudo cp namd2 /usr/local/bin/namd2
cd ..

# JRE + JDK
sudo apt-get install default-jre
sudo apt-get install default-jdk

# Trimmomatic
sudo apt install trimmomatic

# Velvet
sudo apt-get install velvet

# Oases
git clone --recursive git://github.com/dzerbino/oases.git
cd oases
make BUNDLEDZLIB=1
PATH=$PATH:~/Downloads/oases
export PATH
source ~/.bashrc
cd ..

# Bowtie
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.3/bowtie2-2.3.3-linux-x86_64.zip
unzip bowtie2-2.3.3-linux-x86_64.zip
PATH=$PATH:~/Downloads/bowtie2-2.3.3
export PATH
source ~/.bashrc

# FastX Toolkit, Zlib, Cmake, Samtools
sudo apt install fastx-toolkit
sudo apt-get install zlib1g-dev libncurses5-dev
sudo apt-get install cmake samtools libboost-all-dev

# TopHat
wget https://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.1.Linux_x86_64.tar.gz
tar xvfz tophat-2.1.1.Linux_x86_64.tar.gz
PATH=$PATH:~/Downloads/tophat-2.1.1.Linux_x86_64
export PATH
source ~/.bashrc

#Cufflinks
wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz
tar -xvzf cufflinks-2.2.1.Linux_x86_64.tar.gz
PATH=$PATH:~/Downloads/cufflinks-2.2.1.Linux_x86_64
export PATH
source ~/.bashrc
cd ~/Downloads/

# Phylip
sudo apt-get install -y phylip

# FastQC
chmod 755 ~/Downloads/FastQC/fastqc
sudo ln -s ~/Downloads/FastQC/fastqc /usr/local/bin/fastqc


# From Ubuntu Software:
# GNU Image Manipulation Program (GIMP)
# Spotify
# Zoom (download: https://zoom.us/download?os=linux)
# Mendeley Desktop (download: https://www.mendeley.com/guides/download-mendeley-desktop/ubuntu/instructions)
# TeamViewer
#


# System cleanup
sudo apt-get autoclean
sudo apt-get clean
sudo apt-get autoremove
