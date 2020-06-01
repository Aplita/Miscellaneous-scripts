################################################################################
# Install RNAseq programs      #
################################
# Just a list of programs needed for analysis.
# Optimized installation for a workshop.
#
# Written by Ana Paula Vargas.
# Bioinformatics lab - UPCH. Lima, Peru.
# Last updated: 10-ago-18.
###############################################################################


sudo apt-get update
cd
mkdir workshop
cd workshop

# Confirm git is installed
sudo apt install git

# Java Run Environment.
# sudo apt-get install default-jre
wget --no-cookies --no-check-certificate --header "Cookie: gpw_e24=http%3A%2F%2Fwww.oracle.com%2F; oraclelicense=accept-securebackup-cookie" "http://download.oracle.com/otn-pub/java/jdk/8u144-b01/090f390dda5b47b9b721c7dfaa008135/jdk-8u144-linux-x64.tar.gz"
tar -xvzf jdk-8u144-linux-x64.tar.gz
sudo install /usr/bin/java java
cd ~/workshop

# Trimmomatic
sudo apt install trimmomatic

# Velvet.
git clone --recursive git://github.com/dzerbino/velvet.git
cd velvet
make
PATH=$PATH:~/workshop/velvet
export PATH
source ~/.bashrc
sudo apt install velvet

# Oases.
cd ~/workshop
git clone --recursive git://github.com/dzerbino/oases.git
cd oases
make BUNDLEDZLIB=1
PATH=$PATH:~/workshop/oases
export PATH
source ~/.bashrc
cd ~/workshop

# FastX Toolkit: Fastq_to_Fasta
sudo apt install fastx-toolkit

# Bowtie
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.3/bowtie2-2.3.3-linux-x86_64.zip
unzip bowtie2-2.3.3-linux-x86_64.zip
cd bowtie2-2.3.3
PATH=$PATH:~/workshop/bowtie2-2.3.3
export PATH
source ~/.bashrc
cd ~/workshop

# Zlib
sudo apt-get install zlib1g-dev libncurses5-dev

# Cmake
sudo apt-get install cmake

# Samtools
sudo apt-get install samtools

# Boost Libraries
sudo apt-get install libboost-all-dev

# Eigen Libraries
wget http://bitbucket.org/eigen/eigen/get/3.3.4.tar.gz
tar -xvzf 3.3.4.tar.gz
mkdir build_eigen_dir
cd build_eigen_dir
cmake ~/workshop/eigen-eigen-5a0156e40feb
sudo make install
source ~/.bashrc
cd ~/workshop/

# TopHat
wget https://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.1.Linux_x86_64.tar.gz
tar xvfz tophat-2.1.1.Linux_x86_64.tar.gz
PATH=$PATH:~/workshop/tophat-2.1.1.Linux_x86_64
export PATH
source ~/.bashrc
cd ~/workshop/


# Cufflinks
wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz
tar -xvzf cufflinks-2.2.1.Linux_x86_64.tar.gz
PATH=$PATH:~/workshop/cufflinks-2.2.1.Linux_x86_64
export PATH
source ~/.bashrc
