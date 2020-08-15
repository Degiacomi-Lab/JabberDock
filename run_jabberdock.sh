#!/bin/bash
# Starting from minimal Ubuntu 18.4 install :)

# get c compiler
sudo apt install -y build-essential cmake wget

# get gromacs
wget -q --show-progress http://ftp.gromacs.org/pub/gromacs/gromacs-2020.2.tar.gz
tar xzf gromacs-2020.2.tar.gz
cd gromacs-2020.2/
mkdir build
cd build
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON >> gromacs-build.log 2>&1
make -j24 >> gromacs-build.log 2>&1
sudo make install >> gromacs-build.log 2>&1
cd ~

# get conda with python 2.7
wget -q --show-progress https://repo.anaconda.com/miniconda/Miniconda2-py27_4.8.3-Linux-x86_64.sh 
bash Miniconda2-py27_4.8.3-Linux-x86_64.sh -b -p $HOME/miniconda 
source $HOME/miniconda/bin/activate
conda init

# install python packages
conda install -y numpy scipy cython pandas scikit-learn scikit-image matplotlib mpi4py dill
conda install -y -c conda-forge vmd

# TODO make environment.yml with these versions
# numpy-1.16.6
# scipy-1.2.1
# cython-0.29.14 
# pandas-0.22.0
# scikit-learn-0.20.3
# scikit-image
# matplotlib-2.2.3
# dill-0.3.2
# mpi4py-3.0.3
# mpich-3.3.2 
# vmd-1.9.3

# NOTE fixes in this version of JabberDock are probably needed for gromacs > 5.x
git clone https://github.com/aizvorski/JabberDock

# NOTE POW_v2 and biobox have to be linked in the home directory, location is hardwired in a few places
ln -s JabberDock/POW_v2 ~/POW_v2
ln -s JabberDock/biobox ~/biobox

# NOTE this doesn't install with setuptools, but it does cythonize
cd JabberDock && python setup.py install && cd ~
cd biobox && python setup.py install && cd ~

# NOTE the python path needs to be set too; eg import biobox works because home directory is in the python path
echo 'export PATH=$PATH:/home/ubuntu/JabberDock/auto_scripts' >> ~/.bashrc
echo 'export PYTHONPATH=/home/ubuntu/JabberDock/:/home/ubuntu/' >> ~/.bashrc
source ~/.bashrc
source $HOME/miniconda/bin/activate

mkdir jd_tutorial
cd jd_tutorial
cp ~/JabberDock/tutorial/*pdb ./

jabberdock.py -ir 1dfj_0.pdb -il 1dfj_1.pdb -np 12 -ntomp 4
