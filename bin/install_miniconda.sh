#!/bin/bash

mkdir -p tmp
wget -P tmp \
     https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash tmp/Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
. $HOME/miniconda3/bin/activate
conda init
