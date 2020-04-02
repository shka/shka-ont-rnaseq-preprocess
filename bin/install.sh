#!/bin/bash -eu

conda install --file environment.yaml
R CMD BATCH bin/install.R
emacs -l bin/install.el --batch --kill
