#!/bin/bash -eu

conda env update --file environment.yml
R CMD BATCH bin/install.R
emacs -l bin/install.el --batch --kill
