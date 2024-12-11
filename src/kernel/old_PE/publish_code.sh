#!/bin/bash
while getopts m: flag
do
    case "${flag}" in
        m) message=${OPTARG};;
    esac
done

git add .
if [$message == ""]; then
   git commit -m 'new version of PE_2D'
else
  git commit -m "$message"
fi

git push
ssh haswell '( cd Documents/DEV/wf_phd ; git pull )'
ssh haswell '(module purge ; module load  HDF5/1.10.1-intel-2018a ; cd Documents/DEV/wf_phd/src/kernel/New_PE; make intel)'
