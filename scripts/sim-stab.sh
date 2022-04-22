#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=0
#SBATCH --partition=bdw
#SBATCH --mail-user=mazzaris@mis.mpg.de
#SBATCH --mail-type=ALL


julia phase-diagram.jl 
