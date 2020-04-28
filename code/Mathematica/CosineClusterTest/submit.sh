#!/bin/sh

# Give the job a name
#$ -N ClusterTest

# set the shell
#$ -S /bin/sh

# set working directory on all host to
# directory where the job was started
#$ -cwd

# send all process STDERR (fd 2) to this file
#$ -o output.txt

# send all process STDERR (fd 3) to this file
#$ -e errors.txt

# email information
#$ -m e
 
# Just change the email address. You will be emailed when the job has finished.
#$ -M mark.novak@oregonstate.edu

# generic parallel environment with 1 core(s) requested
#$ -pe orte 1

# Load a module, if needed
module load mathematica/12.0

# Commands
math -script ClusterTest.m