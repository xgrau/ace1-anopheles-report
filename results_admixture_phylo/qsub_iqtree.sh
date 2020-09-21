#!/bin/bash
#$ -V
#$ -cwd
#$ -M xavier.graubove@crg.eu
#$ -m a
#$ -q long-sl7
#$ -l virtual_free=50G,h_rt=720:00:00
#$ -o tmp/
#$ -e tmp/

# iqtree
iqtree -s $1 -nt AUTO -ntmax $3 -m TEST -mset GTR -pre $2 -bb 1000 -nm 100000 -cptime 3600
