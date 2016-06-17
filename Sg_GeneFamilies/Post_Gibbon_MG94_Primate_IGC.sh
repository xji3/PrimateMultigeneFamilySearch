#!/bin/bash
sbatch -o ./Primate-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/ENSG00000160072_ENSG00000197785_MG94_nonclock.sh 
sbatch -o ./Primate-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/ENSG00000203950_ENSG00000212747_MG94_nonclock.sh 
