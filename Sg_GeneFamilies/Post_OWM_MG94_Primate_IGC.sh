#!/bin/bash
sbatch -o ./Primate-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/ENSG00000120952_ENSG00000182330_MG94_nonclock.sh 
sbatch -o ./Primate-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/ENSG00000182330_ENSG00000187545_MG94_nonclock.sh 
sbatch -o ./Primate-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/ENSG00000204481_ENSG00000204510_MG94_nonclock.sh 
