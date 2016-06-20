#!/bin/bash
sbatch -o ./Primate-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/ENSG00000120952_ENSG00000204505_MG94_nonclock.sh 
sbatch -o ./Primate-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/ENSG00000120952_ENSG00000229571_MG94_nonclock.sh 
sbatch -o ./Primate-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/ENSG00000120952_ENSG00000274764_MG94_nonclock.sh 
sbatch -o ./Primate-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/ENSG00000120952_ENSG00000280267_MG94_nonclock.sh 
sbatch -o ./Primate-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/ENSG00000204481_ENSG00000204505_MG94_nonclock.sh 
sbatch -o ./Primate-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/ENSG00000204481_ENSG00000229571_MG94_nonclock.sh 
sbatch -o ./Primate-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/ENSG00000204481_ENSG00000274764_MG94_nonclock.sh 
sbatch -o ./Primate-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/ENSG00000204481_ENSG00000280267_MG94_nonclock.sh 
