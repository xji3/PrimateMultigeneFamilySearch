#!/bin/bash
sbatch -o ./Primate-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/Force_ENSG00000088832_ENSG00000198225_MG94_nonclock.sh 
sbatch -o ./Primate-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/Force_ENSG00000134590_ENSG00000212747_MG94_nonclock.sh 
sbatch -o ./Primate-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/Force_ENSG00000166840_ENSG00000255151_MG94_nonclock.sh 
sbatch -o ./Primate-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/Force_ENSG00000189186_ENSG00000226372_MG94_nonclock.sh 
sbatch -o ./Primate-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/Force_ENSG00000241294_ENSG00000243264_MG94_nonclock.sh 
