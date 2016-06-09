#!/bin/bash
sbatch -o ./OutFolder/Primate-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/ENSG00000170627_ENSG00000124196_MG94_nonclock.sh 
sbatch -o ./OutFolder/Primate-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/ENSG00000134107_ENSG00000123095_MG94_nonclock.sh 
