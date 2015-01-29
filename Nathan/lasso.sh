#!/bin/bash
PHENOARRAY[0]=CGIS-combined-european
#PHENOARRAY[1]=CGIS-er_oros-european
#PHENOARRAY[2]=CGIS-palmitate-european
#PHENOARRAY[3]=nPANSS-combined-european
#PHENOARRAY[4]=nPANSS-er_oros-european
#PHENOARRAY[5]=nPANSS-palmitate-european
#PHENOARRAY[6]=pPANSS-combined-european
#PHENOARRAY[7]=pPANSS-er_oros-european
#PHENOARRAY[8]=pPANSS-palmitate-european
#PHENOARRAY[9]=PANSS-combined-european
#PHENOARRAY[10]=PANSS-er_oros-european
#PHENOARRAY[11]=PANSS-palmitate-european
#for j in {0..11}
for j in {0..0}
do
pheno=${PHENOARRAY[${j}]}
#mkdir /gpfs/group/schork/nwineing/jnj/lasso/${pheno}
for i in {1..10}
do
qsub \
-v i=${i},pheno=${pheno} \
-N cv-${i}-${pheno} \
-o /gpfs/group/schork/nwineing/jnj/lasso/temp/${pheno}-${i}.o \
-e /gpfs/group/schork/nwineing/jnj/lasso/temp/${pheno}-${i}.e \
/gpfs/group/schork/nwineing/jnj/lasso/lasso.qsub

done # Close N-fold Cross-Validation Loop

done # Close Cohort/Phenotype Loop


## NOTES ON SCRIPT ##
# Loops through desired cohorts/phenotypes
# Loops through 10 iterations
# Submits new Qsub job for each iteration of each cohort/phenotype
  # by calling "lasso.qsub" and specifying job parameters