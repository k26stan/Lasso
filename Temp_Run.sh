RUN_Lasso

################################################################################
## NEW GAMEPLAN (AUG 20, 2015) #################################################
################################################################################
## LASSO Cross-Validation Pipeline Outline
 # Split Cohort into C groups
 # Leaving 1 group (Test Set) out (for each of C groups)
   # Run Single Locus on Remaining (Training) Set
   # Build Lasso Model
   # Assess on Test Set
 # Weight Models by Correlation of Test Set
 # Find Weighted Average of Models
   # Iteratively Building 1 Paramater at a time (?)
   # Or Try Thresholding at Given BETA Value to limit # Terms

################################################################################
################################################################################
################################################################################

## Names/Paths
echo \### Defining Set Variables and Paths at `date` \###
DATE=20150820
HOME_DIR=/projects/janssen/Psych/Lasso
cd ${HOME_DIR}

## Files
# VAR_FILE=Merged_CND_DEL_CGIS_FULL_BL_CGIS_AGE_DIAG.bed
VAR_FILE=Merged_CND_DEL_p01.Both.bed
VAR_DIR=/projects/janssen/Psych/Lasso
PHENO_DIR=/projects/janssen/Psych/Pheno
PHENO_FILE=Full_Table.txt
PHENO_NAME_LIST=PHENO_NAMES.txt
COV_FILE=Full_Table.txt
COVS=`echo BL_PANSS AGE_DIAG`
PC_COUNT=2
START_STEP=1

COVS=`echo "$COVS" | sed 's/ /QQQ/g'`

########################################
## Run The Script
/projects/janssen/Psych/Scripts/Lasso/Lasso.sh \
${DATE} \
${HOME_DIR} \
${VAR_FILE} \
${VAR_DIR} \
${PHENO_DIR} \
${PHENO_FILE} \
${PHENO_NAME_LIST} \
${COV_FILE} \
${COVS} \
${PC_COUNT} \
${START_STEP}

################################################################################
################################################################################
################################################################################
################################################################################

## Names/Paths
echo \### Defining Set Variables and Paths at `date` \###
DATE=20150127
HOME_DIR=/projects/janssen/Psych/Lasso
cd ${HOME_DIR}

## Files
VAR_FILE=PsychChip_R092670-PSY-3006_R092670_arm.bed
VAR_DIR=/projects/janssen/Psych/Data/Genotyped/
PHENO_DIR=/projects/janssen/Psych/Pheno
PHENO_FILE=Full_Table.txt
PHENO_NAME_LIST=PHENO_NAMES.txt
COV_FILE=Full_Table.txt
COVS=`echo BL_PANSS AGE_DIAG`
PC_COUNT=2
START_STEP=1

COVS=`echo "$COVS" | sed 's/ /QQQ/g'`

########################################
## Run The Script
/projects/janssen/Psych/Scripts/Lasso/Lasso.sh \
${DATE} \
${HOME_DIR} \
${VAR_FILE} \
${VAR_DIR} \
${PHENO_DIR} \
${PHENO_FILE} \
${PHENO_NAME_LIST} \
${COV_FILE} \
${COVS} \
${PC_COUNT} \
${START_STEP}
