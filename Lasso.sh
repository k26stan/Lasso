## Run Lasso Regression on Janssen Psych Cohorts ##
## Based on Nathan's scripts ##
## January 27, 2015 ##
## Kristopher Standish ##

## Gameplan
#1-Specify Files/Paths/etc...
#2-Convert Bed to Raw/012 file
#4-Run Lasso? (loop within R script...only have to load genotypes once!)

##########################################################################
## 1 ## Set up Paths #####################################################
##########################################################################
 # Use Bash
 # Take in arguments, set up directories/paths for files/tools
echo \### 1 - `date` \###
echo \### Define Set Variables and Paths \###

###########################################################
## Manually Input Parameters ##

## Names/Paths for Output
DATE=$1
HOME_DIR=$2

## Parameters and Files
VAR_FILE=$3
VAR_DIR=$4
PHENO_DIR=$5
PHENO_FILE=$6
PHENO_NAME_LIST=$7 # Which Phenotype Files are you using?
COV_FILE=$8 # Path to Covariate File or "F"
COVS=$9 # Which Covariates to Include?
PC_COUNT=${10} # How many PCs to Include as Covariates?
START_STEP=${11} # Which Step do you want to start on?

###########################################################
## Constant Paths ##

## Public Tools
PLINK=/projects/janssen/Tools/plink_linux_x86_64/plink

## Custom Scripts
LASSO_R=/projects/janssen/Psych/Scripts/Lasso/Lasso.R

## Set Specific Paths
VAR_PATH=${VAR_DIR}/${VAR_FILE}
COV_PATH=${PHENO_DIR}/${COV_FILE}
PHENO_PATH=${PHENO_DIR}/${PHENO_FILE}

## Make new folder for Today's adventures
OUT_DIR=${HOME_DIR}/${DATE}_${PHENO_NAME_LIST%%.txt}
mkdir ${OUT_DIR}
cd ${OUT_DIR}

###########################################################
## Pull some Info out of Parameters ##

## Determine if I'm using Covariates
if [ -e ${COV_PATH} ]
then
USE_COVARS=TRUE
else
USE_COVARS=FALSE
fi

## Specify list of Covariates to include (for command and for filename)
if [[ $USE_COVARS == TRUE ]]
then

if [ $PC_COUNT -eq 0 ]
then
COVS_COMMAND=`echo "${COVS}" | sed 's/QQQ/,/g'`
COVS_FILENAME=`echo "${COVS}" | sed 's/QQQ/_/g'`
else
PCS=`seq 1 ${PC_COUNT}`
PCS_COMMAND=`echo "PC"${PCS} | sed 's/ /QQQPC/g'`
COVS_COMMAND=`echo "${COVS}QQQ${PCS_COMMAND}" | sed 's/QQQ/,/g'`
COVS_FILENAME=`echo "${COVS}QQQ${PCS_COMMAND}" | sed 's/QQQ/_/g'`
fi

## Incorporate Country/Site of Study as Binary Covariate (if Included)
if [[ $COVS == *COUN* ]]
then
COVS_COMMAND=`echo $COVS_COMMAND | sed 's/COUN/CN_ARG,CN_AUS,CN_COL,CN_HUN,CN_LTU,CN_MEX,CN_MYS,CN_NZL,CN_POL,CN_RUS,CN_UKR/g'`
fi

fi # Close (if USE_COVARS)
## Specify a File to which to Write Updates
UPDATE_FILE=${OUT_DIR}/Update.txt

## Done
if [ "$START_STEP" -le 1 ]; then
echo `date` "1 - Define Set Variables and Paths - DONE" > ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## 2 ## Convert Variant File #############################################
##########################################################################
if [ "$START_STEP" -le 2 ]; then
echo \### 2 - `date` \###
echo \### Convert Variant File \###
echo `date` "2 - Convert Variant File" >> ${UPDATE_FILE}

## Convert .bed file to .raw (012) file
if [ ${VAR_PATH: -4} == ".bed" ] ; then
${PLINK} \
--bfile ${VAR_PATH%%.bed} \
--silent \
--hardy midp \
--recode A \
--memory 8000 \
--threads 1 \
--out ${OUT_DIR}/${VAR_FILE%%.bed}
fi

## Convert .raw file to .csv file
sed 's/ /,/g' ${OUT_DIR}/${VAR_FILE%%.bed}.raw > ${OUT_DIR}/${VAR_FILE%%.bed}.csv

## Done
# if [ "$START_STEP" -le 2 ]; then
echo `date` "2 - Convert Variant File - DONE" > ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## 3 ## Run Lasso Regression #############################################
##########################################################################
if [ "$START_STEP" -le 3 ]; then
echo \### 3 - `date` \###
echo \### Run Lasso Regression \###
echo `date` "3 - Lasso Regression" >> ${UPDATE_FILE}

NUM_ITER=10
## Run Lasso Regression
Rscript ${LASSO_R} \
${OUT_DIR}/${VAR_FILE%%.bed}.csv \
${PHENO_PATH} \
${COV_PATH} \
${PHENO_DIR}/${PHENO_NAME_LIST} \
${COVS_COMMAND} \
${OUT_DIR} \
${NUM_ITER}


## Done
echo `date` "6 - Make Plots - DONE" > ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## END OF DOC ############################################################
##########################################################################












