## Run Lasso Regression ##
## Based on Nathan's Script "lasso.r" ##
## January 27, 2015 ##
## Kristopher Standish ##

###################################################
## PARSE COMMAND LINE #############################
###################################################

LINE <- commandArgs(trailingOnly = TRUE)
# LINE <- c("/projects/janssen/Psych/Lasso/20150127_PHENO_NAMES/PsychChip_R092670-PSY-3006_R092670_arm.raw","/projects/janssen/Psych/Pheno/Full_Table.txt", "/projects/janssen/Pysch/Pheno/Full_Table.txt","/projects/janssen/Psych/Pheno/PHENO_NAMES.txt","BL_PANSS,AGE_DIAG,PC1,PC2","/projects/janssen/Psych/Lasso/20150127_PHENO_NAMES/",10)
# LINE <- c("/projects/janssen/Psych/Lasso/20150127_PHENO_NAMES/TEMP_1000.csv","/projects/janssen/Psych/Pheno/Full_Table.txt", "/projects/janssen/Psych/Pheno/Full_Table.txt","/projects/janssen/Psych/Pheno/PHENO_NAMES.txt","BL_PANSS,AGE_DIAG,PC1,PC2","/projects/janssen/Psych/Lasso/20150127_PHENO_NAMES/",10)
# LINE <- c("/Users/kstandis/Downloads/20150126b_PHENO_NAMES/TEMP_1000.csv","/Users/kstandis/Data/Burn_Psych/Data/20150125_Full_Table.txt", "/Users/kstandis/Data/Burn_Psych/Data/20150125_Full_Table.txt","/Users/kstandis/Downloads/20150126b_PHENO_NAMES/PHENO_NAMES.txt","BL_PANSS,AGE_DIAG,PC1,PC2","/Users/kstandis/Downloads/20150126b_PHENO_NAMES/",10)
PathToGT <- LINE[1]
PathToPheno <- LINE[2]
PathToCov <- LINE[3]
Pheno_Name_List <- LINE[4]
Cov_List <- LINE[5]
PathToSave <- LINE[6]
Num_Iter <- as.numeric( LINE[7] )

PathToUpdate <- paste(PathToSave,"Update.txt",sep="")

###################################################
## LOAD DATA ######################################
###################################################

## Load Library to Plot
library(gplots) # heatmap.2
library(glmnet)

## Load Genotype File
write(paste(Sys.time(),"- Loading GT Data"),PathToUpdate,append=T)
start_time <- proc.time()
# GT <- read.table( PathToGT, header=T, nrow=1 )
GT <- read.csv( PathToGT, header=T )
write(paste(Sys.time(),"- GT Data Loaded"),PathToUpdate,append=T)
print( proc.time()-start_time )

## Load Phenotype File
PH.l <- read.table( PathToPheno, header=T, sep="\t" )

## Load Covariate File
COV.l <- read.table( PathToCov, header=T, sep="\t" )

## Load Phenotype Names File
PHENO_NAMES <- as.character( read.table( Pheno_Name_List )[,1] )
N_PHENO <- length(PHENO_NAMES)

###################################################
## GET ORGANIZED ##################################
###################################################

## Pull out Covariate Data
 # Split Covariate List
Cov_List.sub <- gsub( "PC","PC_", Cov_List )
Cov_List.sp <- strsplit( Cov_List.sub, "," )[[1]]
cov_cols <- which( colnames(COV.l) %in% Cov_List.sp )
COV <- COV.l[ , c(1,cov_cols) ]

## Impute Missing Covariate Data
for ( c in 2:ncol(COV) ) {
    MN <- mean( COV[,c], na.rm=TRUE )
    MISSING <- which(is.na( COV[,c] ))
    if ( length(MISSING)>0 ) { COV[ MISSING, c ] <- MN }
}

## Pull out Phenotype Columns
pheno_cols <- which( colnames(PH.l) %in% PHENO_NAMES )
PH.2 <- PH.l[ , c(1,pheno_cols) ]
 # Impute Missing Phenotypes
for ( c in 2:ncol(PH.2) ) {
    MN <- mean( PH.2[,c], na.rm=TRUE )
    MISSING <- which(is.na( PH.2[,c] ))
    if ( length(MISSING)>0 ) { PH.2[ MISSING, c ] <- MN }
}

## Impute Missing Genotypes
WHICH_COLS <- which( apply( GT, 2, function(x) NA %in% x ) )
for ( c in WHICH_COLS ) {
    MN <- mean( GT[,c], na.rm=TRUE )
    MISSING <- which(is.na( GT[,c] ))
    if ( length(MISSING)>0 ) { GT[ MISSING, c ] <- MN }
}

## Merge Tables
write(paste(Sys.time(),"- Merging Tables"),PathToUpdate,append=T)
MG.1 <- merge( x=PH.2, y=COV, by="IID" )
MG.2 <- merge( x=MG.1, y=GT[,c(1,7:ncol(GT))], by.x="IID",by.y="FID" )
write(paste(Sys.time(),"- Tables Merged"),PathToUpdate,append=T)

###################################################
## START LOOPING ##################################
###################################################

## Permute Sample Order
perm_order <- sample( 1:nrow(y), replace=F )
num_test_samps <- floor( nrow(y)/Num_Iter )

## Set alpha parameters
 # alpha=1 -> LASSO
 # alpha=0 -> Ridge
alpha <- .5

## Specify Lambda Bins
L.bins <- c("BEST","SE","MIN","MAX")

## Set up Variables to Store Data/Results
BETA <- PRED <- LAMBDA <- list()

## Loop through Phenotypes
for ( p in 1:length(PHENO_NAMES) ) {
# for ( p in 2 ) {
	## Pull out Phenotype Data
	pheno <- PHENO_NAMES[p]
	write( paste(Sys.time(),"- ## Running phenotype:",p,"(",pheno,") of",length(PHENO_NAMES)),PathToUpdate,append=T )
	pheno_rm <- setdiff( PHENO_NAMES, pheno )
	pheno_cols_rm <- which( colnames(MG.2) %in% pheno_rm )
	MG <- MG.2[ , -pheno_cols_rm ]
	## Define Matrices
	y <- data.matrix( MG[,pheno] )
	cov_cols <- which( colnames(MG) %in% colnames(COV) )
	x.cov <- data.matrix( MG[,cov_cols[ 2:length(cov_cols) ] ] )
	rm_cols <- which( colnames(MG) %in% c( pheno,"IID") )
	x.all <- data.matrix( MG[,-rm_cols] )
	rownames(x.cov) <- rownames(x.all) <- MG$IID
	## Set up some Variables to store data
	LAMBDA[[pheno]] <- array( ,c(Num_Iter,8) ) ; colnames(LAMBDA[[pheno]]) <- c("Best_C","SE_C","Min_C","Max_C","Best_A","SE_A","Min_A","Max_A")
	BETA[[pheno]] <- list()
	BETA[[pheno]]$cov <- BETA[[pheno]]$all <- list() # array( ,c(ncol(x.all),Num_Iter) )
	PRED[[pheno]] <- list()
	PRED[[pheno]]$cov <- PRED[[pheno]]$all <- list() 
	# PRED[[pheno]]$MAX <- PRED[[pheno]]$MIN <- PRED[[pheno]]$BEST <- PRED[[pheno]]$SE <- list() # array( ,c(nrow(x.all.ts),Num_Iter) )
	## Loop Through Iterations
	for ( i in 1:Num_Iter ) {
		iter <- paste("I",i,sep="_")
		## Specify Training/Test Sets
		which_samps <- num_test_samps*(i-1)+1:num_test_samps
		test_set <- perm_order[ which_samps ]
		train_set <- perm_order[ -which_samps ]
		 # Pull out Data for each
		y.tr <- y[ train_set, ]
		x.cov.tr <- x.cov[ train_set, ]
		x.all.tr <- x.all[ train_set, ]
		y.ts <- y[ test_set, ]
		x.cov.ts <- x.cov[ test_set, ]
		x.all.ts <- x.all[ test_set, ]

		## Run Lasso Regression
		 # On Covariates Alone
		print("Fitting Covariates")
		fit.cov <- glmnet(x=x.cov.tr, y=y.tr, alpha=alpha)
		cv.fit.cov <- cv.glmnet(x=x.cov.tr, y=y.tr, alpha=alpha, nfolds=20)
		L.cov.best <- cv.fit.cov$lambda.min
		L.cov.se <- cv.fit.cov$lambda.1se
		L.cov.min <- min(fit.cov$lambda)
		L.cov.max <- max(fit.cov$lambda)
		B.cov <- data.matrix( coef(fit.cov, s=c(L.cov.best,L.cov.se,L.cov.min,L.cov.max) ) )
		pred.cov <- predict(fit.cov, newx=x.cov.ts, s=c(L.cov.best,L.cov.se,L.cov.min,L.cov.max) )
		colnames(B.cov) <- colnames(pred.cov) <- L.bins
		# pairs( data.frame( pred.cov.best, pred.cov.se, pred.cov.max, pred.cov.min ) )

		 # On Covariates + SNPs
		print("Fitting Covariates + SNPs")
		write(paste(Sys.time(),"- Fitting Covariates + SNPs"),PathToUpdate,append=T)
		fit.all <- glmnet(x=x.all.tr, y=y.tr, alpha=alpha)
		cv.fit.all <- cv.glmnet(x=x.all.tr, y=y.tr, alpha=alpha, nfolds=20)
		L.all.best <- cv.fit.all$lambda.min
		L.all.se <- cv.fit.all$lambda.1se
		L.all.min <- min(fit.all$lambda)
		L.all.max <- max(fit.all$lambda)
		B.all <- data.matrix( coef(fit.all, s=c(L.all.best,L.all.se,L.all.min,L.all.max) ) )
		pred.all <- predict(fit.all, newx=x.all.ts, s=c(L.all.best,L.all.se,L.all.min,L.all.max) )
		colnames(B.all) <- colnames(pred.all) <- L.bins
		# pairs( data.frame( pred.all.best, pred.all.se, pred.all.max, pred.all.min ) )

		## Filter out unused Beta values
		which_zero <- which( rowSums(B.all)==0 )
		B.all.2 <- B.all[ -which_zero, ]

		## Compile Outputs
		print("Compiling Data")
		LAMBDA[[pheno]][i,] <- c( L.cov.best, L.cov.se, L.cov.min, L.cov.max, L.all.best, L.all.se, L.all.min, L.all.max )
		BETA[[pheno]]$cov[[iter]] <- data.frame( B.cov )
		BETA[[pheno]]$all[[iter]] <- data.frame( B.all.2 )
		PRED[[pheno]]$cov[[iter]] <- data.frame( y.ts, pred.cov )
		PRED[[pheno]]$all[[iter]] <- data.frame( y.ts, pred.all )

		## Status Update
		print(paste( "Done with iter:",i,"of",Num_Iter) )
		write(paste(Sys.time(),"- Done with iter:",i,"of",Num_Iter),PathToUpdate,append=T)
	} # Close Iteration Loop
	print(paste( "## Done with phenotype:",p,"(",pheno,") of",length(PHENO_NAMES) ))
} # Close Phenotype Loop

write(paste(Sys.time(),"- Done fitting models"),PathToUpdate,append=T)

# fit1 <- fit.cov
# fit1 <- fit.all
# COLS.list <- c("firebrick2","chocolate2","gold2","springgreen2","steelblue2","slateblue3")
# COLS <- colorRampPalette(COLS.list)(nrow(fit1$beta))
# plot( 0,0,type="n", xlim=range(fit1$lambda), ylim=range(fit1$beta) )
# for ( l in 1:nrow(fit1$beta) ) {
# 	points( fit1$lambda, fit1$beta[l,], col=COLS[l], type="o", pch=20 )
# }

# n <- ncol(coef(fit.cov))
# plot( 0,0,type="n", xlim=c(1,1.1*n), ylim=range(coef(fit.cov)) )
# for ( i in 1:nrow(coef(fit.cov)) ) { points( 1:n, coef(fit.cov)[i,], col=i, type="o", pch=20 ) }
# text( n, coef(fit.cov)[,n], label=rownames(coef(fit.cov)), col=1:nrow(coef(fit.cov)), pos=4 )

# n <- ncol(coef(fit.all))
# plot( 0,0,type="n", xlim=c(1,1.1*n), ylim=range(coef(fit.all)) )
# for ( i in 1:nrow(coef(fit.all)) ) { points( 1:n, coef(fit.all)[i,], col=i, type="o", pch=20 ) }
# text( n, coef(fit.all)[,n], label=rownames(coef(fit.all)), col=1:nrow(coef(fit.all)), pos=4 )

###################################################
## WRITE DATA #####################################
###################################################

## Write Data
save( BETA, file=paste(PathToSave,"BETA.Rdata",sep="") )
save( PRED, file=paste(PathToSave,"PRED.Rdata",sep="") )
save( LAMBDA, file=paste(PathToSave,"LAMBDA.Rdata",sep="") )

write(paste(Sys.time(),"- Done saving Rdata"),PathToUpdate,append=T)

###################################################
## PLOT RESULTS ###################################
###################################################

## Compile Predicted Values into one big data frame
 # 9 columns: 1 actual, 4 cov pred, 4 cov+SNP pred
PRED.2 <- list()
for ( p in 1:length(PRED) ) {
	pheno <- names(PRED)[p]
	PRED.2[[pheno]] <- array( , c(0,9) )
	colnames(PRED.2[[pheno]]) <- c( "OBS",paste(L.bins,"C",sep="_"),paste(L.bins,"A",sep="_") )
	for ( i in 1:Num_Iter ) {
		iter <- paste("I",i,sep="_")
		TEMP <- data.frame( PRED[[pheno]]$cov[[iter]], PRED[[pheno]]$all[[iter]][,L.bins] )
 		colnames(TEMP) <- c( "OBS",paste(L.bins,"C",sep="_"),paste(L.bins,"A",sep="_") )
		PRED.2[[pheno]] <- rbind( PRED.2[[pheno]], TEMP )
	}
}

## Plot Predicted vs Actual Values
pheno <- PHENO_NAMES[2]
par(mfrow=c(2,4))
for ( c in 2:ncol(PRED.2[[pheno]]) ) {
	MOD <- lm( PRED.2[[pheno]][,1] ~ PRED.2[[pheno]][,c] )
	NAMES <- colnames( PRED.2[[pheno]] )
	plot( PRED.2[[pheno]][,1] ~ PRED.2[[pheno]][,c], pch="+", col="steelblue2", xlab=NAMES[c], ylab=NAMES[1], main="Predicted vs Observed Values" )
	abline( MOD, col="steelblue4" )
	P_VAL <- summary(MOD)$coefficients[2,4]
	COR <- cor( PRED.2[[pheno]][,1], PRED.2[[pheno]][,c] )
	if ( P_VAL < .05 ) { STAR <- "*" }else{ STAR <- "" }
	text( min(PRED.2[[pheno]][,c]), min(PRED.2[[pheno]][,1]), label=paste(STAR,"P:",formatC(P_VAL,format="e",digits=2)), pos=4 )
	text( min(PRED.2[[pheno]][,c]), min(PRED.2[[pheno]][,1])+4, label=paste("Corr:",round(COR,3)), pos=4 )
}

## Correlation b/n Observed & Predicted
CORS <- array( ,c(9,4) )
colnames(CORS) <- L.bins
rownames(CORS) <- colnames(PRED.2[[1]])
for ( p in 1:length(PHENO_NAMES) ) { CORS[,p] <- cor( PRED.2[[p]] )[,1] }
 # Plot it
COLS <- c("firebrick2","gold2","chartreuse2","deepskyblue2")
YLIM <- c(-1,1)
# plot( 0,0,type="n", xlim=XLIM, ylim=YLIM, main="Correlation Coefficient b/n Predicted & Observed", xlab="" )
barplot( t(CORS), beside=T, col=COLS, ylim=YLIM, main="Correlation Coefficient b/n Predicted & Observed", ylab="Correlation Coef", las=2 )
abline( h=seq(-1,1,.2), lty=2, col="grey50" )
barplot( t(CORS), beside=T, col=COLS, ylim=YLIM, main="Correlation Coefficient b/n Predicted & Observed", ylab="Correlation Coef", las=2, add=T )

## Merge Beta Values from different Iterations
BETA.2 <- list()
for ( p in 1:length(BETA) ) {
	pheno <- names(BETA)[p]
	BETA.2[[pheno]] <- list()
	for ( l in 1:length(L.bins) ) {
		lbin <- L.bins[l]
		BETA.2[[pheno]][[lbin]] <- list()
		BETA.2[[pheno]][[lbin]]$cov <- data.frame( Row.names=rownames(BETA[[pheno]]$cov$I_1), I_1=BETA[[pheno]]$cov$I_1[,lbin] )
		rownames(BETA.2[[pheno]][[lbin]]$cov) <- rownames( BETA[[pheno]]$cov$I_1 )
		BETA.2[[pheno]][[lbin]]$all <- data.frame( Row.names=rownames(BETA[[pheno]]$all$I_1), I_1=BETA[[pheno]]$all$I_1[,lbin] )
		rownames(BETA.2[[pheno]][[lbin]]$all) <- rownames( BETA[[pheno]]$all$I_1 )
		for ( i in 2:Num_Iter ) {
			iter <- paste("I",i,sep="_")
			TEMP <- data.frame( BETA[[pheno]]$cov[[iter]][,lbin] )
			rownames(TEMP) <- rownames(BETA[[pheno]]$cov[[iter]])
			colnames(TEMP) <- iter
			BETA.2[[pheno]][[lbin]]$cov <- merge( BETA.2[[pheno]][[lbin]]$cov, TEMP, by.x="Row.names", by.y="row.names", all=T ) # [,2:(i+2)]
			TEMP <- data.frame( BETA[[pheno]]$all[[iter]][,lbin] )
			rownames(TEMP) <- rownames(BETA[[pheno]]$all[[iter]])
			colnames(TEMP) <- iter
			BETA.2[[pheno]][[lbin]]$all <- merge( BETA.2[[pheno]][[lbin]]$all, TEMP, by.x="Row.names", by.y="row.names", all=T ) # [,2:(i+1)]
		}
		## Remove Beta Values that are 0 for all Iterations
		 # Cov
		which_zero <- which( rowSums(BETA.2[[pheno]][[lbin]]$all[2:(Num_Iter+1)],na.rm=T)==0 )
		if ( length(which_zero) > 0 ) {
			NAMES <- BETA.2[[pheno]][[lbin]]$cov[ -which_zero, "Row.names" ]
			BETA.2[[pheno]][[lbin]]$cov <- BETA.2[[pheno]][[lbin]]$cov[ -which_zero, 2:(Num_Iter+1) ]
		}else{
			NAMES <- BETA.2[[pheno]][[lbin]]$cov[ , "Row.names" ]
			BETA.2[[pheno]][[lbin]]$cov <- BETA.2[[pheno]][[lbin]]$cov[ , 2:(Num_Iter+1) ]
		}
		rownames(BETA.2[[pheno]][[lbin]]$cov) <- NAMES
		 # All
		which_zero <- which( rowSums(BETA.2[[pheno]][[lbin]]$all[2:(Num_Iter+1)],na.rm=T)==0 )
		if ( length(which_zero) > 0 ) {
			NAMES <- BETA.2[[pheno]][[lbin]]$all[ -which_zero, "Row.names" ]
			BETA.2[[pheno]][[lbin]]$all <- BETA.2[[pheno]][[lbin]]$all[ -which_zero, 2:(Num_Iter+1) ]
		}else{
			NAMES <- BETA.2[[pheno]][[lbin]]$all[ , "Row.names" ]
			BETA.2[[pheno]][[lbin]]$all <- BETA.2[[pheno]][[lbin]]$all[ , 2:(Num_Iter+1) ]
		}
		rownames(BETA.2[[pheno]][[lbin]]$all) <- NAMES
	}
}

## Plot Coefficients vs Iteration
COLS.list <- c("firebrick2","chocolate2","gold1","springgreen2","steelblue2","slateblue3")
COLS <- colorRampPalette(COLS.list)(Num_Iter)
for ( p in 1:length(BETA.2) ) {
	pheno <- names(BETA.2)[p]
	par(mfrow=c(2,4))
	for ( set in c("cov","all") ) {
		for ( l in 1:length(L.bins) ) {
			lbin <- L.bins[l]
			par(ask=T)
			# barplot( t(data.matrix( BETA.2[[pheno]][[lbin]][[set]] )), col=COLS, beside=T, main=paste(pheno,":",set,"-",lbin), las=2 )
			boxplot( t(data.matrix( BETA.2[[pheno]][[lbin]][[set]] )), col=COLS, beside=T, main=paste(pheno,":",set,"-",lbin), las=2 )
		}
	}
}





#################################################
## END OF DOC ###################################
#################################################

# library(glmnet)
# set <- commandArgs()[3]                                 ## 'set' was the dataset (e.g. CGIS-er_oros-european)
# X <- commandArgs()[4]                                   ## This is the cross validation number
# x <- as.numeric(X)

# ## READ IN DATA
# pheno <- read.table("CGIS-er_oros-european.pheno",F)
# covar <- read.table("CGIS-er_oros-european.cov",F)
# geno <- read.csv("genotyped-qc-imputed.csv")            ## This is a 1 GB file

# ## IMPUTE MISSING VARIABLES ##
# avg <- mean(pheno$V3, na.rm=TRUE)
# for (i in 1:dim(pheno)[1]) {
#     if (is.na(pheno$V3[i])) { pheno$V3[i] <- avg }    
# }
# for (j in 3:dim(covar)[2]) {
#     avg <- mean(covar[,j], na.rm=TRUE)
#     for (i in 1:dim(covar)[1]) {
#         if (is.na(covar[i,j])) { covar[i,j] <- avg }    
#     }
# }

# ## MERGE DATASETS ##
# p1 <- pheno[,2:dim(pheno)[2]]
# c1 <- covar[,2:dim(covar)[2]]
# p2 <- merge(p1, c1, by="V2")
# a1 <- merge(p2, geno, by.x="V2", by.y="ID")

# ## DEFINE MATRICES ##
# y <- data.matrix(a1[,2])
# xf <- data.matrix(a1[,3:dim(p2)[2]])
# xl <- data.matrix(a1[(dim(p2)[2]+1):dim(a1)[2]])

# ## SET DATASETS ##
# set.seed(12345678)
# perm <- sample(1:dim(y)[1],replace=FALSE)
# incr <- floor(dim(xl)[1]/10)
# if (x==10) {
# 	k <- perm[(((x-1)*incr)+1):dim(xl)[1]]
# }else{
# 	k <- perm[(((x-1)*incr)+1):(x*incr)]
# }

# ## TRAINING DATASET ##
# y.tr <- data.matrix(y[-k,])
# xf.tr <- data.matrix(xf[-k,])
# xl.tr <- data.matrix(xl[-k,])
# xall.tr <- data.matrix(data.frame(xf.tr, xl.tr))

# ## TESTING DATASET ##
# y.ts <- data.matrix(y[k,])
# xf.ts <- data.matrix(xf[k,])
# xl.ts <- data.matrix(xl[k,])
# xall.ts <- data.matrix(data.frame(xf.ts, xl.ts))

# ## LASSO - COVARIATES ##
# fit.cov <- glmnet(x=xf.tr, y=y.tr)
# cv.fit.cov <- cv.glmnet(x=xf.tr, y=y.tr)
# #m.cov <- coef(fit.cov, s=cv.fit.cov$lambda.min)
# m.cov <- coef(fit.cov, s=cv.fit.cov$lambda.1se)
# pred.cov <- m.cov[1] + (xf.ts %*% m.cov[2:dim(m.cov)[1]])

# ## LASSO - ALL ##
# fit.all <- glmnet(x=xall.tr, y=y.tr)
# cv.fit.all <- cv.glmnet(x=xall.tr, y=y.tr)
# #m.all <- coef(fit.all, s=cv.fit.all$lambda.min)
# m.all <- coef(fit.all, s=cv.fit.all$lambda.1se)
# pred.all <- m.all[1] + (xall.ts %*% m.all[2:dim(m.all)[1]])

# pred.out <- data.frame(y.ts, pred.cov, pred.all)

# names(pred.out) <- c("observed","covariates","covariates-snps")
# #write.csv(pred.out, paste("/gpfs/group/schork/nwineing/jnj/lasso/", set, "/lasso-pred-", X, ".csv", sep=""), quote=FALSE, row.names=FALSE)
# write.csv(pred.out, paste(set, "/1se-lasso-pred-", X, ".csv", sep=""), quote=FALSE, row.names=FALSE)
# q()



















