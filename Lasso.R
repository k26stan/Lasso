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
# LINE <- c("/Users/kstandis/Downloads/20150127_Lasso/TEMP_1000.csv","/Users/kstandis/Data/Burn_Psych/Data/20150125_Full_Table.txt", "/Users/kstandis/Data/Burn_Psych/Data/20150125_Full_Table.txt","/Users/kstandis/Downloads/20150126b_PHENO_NAMES/PHENO_NAMES.txt","BL_PANSS,AGE_DIAG,PC1,PC2","/Users/kstandis/Downloads/20150126b_PHENO_NAMES/",10)
# LINE <- c("/projects/janssen/Psych/Lasso/20150817_PHENO_NAMES/Merged_CND_DEL_p01.Both.csv","/projects/janssen/Psych/Pheno/Full_Table.txt", "/projects/janssen/Psych/Pheno/Full_Table.txt","/projects/janssen/Psych/Pheno/PHENO_NAMES.txt","BL_PANSS,AGE_DIAG,PC1,PC2","/projects/janssen/Psych/Lasso/20150817_PHENO_NAMES/",10)
PathToGT <- LINE[1]
PathToPheno <- LINE[2]
PathToCov <- LINE[3]
Pheno_Name_List <- LINE[4]
Cov_List <- LINE[5]
PathToSave <- LINE[6]
Num_Iter <- as.numeric( LINE[7] )

## Path To Bim File
TEMP <- unlist(strsplit( PathToGT, "/" ))
TEMP <- paste( TEMP[-(length(TEMP)-1)], collapse="/" )
PathToGTKey <- gsub( "csv","bim", TEMP )

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
print(paste( "Genotypes Loaded:", (proc.time()-start_time)[3] ))
BIM <- read.table( PathToGTKey, header=F, sep="\t" )
colnames(BIM) <- c("CHR","SNP","X","POS","REF","ALT")

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
perm_order <- sample( 1:nrow(MG.2), replace=F )
num_test_samps <- floor( nrow(MG.2)/Num_Iter )

## Set alpha parameters
 # alpha=1 -> LASSO
 # alpha=0 -> Ridge
alpha <- .8

## Specify Lambda Bins
L.bins <- c("BEST","SE","MIN","MAX")

## Set up Variables to Store Data/Results
BETA <- PRED <- LAMBDA <- L <- CV <- list()

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
	L[[pheno]] <- list()
	L[[pheno]]$cov <- L[[pheno]]$all <- list()
	CV[[pheno]] <- list()
	CV[[pheno]]$cov <- CV[[pheno]]$all <- list()
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
		L[[pheno]]$all[[iter]] <- fit.all
		L[[pheno]]$cov[[iter]] <- fit.cov
		CV[[pheno]]$all[[iter]] <- cv.fit.all
		CV[[pheno]]$cov[[iter]] <- cv.fit.cov
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
save( L, file=paste(PathToSave,"L.Rdata",sep="") )
save( CV, file=paste(PathToSave,"CV.Rdata",sep="") )
load( file=paste(PathToSave,"BETA.Rdata",sep="") )
load( file=paste(PathToSave,"PRED.Rdata",sep="") )
load( file=paste(PathToSave,"LAMBDA.Rdata",sep="") )
load( file=paste(PathToSave,"L.Rdata",sep="") )
load( file=paste(PathToSave,"CV.Rdata",sep="") )

write(paste(Sys.time(),"- Done saving Rdata"),PathToUpdate,append=T)

# load( "Data/Burn_Psych/Plots/20150817_Lasso/Alpha.8/L.Rdata" )
# load( "Data/Burn_Psych/Plots/20150817_Lasso/Alpha.8/LAMBDA.Rdata" )
# load( "Data/Burn_Psych/Plots/20150817_Lasso/Alpha.8/BETA.Rdata" )
# load( "Data/Burn_Psych/Plots/20150817_Lasso/Alpha.8/PRED.Rdata" )
# load( "Data/Burn_Psych/Plots/20150817_Lasso/Alpha.8/CV.Rdata" )

###################################################
## PLOT RESULTS ###################################
###################################################

## Compile Predicted Values into one big data frame
 # 9 columns: 1 observed, 4 cov-only pred, 4 cov+SNP pred
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
COLS <- c("firebrick2","gold2","chartreuse2","deepskyblue2")
for ( p in 1:length(PRED) ) {
	pheno <- names(PRED)[p]
	LIM <- range( unlist(lapply( PRED[[pheno]]$all, range )) )
	jpeg( paste(PathToSave,"1-Pred_v_Obs.",pheno,".jpeg",sep=""), height=1200,width=2400,pointsize=30 )
	par(mfrow=c(2,4))
	for ( c in 2:ncol(PRED.2[[pheno]]) ) {
		MOD <- lm( PRED.2[[pheno]][,1] ~ PRED.2[[pheno]][,c] )
		NAMES <- colnames( PRED.2[[pheno]] )
		plot( PRED.2[[pheno]][,1] ~ PRED.2[[pheno]][,c], pch="+", col=COLS[p], xlab=NAMES[c], ylab=NAMES[1], main=paste("Predicted vs Observed -",pheno),xlim=LIM,ylim=LIM )
		abline( 0,1, col="grey50",lty=3,lwd=1 )
		abline( MOD, col=gsub("2","4",COLS[p]),lty=2,lwd=2 )
		P_VAL <- summary(MOD)$coefficients[2,4]
		COR <- cor( PRED.2[[pheno]][,1], PRED.2[[pheno]][,c] )
		if ( P_VAL < .05 ) { STAR <- "*" }else{ STAR <- "" }
		# text( quantile(LIM,0), quantile(LIM,.95), label=paste(STAR,"P:",formatC(P_VAL,format="e",digits=2)), pos=4 )
		text( quantile(LIM,0), quantile(LIM,.9), label=paste("Corr:",round(COR,3)), pos=4 )
	}
	dev.off()
}

## Correlation b/n Observed & Predicted
CORS <- array( ,c(9,4) )
colnames(CORS) <- names(PRED.2)
rownames(CORS) <- colnames(PRED.2[[1]])
for ( p in 1:length(PRED.2) ) { CORS[,p] <- cor( PRED.2[[p]] )[,1] }
EXPL <- CORS^2
 # Plot it
COLS <- c("firebrick2","gold2","chartreuse2","deepskyblue2")
YLIM <- c(-1,1)
jpeg( paste(PathToSave,"1b-Pred_v_Obs.COR.jpeg",sep=""), height=1800,width=2400,pointsize=30 )
par(mfrow=c(2,1))
 # Correlation Coefficient (Obs vs Pred)
barplot( t(CORS), beside=T, col=COLS, ylim=YLIM, main="Correlation Coefficient b/n Predicted & Observed", ylab="Correlation Coef (R)", las=2 )
abline( h=seq(-1,1,.2), lty=2, col="grey50" )
barplot( t(CORS), beside=T, col=COLS, ylim=YLIM, main="Correlation Coefficient b/n Predicted & Observed", ylab="Correlation Coef (R)", las=2, add=T )
legend( "bottomleft", fill=COLS, legend=colnames(CORS) )
 # Variance Explained (Obs vs Pred)
YLIM <- c(0,1)
barplot( t(EXPL), beside=T, col=COLS, ylim=YLIM, main="Variance Explained by Model", ylab="Correlation Coef (R^2)", las=2 )
abline( h=seq(-1,1,.2), lty=2, col="grey50" )
barplot( t(EXPL), beside=T, col=COLS, ylim=YLIM, main="Variance Explained by Model", ylab="Correlation Coef (R^2)", las=2, add=T )
# legend( "bottomleft", fill=COLS, legend=colnames(CORS) )
dev.off()

## Plot Lambda Values for each Model
COLS <- c("firebrick2","gold2","chartreuse2","deepskyblue2")
jpeg( paste(PathToSave,"2-LAMBDA_Boxplot.jpeg",sep=""), height=1800,width=2400,pointsize=30 )
par(mfrow=c(2,2))
for ( p in 1:length(BETA) ) {
	pheno <- names(BETA)[p]
	boxplot( log10(LAMBDA[[pheno]]), col=COLS[p], las=2, main=paste("Lambda vs Model Fit -",pheno),ylab="log10(Lambda) (per model in 10xCV)" )
}
dev.off()

# ## Distribution of Response Phenotype
#  # Show Distribution of Raw Phenotype & Residuals
# BIN.SIZES <- c(.5,1,2,5,10)
# par(mfrow=c(1,4))
# for ( p in 1:length(PRED.2) ) {
# 	pheno <- names(PRED.2)[p]
# 	MOD <- lm( PRED.2[[pheno]][,"OBS"] ~ PRED.2[[pheno]][,"SE_A"] )
# 	XLIM <- range( PRED.2[[pheno]], resid(MOD) )
# 	BIN.SIZE <- BIN.SIZES[ which.min( abs(diff(XLIM)/20 - BIN.SIZES) ) ]
# 	BRKS <- seq(floor(XLIM[1]),ceiling(XLIM[2])+BIN.SIZE,BIN.SIZE)
# 	YLIM <- c( 0, max( hist( PRED.2[[pheno]]$OBS,breaks=BRKS,plot=F )$counts, hist( resid(MOD),breaks=BRKS,plot=F )$counts ) )
# 	# Observed Hist
# 	hist( PRED.2[[pheno]]$OBS, xlim=XLIM,ylim=YLIM,breaks=BRKS,col="grey50",density=20,angle=45,xlab=paste("Patient",pheno),main="Phenotypic Distribution" )
# 	# hist( PRED.2[[pheno]]$SE_A, xlim=XLIM,ylim=YLIM,breaks=BRKS,col=COLS[p],density=30,angle=-45,xlab=paste("Patient",pheno),main="Phenotypic Distribution", add=T )
# 	hist( resid(MOD), xlim=XLIM,ylim=YLIM,breaks=BRKS,col=COLS[p],density=30,angle=-45,xlab=paste("Patient",pheno),main="Phenotypic Distribution", add=T )
# }

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
		for ( set in c("cov","all") ) {
			which_zero <- which( rowSums(BETA.2[[pheno]][[lbin]]$all[2:(Num_Iter+1)],na.rm=T)==0 )
			if ( length(which_zero) > 0 ) {
				NAMES <- BETA.2[[pheno]][[lbin]][[set]][ -which_zero, "Row.names" ]
				BETA.2[[pheno]][[lbin]][[set]] <- BETA.2[[pheno]][[lbin]][[set]][ -which_zero, 2:(Num_Iter+1) ]
			}else{
				NAMES <- BETA.2[[pheno]][[lbin]][[set]][ , "Row.names" ]
				BETA.2[[pheno]][[lbin]][[set]] <- BETA.2[[pheno]][[lbin]][[set]][ , 2:(Num_Iter+1) ]
			}
			rownames(BETA.2[[pheno]][[lbin]][[set]]) <- NAMES
			## Sort by ?RowSums?
			BETA.2[[pheno]][[lbin]][[set]][which(is.na(BETA.2[[pheno]][[lbin]][[set]]),arr.ind=T)] <- 0
			ORDER <- order( abs(rowSums(BETA.2[[pheno]][[lbin]][[set]])), decreasing=T )
			BETA.2[[pheno]][[lbin]][[set]] <- BETA.2[[pheno]][[lbin]][[set]][ORDER,]
		}
	}
}
# tail( BETA.2$DEL_PANSS$SE$all )

## Plot Map of Genome w/ SNPs Included in Model
COLS <- c("firebrick2","gold2","chartreuse2","deepskyblue2")
for ( p in 1:length(BETA) ) {
	pheno <- names(BETA)[p]
	jpeg( paste(PathToSave,"3a-Param_Genome.",pheno,".jpeg",sep=""), height=2400,width=3000,pointsize=30 )
	par(mfrow=c(2,2))
	for ( l in 1:length(L.bins) ) {
		lbin <- L.bins[l]
		plot( BIM$CHR, BIM$POS, pch=20,col="black", main=paste("Variants Included in:",lbin,"-",pheno),xlab="Chromosome",ylab="Position",xaxt="n" )
		axis( 1, at=1:24, labels=c(1:22,"X","Y") )
		for ( i in 1:Num_Iter ) {
			iter <- paste("I",i,sep="_")
			WHICH <- which( BETA.2[[pheno]][[lbin]]$all[[iter]] != 0 )
			WHICH <- WHICH[order(abs(BETA.2[[pheno]][[lbin]]$all[[iter]][WHICH]),decreasing=T)]
			WHICH.SNP <- sapply( strsplit( rownames(BETA.2[[pheno]][[lbin]]$all)[WHICH], "_" ), "[",1 )
			COLS.ramp <- c( colorRampPalette(c(gsub("2","4",COLS[p]),gsub("2","1",COLS[p])))(20),colorRampPalette(c(gsub("2","1",COLS[p]),"white"))(1.2*length(WHICH)) )[1:length(WHICH)]
			WHICH.SNP.IN <- match( WHICH.SNP, BIM$SNP ) # which( BIM$SNP %in% WHICH.SNP )
			points( .05*i+BIM$CHR[WHICH.SNP.IN], BIM$POS[WHICH.SNP.IN], pch="+",col=COLS.ramp )
		}
		points( BIM$CHR, BIM$POS, pch=20,col="black" )
	}
	dev.off()
}

## Boxplot Number of Predictors for MIN/MAX/BEST/SE
 # Compile # Parameters
N_PARAM <- array( ,c(2*Num_Iter*length(BETA.2)*length(L.bins),5) )
colnames(N_PARAM) <- c("Pheno","Lambda","Set","N_Param","Dev")
row <- 1
for ( p in 1:length(BETA.2) ) {
	pheno <- names(BETA.2)[p]
	for ( set in c("cov","all") ) {
		for ( l in 1:length(L.bins) ) {
			lbin <- L.bins[l]
			for ( i in 1:Num_Iter ) {
				iter <- paste("I",i,sep="_")
				N_PARAM[row,"Pheno"] <- pheno
				N_PARAM[row,"Lambda"] <- lbin
				N_PARAM[row,"Set"] <- set
				N_PARAM[row,"N_Param"] <- length(which( BETA.2[[pheno]][[lbin]][[set]][[iter]]!=0 ))
				N_PARAM[row,"Dev"] <- "Later"
				row <- row+1
			}
		}
	}
}
 # Boxplot of # Predictors in Model
jpeg( paste(PathToSave,"3-Num_Param.jpeg",sep=""), height=1200,width=1800,pointsize=30 )
par(mar=c(8,5,5,2))
# par(mfrow=c(2,1)) ; for ( set in c("cov","all") ) {
for ( set in "all" ) {
	YMAX <- max( as.numeric(as.character( N_PARAM[which(N_PARAM[,"Set"]==set),"N_Param"] )) )
	boxplot( as.numeric(as.character(N_Param)) ~ Pheno+Lambda, data=N_PARAM, subset=N_PARAM[,"Set"]==set, 
		main=paste("# Parameters per Model Fit -",set),ylab="# Parameter (per model in 10xCV)", col=COLS,las=2 )
	abline( h=seq(0,YMAX+200,200),lty=3,col="grey50",lwd=1 )
	boxplot( as.numeric(as.character(N_Param)) ~ Pheno+Lambda, data=N_PARAM, subset=N_PARAM[,"Set"]==set, 
		main=paste("# Parameters per Model Fit -",set),ylab="# Parameter (per model in 10xCV)", col=COLS,las=2,add=T )
}
legend( "bottomleft", fill=COLS, legend=unique(N_PARAM[,"Pheno"]) )
dev.off()

## Number of Predictors Included in Model
 # & Deviance Explained
COLS.list <- c("firebrick2","chocolate2","gold1","springgreen2","steelblue2","slateblue3")
COLS.iter <- colorRampPalette(COLS.list)(Num_Iter)
jpeg( paste(PathToSave,"3-Param_v_Lambda.jpeg",sep=""), height=1500,width=3000,pointsize=36 )
par(mar=c(8,5,5,5))
par(mfrow=c(2,4))
 # Penalty & # Parameters vs Fit
for ( p in 1:length(BETA.2) ) {
	pheno <- names(BETA.2)[p]
	XLIM <- log10(range( unlist(lapply( L[[pheno]]$all, function(x) x$lambda )) ))
	YLIM <- c( 0, 1 )
	YLIM.2 <- c( 0, max( unlist(lapply( L[[pheno]]$all, function(x) x$df )) ) )
	plot( 0,0,type="n", xlim=XLIM,ylim=YLIM, main=paste("Model Fits vs Penalty -",pheno),ylab="% Deviance",xlab="log10(Lambda) (Penalization)",yaxt="n")
	for ( i in 1:Num_Iter ) {
		iter <- paste("I",i,sep="_")
		points( log10(L[[pheno]]$all[[iter]]$lambda), L[[pheno]]$all[[iter]]$dev.ratio, type="l",col=COLS.iter[i] )
		points( log10(L[[pheno]]$all[[iter]]$lambda), L[[pheno]]$all[[iter]]$df / YLIM.2[2], type="l",col=COLS.iter[i],lty=2 )
		axis( 2, at=seq(0,1,.2), las=2)
		axis( 4, at=seq(0,1,.1), label=round(seq(0,YLIM.2[2],length.out=11)), las=2 )
		mtext("# Parameters", side=4, line=3, cex=.7,las=0)
		abline( v=log10(LAMBDA[[pheno]][i,grep("_A",colnames(LAMBDA[[pheno]]))]), col=COLS.iter[i],lty=1:4 )
		if ( p==1 ) {
			legend( quantile(XLIM,.3),YLIM[2], lty=1:4,title="(Vertical)",legend=grep("_A",colnames(LAMBDA[[pheno]]),value=T),cex=.8)
			legend( quantile(XLIM,.65),YLIM[2], lty=c(1,2),title="(Trace)",legend=c("% Dev","# Par"),cex=.8)
		}
	}
}
 # Cross-Validation Fits (CVfit)
for ( p in 1:length(BETA.2) ) {
	pheno <- names(BETA.2)[p]
	XLIM <- log10(range( unlist(lapply( CV[[pheno]]$all, function(x) x$lambda )) ))
	YLIM <- c( min(unlist(lapply( CV[[pheno]]$all, function(x) x$cvlo ))), max(unlist(lapply( CV[[pheno]]$all, function(x) x$cvup ))) )
	plot( 0,0,type="n", xlim=XLIM,ylim=YLIM, main=paste("Model Fits vs Penalty -",pheno),ylab=CV[[pheno]]$all[[iter]]$name,xlab="log10(Lambda) (Penalization)")
	for ( i in 1:Num_Iter ) {
		iter <- paste("I",i,sep="_")
		points( log10(CV[[pheno]]$all[[iter]]$lambda), CV[[pheno]]$all[[iter]]$cvm, type="l",col=COLS.iter[i] )
		arrows( log10(CV[[pheno]]$all[[iter]]$lambda), CV[[pheno]]$all[[iter]]$cvlo, log10(CV[[pheno]]$all[[iter]]$lambda), CV[[pheno]]$all[[iter]]$cvup, code=3,angle=90,length=.1,col=COLS.iter[i] )
		abline( v=log10(LAMBDA[[pheno]][i,grep("_A",colnames(LAMBDA[[pheno]]))]), col=COLS.iter[i],lty=1:4 )
	}
}
dev.off()

## Plot Coefficients vs Iteration
Num_Par <- 30
# COLS.list <- c("firebrick2","chocolate2","gold1","springgreen2","steelblue2","slateblue3")
# COLS <- colorRampPalette(COLS.list)(Num_Par)
COLS <- c("firebrick2","gold2","chartreuse2","deepskyblue2")
for ( p in 1:length(BETA.2) ) {
	pheno <- names(BETA.2)[p]
	jpeg( paste(PathToSave,"4-BETA_Box.",pheno,".jpeg",sep=""), height=1200,width=2400,pointsize=36 )
	par(mfrow=c(2,4))
	for ( set in c("cov","all") ) {
		for ( l in 1:length(L.bins) ) {
			lbin <- L.bins[l]
			if ( Num_Par=="all" ) {
				Num_Par.iter <- nrow( BETA.2[[pheno]][[lbin]][[set]] )
			}else{ Num_Par.iter <- min( Num_Par, nrow(BETA.2[[pheno]][[lbin]][[set]]) ) }
			# par(ask=T)
			boxplot( t(data.matrix( BETA.2[[pheno]][[lbin]][[set]][1:Num_Par.iter,] )), col=COLS[p], beside=T, main=paste(pheno,":",set,"-",lbin),ylab="Beta Estimate", las=2 )
		}
	}
	dev.off()
}

## Heatmap of Beta Values?
HEAT_COLS <- c( colorRampPalette(COLS.list[1:3])(50),"black",colorRampPalette(COLS.list[4:6])(50))
HEAT_BRKS <- c( seq(-1,0,length.out=51)[-51],c(-1,1)*1e-9,seq(0,1,length.out=51)[-1] )
Num_Par <- "all"
Num_Par <- 201
for ( p in 1:length(BETA.2) ) {
	pheno <- names(BETA.2)[p]
	set="all"
	for ( l in 1:length(L.bins) ) {
	# for ( l in grep("SE",L.bins) ) {
		lbin <- L.bins[l]
		TEMP <- BETA.2[[pheno]][[lbin]][[set]]
		TEMP[which(is.na(TEMP),arr.ind=T)] <- 0
		if ( Num_Par=="all" ) {
			Num_Par.iter <- nrow( TEMP )
		}else{ Num_Par.iter <- min( Num_Par, nrow(TEMP) ) }
		COL_COLS <- colorRampPalette(c("white","magenta3"))(Num_Par.iter)
		if ( nrow(TEMP)>1 ) {
			jpeg( paste(PathToSave,"4-BETA_Heat.",pheno,".",lbin,".",Num_Par,".jpeg",sep=""), height=1200,width=2400,pointsize=30 )
			# heatmap.2( t(data.matrix( TEMP[1:Num_Par.iter,] )), col=HEAT_COLS, breaks=HEAT_BRKS,
				# scale="none",trace="none",ColSideColors=COL_COLS,main=paste("Beta Values -",lbin) )
			heatmap.2( t(data.matrix( TEMP[1:Num_Par.iter,] )), col=HEAT_COLS, breaks=HEAT_BRKS,
				scale="none",trace="none",ColSideColors=COL_COLS,main=paste("Beta Values -",pheno,"-",lbin),
				)
				# Colv=F,dendrogram="row" )
			dev.off()
		}
	}
}

## Model Averaging
 # Take simple mean of Beta Coefficients to create "Averaged Model"
 # Then calculate R2 for "Averaged Model" while iteratively including additional terms
BETA.3 <- lapply( BETA.2, function(x) lapply( x, function(y) lapply( y, function(z) rowMeans(z) )))
MG.3 <- data.frame( MG.2, 1 )
colnames(MG.3)[ncol(MG.3)] <- "(Intercept)"
MG.3 <- MG.3[ which(MG.3$IID %in% rownames(PRED.2[[1]]) ), ]
Num_Par <- 300
# barplot( BETA.3$DEL_PANSS$SE$all )
PRED.3 <- list()
for ( p in 1:length(BETA.3) ) {
	pheno <- names(BETA.3)[p]
	PRED.3[[pheno]] <- list()
	# PRED.3[[pheno]] <- PRED.2[[pheno]][,-(2:5)]
	# colnames(PRED.3[[pheno]]) <- c("OBS",L.bins)
	for ( l in 1:length(L.bins) ) {
		lbin <- L.bins[l]
		# Pull Columns from MG.3
		# Pull Betas...iteratively...shit...redo some of this...
		# How many Parameters?
		if ( Num_Par=="all" ) {
			Num_Par.iter <- nrow( BETA.3[[pheno]][[lbin]]$all )
		}else{ Num_Par.iter <- min( Num_Par, nrow(data.frame(BETA.3[[pheno]][[lbin]]$all)) ) }
		PRED.3[[pheno]][[lbin]] <- array(,c(nrow(PRED.2[[pheno]]),Num_Par.iter+1) ) # list()
		colnames(PRED.3[[pheno]][[lbin]]) <- c("OBS",paste("Mod",1:Num_Par.iter,sep="_"))
		PRED.3[[pheno]][[lbin]][,"OBS"] <- MG.3[,pheno] # PRED.2[[pheno]][,"OBS"]
		for ( n in 1:Num_Par.iter ) {
			B <- BETA.3[[pheno]][[lbin]]$all[1:n]
			which_vars <- names(B) ; dim(B) <- c(n,1) ; rownames(B) <- which_vars
			X <- c(MG.3[,which_vars],recursive=T) ; dim(X) <- c(length(X)/n,n)
			# X <- data.frame( MG.3[,which_vars] )
			PRED.3[[pheno]][[lbin]][,n+1] <- c( X %*% B )
		}
	}
	print(paste("Done with",pheno))
}
PRED.3.res <- lapply( PRED.3, function(x) lapply( x, function(y) cor(y)[,"OBS"] ))
lapply( PRED.3.res, function(x) lapply( x, function(y) tail(y,1) ))

## Plot Correlation Coefficient of Averaged Model
XLIM <- c(0,Num_Par)
YLIM <- c(-1,1)
COLS <- c("firebrick2","gold2","chartreuse2","deepskyblue2")
jpeg( paste(PathToSave,"5-Pred_AvgMod.jpeg",sep=""), height=1200,width=1600,pointsize=30 )
plot( 0,0,type="n", xlim=XLIM,ylim=YLIM,xlab="# Parameters",ylab="Correlation Coefficient (Obs vs Pred)",main="Quality of Averaged Model" )
abline( h=seq(-1,1,.2),lty=3,col="grey50" )
abline( v=seq(0,Num_Par,25),lty=3,col="grey50" )
for ( p in 1:length(BETA.3) ) {
	pheno <- names(BETA.3)[p]
	PRED.3[[pheno]] <- list()
	for ( l in 1:length(L.bins) ) {
		lbin <- L.bins[l]
		x.temp <- PRED.3.res[[pheno]][[lbin]][-1]
		points( 1:length(x.temp), x.temp, col=COLS[p],lty=l,lwd=3,type="l" )
	}
}
legend("bottomleft",legend=names(BETA.3),fill=COLS)
legend("bottomright",legend=L.bins, lty=1:4,lwd=3 )
dev.off()


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



















