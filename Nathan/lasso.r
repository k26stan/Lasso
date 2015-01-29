library(glmnet)
set <- commandArgs()[3]                                 ## 'set' was the dataset (e.g. CGIS-er_oros-european)
X <- commandArgs()[4]                                   ## This is the cross validation number
x <- as.numeric(X)

## READ IN DATA
pheno <- read.table("CGIS-er_oros-european.pheno",F)
covar <- read.table("CGIS-er_oros-european.cov",F)
geno <- read.csv("genotyped-qc-imputed.csv")            ## This is a 1 GB file

## IMPUTE MISSING VARIABLES ##
avg <- mean(pheno$V3, na.rm=TRUE)
for (i in 1:dim(pheno)[1]) {
    if (is.na(pheno$V3[i])) { pheno$V3[i] <- avg }    
}
for (j in 3:dim(covar)[2]) {
    avg <- mean(covar[,j], na.rm=TRUE)
    for (i in 1:dim(covar)[1]) {
        if (is.na(covar[i,j])) { covar[i,j] <- avg }    
    }
}

## MERGE DATASETS ##
p1 <- pheno[,2:dim(pheno)[2]]
c1 <- covar[,2:dim(covar)[2]]
p2 <- merge(p1, c1, by="V2")
a1 <- merge(p2, geno, by.x="V2", by.y="ID")

## DEFINE MATRICES ##
y <- data.matrix(a1[,2])
xf <- data.matrix(a1[,3:dim(p2)[2]])
xl <- data.matrix(a1[(dim(p2)[2]+1):dim(a1)[2]])

## SET DATASETS ##
set.seed(12345678)
perm <- sample(1:dim(y)[1],replace=FALSE)
incr <- floor(dim(xl)[1]/10)
if (x==10) { k <- perm[(((x-1)*incr)+1):dim(xl)[1]] } else { k <- perm[(((x-1)*incr)+1):(x*incr)] }

## TRAINING DATASET ##
y.tr <- data.matrix(y[-k,])
xf.tr <- data.matrix(xf[-k,])
xl.tr <- data.matrix(xl[-k,])
xall.tr <- data.matrix(data.frame(xf.tr, xl.tr))

## TESTING DATASET ##
y.ts <- data.matrix(y[k,])
xf.ts <- data.matrix(xf[k,])
xl.ts <- data.matrix(xl[k,])
xall.ts <- data.matrix(data.frame(xf.ts, xl.ts))

## LASSO - COVARIATES ##
fit.cov <- glmnet(x=xf.tr, y=y.tr)
cv.fit.cov <- cv.glmnet(x=xf.tr, y=y.tr)
#m.cov <- coef(fit.cov, s=cv.fit.cov$lambda.min)
m.cov <- coef(fit.cov, s=cv.fit.cov$lambda.1se)
pred.cov <- m.cov[1] + (xf.ts %*% m.cov[2:dim(m.cov)[1]])

## LASSO - ALL ##
fit.all <- glmnet(x=xall.tr, y=y.tr)
cv.fit.all <- cv.glmnet(x=xall.tr, y=y.tr)
#m.all <- coef(fit.all, s=cv.fit.all$lambda.min)
m.all <- coef(fit.all, s=cv.fit.all$lambda.1se)
pred.all <- m.all[1] + (xall.ts %*% m.all[2:dim(m.all)[1]])

pred.out <- data.frame(y.ts, pred.cov, pred.all)

names(pred.out) <- c("observed","covariates","covariates-snps")
#write.csv(pred.out, paste("/gpfs/group/schork/nwineing/jnj/lasso/", set, "/lasso-pred-", X, ".csv", sep=""), quote=FALSE, row.names=FALSE)
write.csv(pred.out, paste(set, "/1se-lasso-pred-", X, ".csv", sep=""), quote=FALSE, row.names=FALSE)
q()


