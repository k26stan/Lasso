p1 <- read.csv("1se-lasso-pred-1.csv")
p2 <- read.csv("1se-lasso-pred-2.csv")
p3 <- read.csv("1se-lasso-pred-3.csv")
p4 <- read.csv("1se-lasso-pred-4.csv")
p5 <- read.csv("1se-lasso-pred-5.csv")
p6 <- read.csv("1se-lasso-pred-6.csv")
p7 <- read.csv("1se-lasso-pred-7.csv")
p8 <- read.csv("1se-lasso-pred-8.csv")
p9 <- read.csv("1se-lasso-pred-9.csv")
p10 <- read.csv("1se-lasso-pred-10.csv")

p <- rbind(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10)

top <- data.frame(NA, cor(p[,1],p[,2])^2, cor(p[,1],p[,3])^2)
names(top) <- names(p)

out <- rbind(top,p)

names(out) <- c("observed","covariates","SNPs_covariates")
out[1,1] <- "r2"

write.csv(out, "1se-lasso-prediction.csv", quote=FALSE, row.names=FALSE)
q()

