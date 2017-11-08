# Script to 
#1. Summarize the distribution of deltaA for mean and sd
#2. Take Chr8 data, calculate deltaA, standardize it
#3. Make a grid of proportions of sd's for B_co and B_cc
#4. Run Chr8 according to the sd values of the grid
#5. Make one plot showing cc, co, and BMA for each grid combination
# Created 11/01/2017

rm(list=ls())
#Load libraries and assign parameters/PATHs
library(data.table)
library(Rmpfr) #For extreme values in the BMA cf calculation
chr = 8
pop = "AFR"
LOC.path <- "C:/Users/Lilith Moss/Documents/MATH/RESEARCH/Admixture_Project/Real_Data/"
HPC.path <- "/home/pmd-01/chemenya/ADMIX/PAGE/"
setwd("C:/Users/Lilith Moss/Documents/MATH/RESEARCH/Admixture_Project/Simulation_Results/11.7.2017/")

####################
# READ IN DATA 
####################
#Global Ancestry
Q <- do.call(rbind,read.table(paste0(LOC.path,"ONCO-LAPC_Phase3GlobalRef_chr1-22_EM0_TrioPhased.0.ForwardBackward.",pop,".scoreSumByIndiv.txt.GlobalAncestry"), header=F, sep="")) #Raw Q is read as a list
#Local Ancestry
#Local is 12421(loci) x 2288(People)
Local <- fread(paste0(LOC.path,"ONCO-LAPC_Phase3GlobalRef_chr", chr, "_EM0_TrioPhased.0.ForwardBackward.txt.",pop))
d <- Local
#Outcome 
d.fam <- read.table(paste0(LOC.path,"OncoArray_cleaned_AAPC-AABC-LAPC_for_Impute_V1.fam.LAPC.ageStudy.indexID.txt"),header=T)
#Map File for plots
d.map <- read.table(paste(LOC.path,"ONCO-LAPC_phase3_global_ref.Chr", chr, ".GeneticPos.txt", sep=""), header=T, sep="")

######################
# STANDARDIZE deltaA
######################
#Define outcomes for models
Y <- d.fam[,3]-1
Y <- ifelse(Y==-10, NA, Y)
Z <- ifelse(Y==1, 0, ifelse(Y==0, 1, NA))
#Averages for plotting
avg.L.cases <- apply(d, 1, FUN=function(L) { mean(as.numeric(L)[Y==1], na.rm=T) })
avg.L.controls <- apply(d, 1, FUN=function(L) { mean(as.numeric(L)[Y==0], na.rm=T) })
#Take average local ancestry per locus for all people
L.ave <- as.matrix(apply(Local,1,mean))
#Take average of global ancestry for all individuals
Q.ave <- mean(Q)
#Take deltaA
dA <- (L.ave/2 - Q.ave)
#Get distribution of dA
mdA <- mean(dA)
sdA <- sd(dA)
#Standardize dA
stan.dA <- (dA-mdA)/sdA
#Mean deltaA in controls for prior or alpha
Local.matrix <- t(as.matrix(Local))
all.dat <- cbind(Y,Q,Local.matrix) #[Y Q L1 L2 L3 ... L12421]
controls.dat <- all.dat[Y==0,]
Local.controls.only <- t(controls.dat[,3:12423])
Q.controls.only <- controls.dat[,2]
L.ave.controls.only <- as.matrix(apply(Local.controls.only,1,mean))
Q.ave.controls.only <- mean(Q.controls.only)
dA.controls.only <- (L.ave.controls.only/2 - Q.ave.controls.only)
mdA.controls.only <- mean(dA.controls.only) #Mean deviation in controls only
sdA.controls.only <- sd(dA.controls.only)
stan.dA.controls <- (dA.controls.only-mdA.controls.only)/sdA.controls.only
#Mean deltaA in cases 
cases.dat <- all.dat[Y==1,]
Local.cases.only <- t(cases.dat[,3:12423])
Q.cases.only <- cases.dat[,2]
L.ave.cases.only <- as.matrix(apply(Local.cases.only,1,mean))
Q.ave.cases.only <- mean(Q.cases.only)
dA.cases.only <- (L.ave.cases.only/2 - Q.ave.cases.only)
mdA.cases.only <- mean(dA.cases.only) #Mean deviation in controls only
sdA.cases.only <- sd(dA.cases.only)
stan.dA.cases.only <- (dA.cases.only-mdA.cases.only)/sdA.cases.only

#############################
#Model standardized deltaA
#############################
#Define model parameters for CO, CC, and BMA
pmw=c(1,1)
pmw <- pmw/sum(pmw)
n<-length(Y) 
models <- rbind(c(0,1),c(1,1))
numModels <- 2

#Set Constant Parameters
sigma.1 <- 1
sigma.2 <- 1
#Set new grid parameters
alpha.SE.m <- c(0.05,0.10,0.50,0.75)
beta.SE.m <- c(0.05,0.10,0.50,0.75)
SE.grid <- data.frame(data.matrix(expand.grid(alpha.SE.m,beta.SE.m) ))
names(SE.grid) <- c("alpha.SE","beta.SE")

for(i in 1:nrow(SE.grid)){
  # j=0
  #Calculate Variance
  alpha.SE <- SE.grid[i,1]
  beta.SE <- SE.grid[i,2]
  #Run for every row of L (Takes about 12 minutes for 13k rows)
  # system.time(r <- apply(d, 1, FUN=function(L) {
  system.time(r <- apply(d[1:1000,], 1, FUN=function(L) {
    #Counter
    # j=j+1
    # print(j)
    #Run simple CO and CC models on standardized deltaA
    reg.raw <- as.list(rep(0, numModels))
    #Regression on transformed variable P which is standardized to have variance 1
    # P <- (Q-as.numeric((L/2)))/sdA
    P <- (as.numeric((L/2))-Q)/sdA
    #Case-Only
    reg.raw[[1]] <- lm(P ~ -1 + Y) #Case-only
    #Case-Control
    reg.raw[[2]] <- lm(P ~ 1 + Y)  #Case - Control
    #Back transform for interprability
    reg <- reg.raw
    # reg[[1]]$coefficients <- reg.raw[[1]]$coefficients*sdA #11.8.2017
    # reg[[2]]$coefficients <- reg.raw[[2]]$coefficients*sdA #11.8.2017
    ####################
    #Run BMA analysis
    ####################
    #Calculate BMA Result
    betas <- (unlist(lapply(reg, FUN=function(r) { summary(r)$coef["Y","Estimate"] })))
    # betas.se <- sdA*(unlist(lapply(reg, FUN=function(r) { summary(r)$coef["Y","Std. Error"] }))) #11.8.2017
    betas.se <- (unlist(lapply(reg, FUN=function(r) { summary(r)$coef["Y","Std. Error"] })))
    #Prior Covariance matrix for the case-control model
    V.cc <- matrix(c((alpha.SE/sigma.1),0,0,beta.SE/sigma.1),ncol=2,byrow=T)
    #Prior Variance of Case-only matrix using SE(Beta) for the intercept
    cov1 <- sigma.1*matrix(beta.SE/sigma.1) 
    #Prior Variance of Case-control matrix
    cov2 <- sigma.1*V.cc 
    #Prior covariance matricies for both CO and CC models
    Cov_0 <- list(cov1,cov2)
    #Invert the Prior covariance matrix to get prior precision matrix
    lambda_0 <- list(solve(Cov_0[[1]]),solve(Cov_0[[2]]) )
    #Specify Design Matrix (list of design matrices)
    X <- list( matrix(Y),matrix(c(rep(1,length(Y)),Y),ncol=2) )
    #Specify Posterior Precision Matrix from Prior precision for each model
    lambda_n <- list(t(X[[1]])%*%X[[1]]+lambda_0[[1]],t(X[[2]])%*%X[[2]]+lambda_0[[2]])
    #Specify prior mean vector for Betas and Alpha
    mu_0 <- list(matrix(0),matrix(c(mdA.controls.only,0)))
    #Specify posterior mean vector for betas
    mu_n <- list( solve(t(X[[1]])%*%X[[1]]+lambda_0[[1]]) * ( (t(X[[1]])%*%X[[1]]*betas[1]) + (lambda_0[[1]]*mu_0[[1]]) ),
                  solve(t(X[[2]])%*%X[[2]]+lambda_0[[2]]) %*% ( (t(X[[2]])%*%X[[2]] %*% matrix(reg[[2]]$coefficients)) + (lambda_0[[2]] %*% mu_0[[2]]) ) )
    b0 = 1
    a0 = 1/(sigma.1) + 1
    #Specify Posterior hyperparameters for sigma^2
    an <- a0+(n/2)
    # bn <- list( b0+(1/2)*(t(Q)%*%Q + t(mu_0[[1]])%*%lambda_0[[1]]%*%mu_0[[1]] - t(mu_n[[1]])%*%lambda_n[[1]]%*%mu_n[[1]]),
    #             b0+(1/2)*(t(Q)%*%Q + t(mu_0[[2]])%*%lambda_0[[2]]%*%mu_0[[2]] - t(mu_n[[2]])%*%lambda_n[[2]]%*%mu_n[[2]]) )
    bn <- list( b0+(1/2)*(t(P)%*%P + t(mu_0[[1]])%*%lambda_0[[1]]%*%mu_0[[1]] - t(mu_n[[1]])%*%lambda_n[[1]]%*%mu_n[[1]]),
                b0+(1/2)*(t(P)%*%P + t(mu_0[[2]])%*%lambda_0[[2]]%*%mu_0[[2]] - t(mu_n[[2]])%*%lambda_n[[2]]%*%mu_n[[2]]) )
    
    #Calculate large values using multiple precision package (Rmpfr)
    lterm1 <- exp(as(((-n/2)*log(2*pi)),"mpfr")) #1/(2pi)^(n/2)
    lterm2 <- list( exp(as(an*log(bn[[1]]),"mpfr")), exp(as(an*log(bn[[2]]),"mpfr")) )
    lterm3 <- gamma(as(an,"mpfr"))
    lterm4 <- gamma(as(a0,"mpfr"))
    #Calculate Marginal Likelihood
    PrDGivenM1 <-lterm1*sqrt(det(lambda_0[[1]])/det(lambda_n[[1]]))*((b0^a0)/lterm2[[1]])*(lterm3/lterm4) 
    PrDGivenM2 <-lterm1*sqrt(det(lambda_0[[2]])/det(lambda_n[[2]]))*((b0^a0)/lterm2[[2]])*(lterm3/lterm4) 
    PrDGivenM <- c(PrDGivenM1,PrDGivenM2)
    #Calculate Posterior Model Probabilities
    PrMGivenD.new <- PrDGivenM*pmw/sum( PrDGivenM*pmw )
    post.beta <- sum(betas*PrMGivenD.new)
    post.se <- sqrt(sum(((betas.se^2)+(betas^2))*PrMGivenD.new) - (post.beta^2))
    z.score <- post.beta/post.se
    p.value <- 2*pnorm(-abs(z.score))
    #BMA.cf <- as.numeric( c(post.beta, post.se, z.score, p.value))
    BMA.cf <- as.numeric( c(post.beta, post.se, z.score, p.value,
                            as.numeric(PrMGivenD.new[1]),as.numeric(PrMGivenD.new[2])))
    names(BMA.cf) <- c("post.beta", "post.se", "z.score", "p.value",
                       "PrDGivenM1","PrDGivenM2")
    #Back-transform cf
    BMA.cf[1:2] <- sdA*BMA.cf[1:2]
    #AIC Version
    ll <- unlist(lapply(reg, AIC))
    fitness <- ll-log(pmw)
    PrMGivenD.aic <- exp(-fitness+min(fitness))/sum(exp(-fitness+min(fitness)))
    betas <- unlist(lapply(reg, FUN=function(r) { summary(r)$coef["Y","Estimate"] }))
    betas.se <- sdA*(unlist(lapply(reg, FUN=function(r) { summary(r)$coef["Y","Std. Error"] })))
    post.beta.aic <- sum(betas*PrMGivenD.aic)
    post.se.aic <- sqrt(sum(((betas.se^2)+(betas^2))*PrMGivenD.aic) - (post.beta.aic^2))
    z.score.aic <- post.beta.aic/post.se.aic
    p.value.aic <- 2*pnorm(-abs(z.score.aic))
    BMA.aic <- c(post.beta.aic, post.se.aic, z.score.aic, p.value.aic,
                 PrMGivenD.aic[1],PrMGivenD.aic[2])
    names(BMA.aic) <- c("post.beta.aic", "post.se.aic", "z.score.aic", "p.value.aic",
                        "PrDGivenM1","PrDGivenM2")
    #Back-transform aic
    BMA.aic[1:2] <- sdA*BMA.aic[1:2]
    
    #Write Results
    p <- as.numeric((L/2)) - Q
    # Reg <- lm(Q ~ -1 + offset((as.numeric(L)/2))+Y+Z) 11.8.2017
    Reg <- lm(p ~ -1+Y+Z) 
    s.reg <- summary(Reg)$coef
    s.reg.co <- summary(reg[[1]])$coef
    s.reg.cc <- summary(reg[[2]])$coef
    s.reg.co[,1:2] <- sdA*s.reg.co[,1:2]
    s.reg.cc[,1:2] <- sdA*s.reg.cc[,1:2]
    write.table(t(c(s.reg[1,],s.reg[2,])), paste0("Chr", chr, "CaseOnlyControlOnly",pop,"_SE(alpha)_",alpha.SE,"_SE(beta)_",beta.SE,".txt"), append=T, quote=F, sep=" ", row.names=F, col.names=F)
    write.table(t(c(s.reg.cc[1,],s.reg.cc[2,])), paste0("Chr", chr, "CaseControl",pop,"_SE(alpha)_",alpha.SE,"_SE(beta)_",beta.SE,".txt", sep=""), append=T, quote=F, sep=" ", row.names=F, col.names=F)
    write.table(t(c(BMA.cf)), paste0("Chr", chr, "BMA.cfResults",pop,"_SE(alpha)_",alpha.SE,"_SE(beta)_",beta.SE,".txt"), append=T, quote=F, sep=" ", row.names=F, col.names=F)
    write.table(t(c(BMA.aic)), paste0("Chr", chr, "BMA.aicResults",pop,"_SE(alpha)_",alpha.SE,"_SE(beta)_",beta.SE,".txt"), append=T, quote=F, sep=" ", row.names=F, col.names=F) 
  })) #End of Apply
  
  ############################################################
  #Plot Results for this combination of SE(alpha) and SE(beta)
  ############################################################
  pdf(paste("Chr", chr, "_Results_",i,"_",pop,"_SE(alpha)_",alpha.SE,"_SE(beta)_",beta.SE,".pdf", sep=""), width=12, height=8)
    #Read in just calculated results
    r.cc <- as.data.frame(fread(paste0("Chr", chr, "CaseOnlyControlOnly",pop,"_SE(alpha)_",alpha.SE,"_SE(beta)_",beta.SE,".txt")))
    r.bma.cf <- as.data.frame(fread(paste0("Chr", chr, "BMA.cfResults",pop,"_SE(alpha)_",alpha.SE,"_SE(beta)_",beta.SE,".txt")))#[,-1]
    r.bma.aic <- as.data.frame(fread(paste0("Chr", chr, "BMA.aicResults",pop,"_SE(alpha)_",alpha.SE,"_SE(beta)_",beta.SE,".txt")))#[,-1]
    r.caco   <- as.data.frame(fread(paste0("Chr", chr, "CaseControl",pop,"_SE(alpha)_",alpha.SE,"_SE(beta)_",beta.SE,".txt")))
    
    #PLOT 4: Genetic position
    par(fig=c(0,1,0,.4))
    plot(d.map$pos/1000000, avg.L.cases, pch=16, cex=0.5, xlab="Position (MB)", ylab="Local Ancestry", col="blue", ylim=c(min(c(avg.L.cases,avg.L.controls)), max(c(avg.L.cases, avg.L.controls))))
    points(d.map$pos/1000000, avg.L.controls, pch=16, cex=0.5, col="green")
    abline(h=mean(Q[Y==1], na.rm=T)*2, col="blue", lwd=2)
    abline(h=mean(Q[Y==0], na.rm=T)*2, col="green", lwd=2)
    
    #PLOT 3: Posterior Probabilities between CO(m1) and CC(m2)
    par(fig=c(0,1,.4,.8), new=T)
    plot(d.map$pos/1000000, type="n",r.bma.cf[,5], pch=16, cex=0.5, ylim=c(0,1), xlab="", 	ylab="Posterior Probability", col="blue", xaxt="n")
    lines(d.map$pos/1000000, r.bma.cf[,6], pch=16, cex=0.5, ylim=c(0,1), xlab="", 	ylab="-log(p-value)", col="green", xaxt="n")
    lines(d.map$pos/1000000, r.bma.cf[,5], pch=16, cex=0.5, ylim=c(0,1), xlab="", 	ylab="-log(p-value)", col="blue", xaxt="n")
    abline(h=0.5, col="grey", lty=2)
    
    #PLOT 2: Effect sizes and SEs 
    par(fig=c(0,1,.2,.6), new=T)
    plot(d.map$pos/1000000, type="n", r.cc[,1], pch=16, cex=0.5, ylim=c(min(r.cc[,1],r.cc[,5]),max(r.cc[,1],r.cc[,5])), xlab="", 	ylab="Effect Size, SE", col="blue", xaxt="n")
    #Case-only
    lines(d.map$pos/1000000, r.cc[,1], pch=16, cex=0.5, ylim=c(min(r.cc[,1],r.cc[,5]),max(r.cc[,1],r.cc[,5])), xlab="", 	ylab="-log(p-value)", col="blue", xaxt="n")
    lines(d.map$pos/1000000, r.cc[,2], pch=16, lty=2,cex=0.25, ylim=c(min(r.cc[,1],r.cc[,5]),max(r.cc[,1],r.cc[,5])), xlab="", 	ylab="-log(p-value)", col="blue", xaxt="n")
    #Case-control
    lines(d.map$pos/1000000, r.cc[,5], pch=16, cex=0.5, ylim=c(min(r.cc[,1],r.cc[,5]),max(r.cc[,1],r.cc[,5])), xlab="", 	ylab="-log(p-value)", col="green", xaxt="n")
    lines(d.map$pos/1000000, r.cc[,6], pch=16, lty=2, cex=0.25, ylim=c(min(r.cc[,1],r.cc[,5]),max(r.cc[,1],r.cc[,5])), xlab="", 	ylab="-log(p-value)", col="green", xaxt="n")
    #BMA
    lines(d.map$pos/1000000, r.bma.cf[,1], pch=16, cex=0.5, ylim=c(min(r.cc[,1],r.cc[,5]),max(r.cc[,1],r.cc[,5])), xlab="", 	ylab="-log(p-value)", col="pink", xaxt="n")
    lines(d.map$pos/1000000, r.bma.cf[,2], pch=16, lty=2, cex=0.5, ylim=c(min(r.cc[,1],r.cc[,5]),max(r.cc[,1],r.cc[,5])), xlab="", 	ylab="-log(p-value)", col="pink", xaxt="n")
    abline(h=0, col="grey", lty=2)
    
    #PLOT 1:P-values for different methods
    par(fig=c(0,1,.6,1), new=T)
    #Case-Only
    plot(d.map$pos/1000000, -log10(r.cc[,4]), pch=16, cex=0.5, ylim=c(0,max(c(10, -log10(r.cc[,4])))), xlab="", 	ylab="-log(p-value)", col="blue", xaxt="n")
    #Case-Control
    points(d.map$pos/1000000, -log10(r.cc[,8]), pch=16, cex=0.5, ylim=c(0,max(c(10, -log10(r.cc[,8])))), xlab="", 	ylab="-log(p-value)", col="green", xaxt="n")
    #BMA.cf
    points(d.map$pos/1000000, -log10(r.bma.cf[,4]), pch=16, cex=0.5,ylim=c(0,max(c(10, -log10(r.bma.cf[,4])))), xlab="", 	ylab="-log(p-value)", col="pink", xaxt="n")
    # s.co <- r.cc[,8]<0.001
    # s.co1 <- r.bma.cf[,4]<0.001
    # points(d.map$pos[s.co]/1000000, rep(-.5, sum(s.co)), pch=15, cex=1.2, col="red")
    # points(d.map$pos[s.co1]/1000000, rep(-.5, sum(s.co1)), pch=15, cex=1.2, col="red")
    abline(h=5, col="grey", lty=2)
    abline(h=6, col="grey", lty=2)
    title(paste0("Pvalues: Posteriors: Betas & SE: Position  Chr:",chr," Pop:",pop," SE(alpha)=",alpha.SE," SE(Beta)=",beta.SE) )
  dev.off()
} #End of SE(alpha) SE(Beta) combination






reg.raw <- as.list(rep(0, numModels))
#re-define the regressions
P <- (Q-as.numeric((L/2)))/sdA
p <- (Q-as.numeric((L/2)))
#Case-Only
reg.raw[[1]] <- lm(P ~ -1 + Y) #Case-only
#Case-Control
reg.raw[[2]] <- lm(P ~ 1 + Y)  #Case - Control
#Back transform for interprability
reg <- reg.raw
reg[[1]]$coefficients <- reg.raw[[1]]$coefficients*sdA
reg[[2]]$coefficients <- reg.raw[[2]]$coefficients*sdA


tr.co <- summary(lm(P ~ -1 + Y)) #Case-only
tr.cc <- summary(lm(P ~ 1 + Y))  #Case - Control

Utr.co <- summary(lm(p ~ -1 + Y)) #Case-only
Utr.cc <- summary(lm(p ~ 1 + Y))  #Case - Control





