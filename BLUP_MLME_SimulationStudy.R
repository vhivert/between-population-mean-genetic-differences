library(parallel)
date=format(Sys.time(), "%e_%b_%Y")
date=gsub(" ", "", date)#remove whitespace for 1 digit dates

sim<- function(Fst=0.1){
  p0   <- sample(px,size=1e5,prob = wgts,replace = T) # runif(n=M,min=0.01,max=0.99)
  p1   <- rbeta(1e5,shape1=p0*(1-Fst)/Fst,shape2=(1-Fst)*(1-p0)/Fst)
  idx=sample(which(p1>0.01 & p1<0.99),size = M,replace = F) # Maf 1%
  p1 <-p1[idx]
  p2   <- rbeta(M,shape1=p0[idx]*(1-Fst)/Fst,shape2=(1-Fst)*(1-p0[idx])/Fst)
  
  X1  <- do.call("cbind",lapply(1:M,function(j) rbinom(n=N,size=2,prob=p1[j])))
  x1  <- do.call("cbind",lapply(1:M,function(j) rbinom(n=n,size=2,prob=p1[j])))
  x2  <- do.call("cbind",lapply(1:M,function(j) rbinom(n=n,size=2,prob=p2[j])))
  
  Z1 <- do.call("cbind",lapply(1:M,function(j) (X1[,j]-2*p1[j])/sqrt(2*p1[j]*(1-p1[j]))  ))
  
  ## Sim neutral
  B_true <- rnorm(n=M,mean=0,sd=sqrt(h2/(sum(2*p1*(1-p1)))))
  trueMean1 <- sum(2*p1*B_true)
  trueMean2 <- sum(2*p2*B_true)
  
  G1     <- c(X1%*%B_true) # Pop1 discovery genetic values
  Y1     <- rnorm(n=N,mean=G1,sd=sqrt(1-h2)) # Discovery phenotype
  g1     <- c(x1%*%B_true) # Pop1 target genetic values
  y1     <- rnorm(n=n,mean=g1,sd=sqrt(1-h2)) # Target phenotype
  g2     <- c(x2%*%B_true) # Pop2 target genetic values
  
  lambda <- M*(1-h2)/(h2)
  ones <- array(1, c(N)) # a vector of ones
  
  # build coefficient matrix by blocks
  coeff <- array(0, c(M + 1, M + 1))
  coeff[1:1, 1:1] <- t(ones) %*% ones
  coeff[1:1, 2:(M+1)] <- t(ones) %*% Z1
  coeff[2: 2:(M+1), 1] <- t(Z1) %*% ones
  coeff[2:(M+1), 2:(M+1)] <- t(Z1) %*% Z1 + I_M * lambda
  
  # build the right hand side
  rhs = rbind(t(ones) %*% Y1, t(Z1) %*% Y1)
  
  # get BLUP solution
  solution <- solve(coeff, rhs)
  mean_BLUE <- solution[1]
  B_blup <- solution[-1]/ sqrt(2*p1*(1-p1)) #per allele
  
  #PGS <- c(X%*%B_blup)
  pgs1 <- c(x1%*%B_blup)
  pgs2 <- c(x2%*%B_blup)
  
  predMean1 <- sum(2*p1*B_blup)
  predMean2 <- sum(2*p2*B_blup)
  
  slopePGS1  <- cov(y1,pgs1)/var(pgs1)
  slopePGS2  <- cov(g2,pgs2)/var(pgs2)
  slopeBeta  <- as.numeric(cov(B_true,B_blup)/var(B_blup))
  results    <- c(slopePGS1=slopePGS1,
                  slopePGS2=slopePGS2,
                  varPGS1=var(pgs1),
                  varPGS2=var(pgs2),
                  Rsq1=cor(g1,pgs1)^2,
                  Rsq2=cor(g2,pgs2)^2,
                  slopeBeta=slopeBeta,
                  trueMean1=trueMean1,predMean1=predMean1,PGSmean1=mean(pgs1),
                  trueMean2=trueMean2,predMean2=predMean2,PGSmean2=mean(pgs2),
                  estimated_mean1=mean_BLUE)
  return(results)
}

sim_NeutralBLUP <- function(Fst=0.1){
  p0   <- sample(px,size=1e5,prob = wgts,replace = T) # runif(n=M,min=0.01,max=0.99)
  p1   <- rbeta(1e5,shape1=p0*(1-Fst)/Fst,shape2=(1-Fst)*(1-p0)/Fst)
  idx=sample(which(p1>0.01 & p1<0.99),size = M,replace = F) # Maf 1%
  p1 <-p1[idx]
  p2   <- rbeta(M,shape1=p0[idx]*(1-Fst)/Fst,shape2=(1-Fst)*(1-p0[idx])/Fst)
  
  X1  <- do.call("cbind",lapply(1:M,function(j) rbinom(n=N,size=2,prob=p1[j])))
  x1  <- do.call("cbind",lapply(1:M,function(j) rbinom(n=n,size=2,prob=p1[j])))
  x2  <- do.call("cbind",lapply(1:M,function(j) rbinom(n=n,size=2,prob=p2[j])))
  
  ## Sim neutral
  B_true <- rnorm(n=M,mean=0,sd=sqrt(h2/(sum(2*p1*(1-p1)))))
  
  trueMean1 <- sum(2*p1*B_true)
  trueMean2 <- sum(2*p2*B_true)
  
  G1     <- c(X1%*%B_true) # Pop1 Discovery genetic values
  Y1     <- rnorm(n=N,mean=G1,sd=sqrt(1-h2)) # Pop1 Discovery phenotypes
  
  # print(var(G1))
  # print(var(G1/var(Y1)))
  g1     <- c(x1%*%B_true) # Pop1 Target genetic values
  y1     <- rnorm(n=n,mean=g1,sd=sqrt(1-h2)) # phenotypes
  g2     <- c(x2%*%B_true) # Pop2 Target genetic values
  
  lambda <- (1-h2)/(h2/(sum(2*p1*(1-p1))))
  
  ones <- array(1, c(N)) # a vector of ones
  
  # build coefficient matrix by blocks
  coeff <- array(0, c(M + 1, M + 1))
  coeff[1:1, 1:1] <- t(ones) %*% ones
  coeff[1:1, 2:(M+1)] <- t(ones) %*% X1
  coeff[2: 2:(M+1), 1] <- t(X1) %*% ones
  coeff[2:(M+1), 2:(M+1)] <- t(X1) %*% X1 + I_M * lambda
  
  # build the right hand side
  rhs = rbind(t(ones) %*% Y1, t(X1) %*% Y1)
  
  # get BLUP solution
  solution <- solve(coeff, rhs)
  mean_BLUE <- solution[1]
  B_blup <- solution[-1] #per allele
  
  pgs1 <- c(x1%*%B_blup)
  pgs2 <- c(x2%*%B_blup)
  
  predMean1 <- sum(2*p1*B_blup)
  predMean2 <- sum(2*p2*B_blup)
  
  slopePGS1  <- cov(g1,pgs1)/var(pgs1)
  slopePGS2  <- cov(g2,pgs2)/var(pgs2)
  slopeBeta  <- as.numeric(cov(B_true,B_blup)/var(B_blup))
  results    <- c(slopePGS1=slopePGS1,
                  slopePGS2=slopePGS2,
                  varPGS1=var(pgs1),
                  varPGS2=var(pgs2),
                  Rsq1=cor(g1,pgs1)^2,
                  Rsq2=cor(g2,pgs2)^2,
                  slopeBeta=slopeBeta,
                  trueMean1=trueMean1,predMean1=predMean1,PGSmean1=mean(pgs1),
                  trueMean2=trueMean2,predMean2=predMean2,PGSmean2=mean(pgs2),
                  estimated_mean1=mean_BLUE)
  return(results)
}

nrep=500

#1000 causals
N    <- 10000
n    <- 1000
M    <- 1000
I_M  <- diag(M)
h2   <- 0.5
px   <- seq(0.01,0.99,by=0.001)
wgts <- 1/(px*(1-px)); wgts <- wgts/sum(wgts)
Fst=c(0.001,0.01,0.05,0.1,0.2)

res=array(dim=c(nrep,14,length(Fst)))
for(i in 1:length(Fst)){
  res[,,i] <- do.call("rbind",mclapply(1:nrep,function(k) sim(Fst=Fst[i]), mc.cores=15))
}

res_neutralBLUP=array(dim=c(nrep,14,length(Fst)))
for(i in 1:length(Fst)){
  res_neutralBLUP[,,i] <- do.call("rbind",mclapply(1:nrep,function(k) sim_NeutralBLUP(Fst=Fst[i]), mc.cores=15))
}

#500 causals
N    <- 10000
n    <- 1000
M    <- 500
I_M  <- diag(M)

res_500=array(dim=c(nrep,14,length(Fst)))
for(i in 1:length(Fst)){
  res_500[,,i] <- do.call("rbind",mclapply(1:nrep,function(k) sim(Fst=Fst[i]), mc.cores=15))
}


#100 causals
N    <- 10000
n    <- 1000
M    <- 100
I_M  <- diag(M)
h2   <- 0.5

res_100=array(dim=c(nrep,14,length(Fst)))
for(i in 1:length(Fst)){
  res_100[,,i] <- do.call("rbind",mclapply(1:nrep,function(k) sim(Fst=Fst[i]), mc.cores=15))
}

#10 causals
N    <- 10000
n    <- 1000
M    <- 10
I_M  <- diag(M)

res_10=array(dim=c(nrep,14,length(Fst)))
for(i in 1:length(Fst)){
  res_10[,,i] <- do.call("rbind",mclapply(1:nrep,function(k) sim(Fst=Fst[i]), mc.cores=15))
}

save.image(file = paste0("BLUP_simulations_",date,".RData"))

######################################################################################################################
path.to.project="./"

## Supplementary Figure 12
FigureS12<-function(save=T){
  load(paste0("~/Documents/Projects/NewProjects/BtwPopInference/documents/Article/Data/BLUP_simulations_",date,".RData"))
  data=list(M1000=res,M1000_neut=res_neutralBLUP,M500=res_500,M100=res_100,M10=res_10)
  Fst=c(0.001,0.01,0.05,0.1,0.2)
  col=c("black","grey","blue","green","red")
  if(save){png(paste0(path.to.project,"FigureS12.png"),width = 12,height = 8.3,units = "in",res = 600)}
  layout(matrix(c(1,2,3),1,3,byrow = TRUE),widths = c(2.5,2.5,0.5),heights = c(2.5),TRUE)
  par(mar=c(5.1, 5.1, 4.1, 2.1))
  
  tmp=do.call(cbind,lapply(data, function(x) x[,1,]))
  l=seq(1,25,by=5)
  tmp=cbind(tmp[,l],tmp[,l+1],tmp[,l+2],tmp[,l+3],tmp[,l+4])
  all.ticks <- seq(1,ncol(tmp))
  my.labels <- Fst
  ticks=c(2.5,7.5,12.5,17.5,22.5)
  
  boxplot(tmp,
          xaxt = "n",
          ylab = bquote("slope"~G[1]~"~"~PGS[1]),
          xlab=bquote(F[ST]),
          border = "white",
          cex.lab = 2,
          cex.axis=1.5,
          at = all.ticks)
  boxplot(tmp,
          xaxt = "n",
          yaxt = "n",
          col = col,
          at = all.ticks,
          lwd = 0.5,
          outline = FALSE,
          add = TRUE)
  abline(h=1,lty=2,lwd=2)
  
  rect(ytop = 2,ybottom = -2,xleft = 5.5,xright = 10.5,col = rgb(0.5,0.5,0.5,alpha = 0.25),border = NA)
  rect(ytop = 2,ybottom = -2,xleft = 15.5,xright = 20.5,col = rgb(0.5,0.5,0.5,alpha = 0.25),border = NA)
  axis(1,at = seq(2.5,25,by=5),labels = my.labels)
  mtext("A.", side=3, adj=0, line=1.2, cex=2, font=2);
  
  tmp=do.call(cbind,lapply(data, function(x) x[,2,]))
  l=seq(1,25,by=5)
  tmp=cbind(tmp[,l],tmp[,l+1],tmp[,l+2],tmp[,l+3],tmp[,l+4])
  
  boxplot(tmp,
          xaxt = "n",
          ylab = bquote("slope"~G[2]~"~"~PGS[2]),
          xlab=bquote(F[ST]),
          border = "white",
          cex.lab = 2,
          cex.axis=1.5,
          at = all.ticks)
  boxplot(tmp,
          xaxt = "n",
          yaxt = "n",
          col = col,
          at = all.ticks,
          lwd = 0.5,
          outline = FALSE,
          add = TRUE)
  abline(h=1,lty=2,lwd=2)
  
  rect(ytop = 2,ybottom = -2,xleft = 5.5,xright = 10.5,col = rgb(0.5,0.5,0.5,alpha = 0.25),border = NA)
  rect(ytop = 2,ybottom = -2,xleft = 15.5,xright = 20.5,col = rgb(0.5,0.5,0.5,alpha = 0.25),border = NA)
  axis(1,at = seq(2.5,25,by=5),labels = my.labels)
  mtext("B.", side=3, adj=0, line=1.2, cex=2, font=2);
  
  par(mar = c(0,0,0,0))
  plot.new()
  legend(x = "center",title = "N SNPs",legend = c("1K","1K Neut","500","100","10"),fill = col)
  if(save){dev.off()}
}
FigureS12(save = T)

## Supplementary Figure 13
FigureS13<-function(save=T){
  load(paste0("~/Documents/Projects/NewProjects/BtwPopInference/documents/Article/Data/BLUP_simulations_",date,".RData"))
  data=list(M1000=res,M1000_neut=res_neutralBLUP,M500=res_500,M100=res_100,M10=res_10)
  Fst=c(0.001,0.01,0.05,0.1,0.2)
  
  if(save){png(paste0(path.to.project,"FigureS13.png"),width = 16,height = 12,units = "in",res = 600)}
  par(mfrow=c(2,5));par(mar=c(5.1,5.1,4.1,2.1))
  
  tmp=do.call(cbind,lapply(data, function(x) x[,8,]))
  xlim=c(min(tmp),max(tmp))
  tmp=do.call(cbind,lapply(data, function(x) x[,9,]))
  ylim=c(min(tmp),max(tmp))
  
  for(f in 1:5){
    #True vs Predicted mean diff
    true=do.call(cbind,lapply(data, function(x) x[,8,f]))
    est=do.call(cbind,lapply(data, function(x) x[,9,f]))
    slope=sapply(1:5,function(x){model=lm(true[,x]~est[,x]);return(coef(model)[2])})
    se=sapply(1:5,function(x){model=lm(true[,x]~est[,x]);return(summary(model)$coefficients[2,2])})
    eq=paste0("(",round(slope,digits=3),"+-",round(se,digits=3),")")
    col=rep(c("black","grey","blue","green","red"),each=nrep)
    plot(x=true,y=est,lty=2,pch=16,ylab=bquote(bar(PGS[1])==sum(2*p[1]*beta[BLUP],1,M)),xlab=bquote(bar(g[1])==sum(2*p[1]*beta,1,M)),col=col,ylim=ylim,xlim=xlim)
    mtext(bquote(F[ST]==.(Fst[f])), side=3, adj=0, line=1.2, cex=1, font=2)
    abline(0,1)
    legend(x = "bottomright",title = "N SNPs",legend = paste(c("1K","1K Neut","500","100","10"),eq),col = c("black","grey","blue","green","red"),pch = 16,bty = "n")
  }
  tmp=do.call(cbind,lapply(data, function(x) x[,11,]))
  xlim=c(min(tmp),max(tmp))
  tmp=do.call(cbind,lapply(data, function(x) x[,12,]))
  ylim=c(min(tmp),max(tmp))
  for(f in 1:5){
    #True vs Predicted mean diff
    true=do.call(cbind,lapply(data, function(x) x[,11,f]))
    est=do.call(cbind,lapply(data, function(x) x[,12,f]))
    slope=sapply(1:5,function(x){model=lm(true[,x]~est[,x]);return(coef(model)[2])})
    se=sapply(1:5,function(x){model=lm(true[,x]~est[,x]);return(summary(model)$coefficients[2,2])})
    eq=paste0("(",round(slope,digits=3),"+-",round(se,digits=3),")")
    col=rep(c("black","grey","blue","green","red"),each=nrep)
    plot(x=true,y=est,lty=2,pch=16,ylab=bquote(bar(PGS[2])==sum(2*p[2]*beta[BLUP],1,M)),xlab=bquote(bar(g[2])==sum(2*p[2]*beta,1,M)),col=col,ylim=ylim,xlim=xlim)
    mtext(bquote(F[ST]==.(Fst[f])), side=3, adj=0, line=1.2, cex=1, font=2);
    abline(0,1)
    legend(x = "bottomright",title = "N SNPs",legend = paste(c("1K","1K Neut","500","100","10"),eq),col = c("black","grey","blue","green","red"),pch = 16,bty = "n")
  }
  if(save){dev.off()}
}
FigureS13(save = T)






