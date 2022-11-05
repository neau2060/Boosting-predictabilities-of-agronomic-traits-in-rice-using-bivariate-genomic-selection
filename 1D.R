###### simplified 1d function ######
gsm.1d.simple <- function(mydata, mykin){
### mydata is a phenotypic vector
### mykin is a kinship matrix 
### Negative nloglikelihood Function
nloglik.REML.1d <- function(theta){

    lambda<-exp(theta)
    logdt<-sum(log(lambda*vv+1))
    h<-1/(lambda*vv+1)
    yy<-sum(yu*h*yu)
    yx<-matrix(0,s,1)
    xx<-matrix(0,s,s)
    for(i in 1:s){
      yx[i]<-sum(yu*h*xu[,i])
      for(j in 1:s){
        xx[i,j]<-sum(xu[,i]*h*xu[,j])
      }
    }
    loglike<- -0.5*logdt-0.5*(n-s)*log(yy-t(yx)%*%solve(xx)%*%yx)-0.5*log(det(xx))
    return(-loglike)
}


### Henderson Equation and Estimation
henderson.FUN <- function(t, y, kin){

  n <- length(y)
  x.design <- matrix(1, n, 1)

  phi <- t[1]
  sigma <- t[2]
  lambda <- phi/sigma

  TL <- t(x.design)%*%x.design
  BL <- x.design
  TR <- t(x.design)
  BR <- diag(n) + solve(kin)/lambda

  #v.est <- ginv(cbind(rbind(TL, BL), rbind(TR, BR)))
  
  com2x2 <- cbind(rbind(TL, BL), rbind(TR, BR))
  v.com2x2 <- eigen(com2x2)[[1]]
  if (min(abs(v.com2x2))>1e-8){   
  v.est <- solve(com2x2)
  }
  else{
  v.est <- solve(com2x2+diag(com2x2)*1e-8)
  }
  




  est <- v.est%*%matrix(c(t(x.design)%*%y, y), n+1, 1)
  
  beta_hat <- est[1, 1]
  eta_hat <- est[-1, 1]

  return(list(beta=beta_hat, eta=eta_hat, v=v.est))
}


x<-matrix(1,length(mydata),1)
s<-ncol(x)
n<-length(mydata)
uu<-eigen(mykin)[[2]]
vv<-eigen(mykin)[[1]]
xu<-t(uu)%*%x
yu<-t(uu)%*%mydata


  
### 1st search
Initial.value <- function(v0, v1, interval){
coef.comb <- as.matrix(seq(v0, v1, interval))


nloglik.REML.1d.batch <- function(v){
  theta <- v
  return(nloglik.REML.1d(theta))
}
nloglik.batch <- apply(coef.comb, 1, nloglik.REML.1d.batch)
nloglik.batch<-nloglik.batch[nloglik.batch!="NaN"]
coef.comb<-coef.comb[nloglik.batch!="NaN",1]

sel.id <- which(nloglik.batch==min(nloglik.batch))[1]
theta <- coef.comb[sel.id]
return(theta)
}

Hat.test <- function(lambda, yu, xu, vv, mydata, mykin){
n<-nrow(yu)
s<-ncol(xu)
h<-1/(lambda*vv+1)
yy<-sum(yu*h*yu)
yx<-matrix(0,s,1)
xx<-matrix(0,s,s)
for(i in 1:s){
  yx[i]<-sum(yu*h*xu[,i])
  for(j in 1:s){
    xx[i,j]<-sum(xu[,i]*h*xu[,j])
  }
}

sigma2 <- (yy-t(yx)%*%solve(xx)%*%yx)/(n-s)
par <- c(lambda*sigma2, sigma2)
var.para <- par

junk2 <- henderson.FUN(par, mydata, mykin)
beta.est <- junk2$beta
eta.est <- junk2$eta
v.est <- junk2$v
  
v.phi <- var.para[1]*mykin
v.sigma <- var.para[2]*diag(n)
v <- v.phi+v.sigma
  
hatMatrix <- v.phi%*%solve(v)
residual <- mydata - rep(beta.est, n) - eta.est  
rr <- cor(mydata,(rep(beta.est, n) + eta.est))

press<-sum((residual/(1-diag(hatMatrix)))^2)
ss <- sum((mydata-mean(mydata))**2)
predic.HAT <- 1-press/ss

gal <- list(predic.HAT=predic.HAT, var.para=var.para)
return(gal)
}

print("1D")
theta <- Initial.value(-1, 1, 0.1)
print(theta)
junk <- optim(par=theta,fn=nloglik.REML.1d,NULL,hessian = TRUE, method="L-BFGS-B",lower=-10,upper=10)
lambda <- exp(junk$par) 
gal <- Hat.test(lambda, yu, xu, vv, mydata, mykin)
var.para <- gal$var.para
predic.HAT <- gal$predic.HAT

res <- list(var.para=var.para, predic.HAT=predic.HAT)

### var.para is the estimated variance components
### predic.HAT is the estimated predictabilities of HAT
return(res)

}
