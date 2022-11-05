###### simplified 1d function ######
gsm.1d.simple <- function(mydata, mykin){
  
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
    loglike<- -0.5*logdt-0.5*(n-s)*log(yy-t(yx)%*%ginv(xx)%*%yx)-0.5*log(det(xx))
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
  BR <- diag(n) + ginv(kin)/lambda

  #v.est <- ginv(cbind(rbind(TL, BL), rbind(TR, BR)))
  
  com2x2 <- cbind(rbind(TL, BL), rbind(TR, BR))
  v.com2x2 <- eigen(com2x2)[[1]]
  if (min(abs(v.com2x2))>1e-8){   
  v.est <- ginv(com2x2)
  }
  else{
  v.est <- ginv(com2x2+diag(com2x2)*1e-8)
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
coef.space <- seq(v0, v1, interval)
coef.comb <- NULL
for(cg in coef.space){
  coef.comb <- rbind(coef.comb, cg)
}
nloglik.REML.1d.batch <- function(v){
  theta <- v
  return(nloglik.REML.1d(theta))
}
nloglik.batch <- apply(coef.comb, 1, nloglik.REML.1d.batch)
pos<-which(nloglik.batch!="NaN")
if ((length(pos))>0)
{nloglik.batch<-nloglik.batch[pos]
coef.comb<-coef.comb[pos,]
}
nloglik.optimal <- min(nloglik.batch)
sel.id <- which(nloglik.batch==nloglik.optimal)[1]
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

sigma2 <- (yy-t(yx)%*%ginv(xx)%*%yx)/(n-s)
par <- c(lambda*sigma2, sigma2)
var.para <- par

junk2 <- henderson.FUN(par, mydata, mykin)
beta.est <- junk2$beta
eta.est <- junk2$eta
v.est <- junk2$v
  
v.phi <- var.para[1]*mykin
v.sigma <- var.para[2]*diag(n)
v <- v.phi+v.sigma
  
hatMatrix <- v.phi%*%ginv(v)
residual <- mydata - rep(beta.est, n) - eta.est  
rr <- cor(mydata,(rep(beta.est, n) + eta.est))

press <- 0
for(i in 1:n){
  residual.tmp <- residual[i]
  hatMatrix.tmp <- hatMatrix[i, i]
  residual_pred.tmp <- residual.tmp/(1-hatMatrix.tmp)
  press <- press + residual_pred.tmp**2
}
ss <- sum((mydata-mean(mydata))**2)
predic.HAT <- 1-press/ss

gal <- list(predic.HAT=predic.HAT, var.para=var.para)
return(gal)
}


theta <- Initial.value(-1, 1, 0.1)
junk <- optim(par=theta,fn=nloglik.REML.1d,NULL,hessian = TRUE, method="L-BFGS-B",lower=-10,upper=10)
lambda <- exp(junk$par) 
gal <- Hat.test(lambda, yu, xu, vv, mydata, mykin)
var.para <- gal$var.para
predic.HAT <- gal$predic.HAT

res <- list(var.para=var.para, predic.HAT=predic.HAT)


return(res)

}

  

  



  
  
 ###### 2d function with HAT ######
gsmHAT.2d <- function(mydata, mykin){

### mydata is the phenotpic vector of the two traits with the first n elements being the phenotypes for the 1st trait and residual elements being the pehnotypes for the 2nd trait
### mykin is the kinship matrix

### Negative nloglikelihood Function with 6 variance compoments
nloglik.REML.2d <- function(t){
  
  phi.11 <- t[1]
  phi.12 <- t[2]
  phi.22 <- t[3]
  
  sigma.11 <- t[4]
  sigma.12 <- t[5]
  sigma.22 <- t[6]


  B11<-phi.11*vv+sigma.11
  B12<-phi.12*vv+sigma.12
  B22<-phi.22*vv+sigma.22

  I.inv<-(1/B22)
  I.block<-1/(B11-B12*I.inv*B12)
  A11<-diag(I.block)
  A12<-diag(-(I.block*(B12)*I.inv))
  A21<-diag(-(I.inv*(B12)*I.block))
  A22<-diag(I.inv+I.inv*(B12)*I.block*(B12)*I.inv)
  
  v_inv<-rbind(cbind(A11,A12),cbind(A21,A22))
  v_det<-sum(log(abs(B11)))+sum(log(abs(B22-((B12*B12)*(1/B11)))))

  beta <- ginv(t(XU)%*%(v_inv)%*%XU)%*%(t(XU)%*%(v_inv)%*%YU)
  nloglik = 0.5*(v_det + unlist(determinant(t(XU)%*%(v_inv)%*%XU))[[1]] + t(YU-XU%*%beta)%*%(v_inv)%*%(YU-XU%*%beta))
  
  return(nloglik)
}

### Henderson Equation and Estimation
henderson.FUN <- function(t, y, kin.inv, x.design){
  n <- length(y)/2

  phi.11 <- t[1]
  phi.12 <- t[2]
  phi.22 <- t[3]
  
  sigma.11 <- t[4]
  sigma.12 <- t[5]
  sigma.22 <- t[6]

  G.matrix <- matrix(c(phi.11, phi.12, phi.12, phi.22), 2, 2)
  R.matrix <- matrix(c(sigma.11, sigma.12, sigma.12, sigma.22), 2, 2)


  
  TL <- t(x.design)%*%kronecker(ginv(R.matrix), diag(n))%*%x.design
  BL <- kronecker(ginv(R.matrix), diag(n))%*%x.design
  TR <- t(x.design)%*%kronecker(ginv(R.matrix), diag(n))
  BR <- kronecker(ginv(R.matrix), diag(n)) + kronecker(ginv(G.matrix), kin.inv)

  com2x2 <- cbind(rbind(TL, BL), rbind(TR, BR))
  v.com2x2 <- eigen(com2x2)[[1]]
  if (min(abs(v.com2x2))>1e-8){   
  v.est <- ginv(com2x2)
  }
  else{
  v.est <- ginv(com2x2+diag(com2x2)*1e-8)
  }
  v_inv.est <- cbind(rbind(TL, BL), rbind(TR, BR))
  est <- v.est%*%rbind(TR%*%matrix(y, 2*n, 1), kronecker(ginv(R.matrix), diag(n))%*%matrix(y, 2*n, 1))
  
  beta_hat <- est[1:2, 1]
  eta_hat <- est[-(1:2), 1]

  return(list(beta=beta_hat, eta=eta_hat, v=v.est))
}


Hat.test <- function(var.para, mydata, mykin, kin.inv, UU, vv, x.design){
n <- length(mydata)/2
y.1 <- mydata[1:n]
y.2 <- mydata[-(1:n)]
junk2 <- henderson.FUN(var.para, mydata, kin.inv, x.design)
beta.est <- junk2$beta
eta.est <- junk2$eta
v.est <- junk2$v
  
v.phi <- kronecker(matrix(c(var.para[1], var.para[2], var.para[2], var.para[3]), 2, 2), mykin)

B11<-var.para[1]*vv+var.para[4]
B12<-var.para[2]*vv+var.para[5]
B22<-var.para[3]*vv+var.para[6]

I.inv<-(1/B22)
I.block<-1/(B11-B12*I.inv*B12)
A11<-diag(I.block)
A12<-diag(-(I.block*(B12)*I.inv))
A21<-diag(-(I.inv*(B12)*I.block))
A22<-diag(I.inv+I.inv*(B12)*I.block*(B12)*I.inv)
  
v.inv<-UU%*%rbind(cbind(A11,A12),cbind(A21,A22))%*%t(UU)


hatMatrix <- v.phi%*%v.inv
residual <- mydata - c(rep(beta.est[1], n), rep(beta.est[2], n)) - eta.est
  
press_1 <- 0
press_2 <- 0
for(i in 1:n){
   residual.tmp <- residual[c(i, n+i)]
   hatMatrix.tmp <- hatMatrix[c(i, n+i), c(i, n+i)]
   residual_pred.tmp <- ginv(diag(2)-hatMatrix.tmp)%*%residual.tmp
   press_1 <- press_1 + residual_pred.tmp[1]**2
   press_2 <- press_2 + residual_pred.tmp[2]**2
}
  
ss_1 <- sum((y.1-mean(y.1))**2)
ss_2 <- sum((y.2-mean(y.2))**2)  
  
predic_hat_1 <- 1-press_1/ss_1
predic_hat_2 <- 1-press_2/ss_2
predic.HAT=c(predic_hat_1, predic_hat_2)
return(predic.HAT)

}





n <- length(mydata)/2
y.1 <- mydata[1:n]
y.2 <- mydata[-(1:n)]

kin.eigen<-eigen(mykin)
uu<-kin.eigen[[2]]
vv<-kin.eigen[[1]] 
x.design <- matrix(0, 2*n, 2)
for(i in 1:2){
  x.design[(((i-1)*n+1):(i*n)), i] <- 1
}
UU<-rbind(cbind(uu,uu*0),cbind(uu*0,uu))
XU <- t(UU)%*%x.design
YU <- t(UU)%*%mydata
kin.inv <- uu%*%diag(1/vv)%*%t(uu)



result1 <- gsm.1d.simple(y.1, kin)
result2 <- gsm.1d.simple(y.2, kin)

Initial.value <- function(v0, v1, interval){

coef.space <- seq(v0, v1, interval)
coef.comb <- NULL
cgce.comb <- NULL
for(cg in coef.space){
	for(ce in coef.space){
		coef.comb <- rbind(coef.comb, c(result1$var.para[1], cg*(result1$var.para[1]*result2$var.para[1])**0.5, result2$var.para[1], result1$var.para[2], ce*(result1$var.para[2]*result2$var.para[2])**0.5, result2$var.para[2]))
		cgce.comb <- rbind(cgce.comb, c(cg, ce))
	}
}
nloglik.REML.2d.batch <- function(v){
	t.init0 <- v
	return(nloglik.REML.2d(t.init0))
}

nloglik.batch <- apply(coef.comb, 1, nloglik.REML.2d.batch)
pos<-which(nloglik.batch!="NaN")
if ((length(pos))>0)
{nloglik.batch<-nloglik.batch[pos]
cgce.comb<-cgce.comb[pos,]
coef.comb<-coef.comb[pos,]
}

nloglik.optimal <- min(nloglik.batch)
sel.id <- which(nloglik.batch==nloglik.optimal)[1]
t.init <- coef.comb[sel.id,]
return(t.init)

}

t.init <- Initial.value(-1, 1, 0.1)
junk <- optim(par=t.init, fn=nloglik.REML.2d, method="L-BFGS-B", lower=c(0.2*result1$var.para[1], -10*(result1$var.para[1]*result2$var.para[1])**0.5, 0.2*result2$var.para[1], 0.2*result1$var.para[2], -10*(result1$var.para[2]*result2$var.para[2])**0.5, 0.2*result2$var.para[2]), upper=c(5*result1$var.para[1], 10*(result1$var.para[1]*result2$var.para[1])**0.5, 5*result2$var.para[1], 5*result1$var.para[2], 10*(result1$var.para[2]*result2$var.para[2])**0.5, 5*result2$var.para[2]))
var.para <- junk$par
Hat1 <- Hat.test(var.para, mydata, mykin, kin.inv, UU, vv, x.design)
Hat <- Hat1
if (Hat1[1]<0 | Hat1[2]<0){
t.init <- Initial.value(-0.1, 0.1, 0.01)
junk <- optim(par=t.init, fn=nloglik.REML.2d, method="L-BFGS-B", lower=c(0.2*result1$var.para[1], -10*(result1$var.para[1]*result2$var.para[1])**0.5, 0.2*result2$var.para[1], 0.2*result1$var.para[2], -10*(result1$var.para[2]*result2$var.para[2])**0.5, 0.2*result2$var.para[2]), upper=c(5*result1$var.para[1], 10*(result1$var.para[1]*result2$var.para[1])**0.5, 5*result2$var.para[1], 5*result1$var.para[2], 10*(result1$var.para[2]*result2$var.para[2])**0.5, 5*result2$var.para[2]))
var.para <- junk$par
Hat2 <- Hat.test(var.para, mydata, mykin, kin.inv, UU, vv, x.design)
Hat <- Hat2
if (Hat2[1]<0 | Hat2[2]<0){
t.init <- Initial.value(-0.01, 0.01, 0.001)
junk <- optim(par=t.init, fn=nloglik.REML.2d, method="L-BFGS-B", lower=c(0.2*result1$var.para[1], -10*(result1$var.para[1]*result2$var.para[1])**0.5, 0.2*result2$var.para[1], 0.2*result1$var.para[2], -10*(result1$var.para[2]*result2$var.para[2])**0.5, 0.2*result2$var.para[2]), upper=c(5*result1$var.para[1], 10*(result1$var.para[1]*result2$var.para[1])**0.5, 5*result2$var.para[1], 5*result1$var.para[2], 10*(result1$var.para[2]*result2$var.para[2])**0.5, 5*result2$var.para[2]))
var.para <- junk$par
### var.para are the variance components
Hat3 <- Hat.test(var.para, mydata, mykin, kin.inv, UU, vv, x.design)
Hat <- Hat3
}
}

All.1D<-c(result1$predic.HAT,result2$predic.HAT)
All.2D<-c(Hat[1],Hat[2])

if (Hat[1]<0 & Hat[2]<0)
{HaT<-All.1D
}
if (Hat[1]*Hat[2]<0)
{Pos.T<-which(ALL.2D<0)
if (All.1D[Pos.T]>0)
{HaT<-All.1D
}
else
{HaT<-All.2D
}
}

if (Hat[1]>0 & Hat[2]>0)
{HaT<-All.2D
}

### Hat is the estimated predictabilities of HAT for the two traits

return(HaT)
}