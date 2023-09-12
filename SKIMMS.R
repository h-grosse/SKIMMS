require(MASS)
require(evmix)
require(Iso)
require(latex2exp)
require(parallel)
require(pracma)
require(readr)
require(stringr)
require(glmnet)
require(caret)
require(magrittr)
require(stringr)
require(MASS)
require(evmix)
require(Iso)
require(latex2exp)
require(parallel)
require(pracma)
require(plotly)

         ###############
         #### Notes ####
         ###############

# 1) optimtauS is optimising over ||beta||_1<=1 and rescaling
#    -> try to make it optimise over ||beta||_2=1
#
# 2) The noise e_i are generated from a standard normal
#    -> try to incorporate the 10% outliers following the Cauchy distribution
#
# To add: A) 5-fold CV 6) parallel computations?
#

createdata <- function(n, p, rho){
  S <- vector("list",length(rho))
  
  for (i in 1:length(rho)) {
    S[[i]] <- matrix(rho[i], nrow = p, ncol = p)
    diag(S[[i]]) <- rep(1,p)
  }  
  
  X <- vector("list",length(n))
  
  for (k in 1:length(n)) {
    X[[k]] <- vector("list",length(rho))
    
    for (i in 1:length(rho)) {
      X[[k]][[i]] <- matrix(NA, nrow = p, ncol = n[k])
      for (j in 1:n[k]) {
        X[[k]][[i]][,j] <- mvrnorm(n = 1, mu = rep(0,p), Sigma = S[[i]])
      }
    }
  }
  
  nbex <- 10
  #--------------------------------------------------
  beta <- vector("list",nbex)
  trueindex <- vector("list",nbex)
  truecoef <- vector("list",nbex)
  beta[[1]] <- rep(0, p)
  trueindex[[1]] <- c(1,2,3,4)
  truecoef[[1]] <- c(3,5,3,5)
  truecoef[[1]] <- truecoef[[1]]/norm(truecoef[[1]], type = "2")
  beta[[1]][trueindex[[1]]] <- truecoef[[1]]
  #-------------------------------------------------
  beta[[2]] <- rep(0, p)
  trueindex[[2]] <- c(1,2,3,4,5,6)
  truecoef[[2]] <- c(3,3,-4,-4,5,5)
  truecoef[[2]] <- truecoef[[2]]/norm(truecoef[[2]], type = "2")
  beta[[2]][trueindex[[2]]] <- truecoef[[2]]
  #-------------------------------------------------
  beta[[3]] <- rep(0, p)
  trueindex[[3]] <- c(1,2,3,4,5,6,7,8,9)
  truecoef[[3]] <- c(3,3,3,-4,-4,-4,5,5,5)
  truecoef[[3]] <- truecoef[[3]]/norm(truecoef[[3]], type = "2")
  beta[[3]][trueindex[[3]]] <- truecoef[[3]]
  #-------------------------------------------------
  beta[[4]] <- rep(0, p)
  trueindex[[4]] <- c(1,2,3,4,5,6)
  truecoef[[4]] <- c(2,-4,-4,2,3,3)
  truecoef[[4]] <- truecoef[[4]]/norm(truecoef[[4]], type = "2")
  beta[[4]][trueindex[[4]]] <- truecoef[[4]]
  #-------------------------------------------------
  beta[[5]] <- rep(0, p)
  trueindex[[5]] <- c(1,2,3,4,5,6)
  truecoef[[5]] <- c(2,-4,-4,2,3,3)
  truecoef[[5]] <- truecoef[[5]]/norm(truecoef[[5]], type = "2")
  beta[[5]][trueindex[[5]]] <- truecoef[[5]]
  #-------------------------------------------------
  beta[[6]] <- rep(0, p)
  trueindex[[6]] <- c(1,2,3,4)
  truecoef[[6]] <- c(3,-4,-4,5)
  truecoef[[6]] <- truecoef[[6]]/norm(truecoef[[6]], type = "2")
  beta[[6]][trueindex[[6]]] <- truecoef[[6]]
  #-------------------------------------------------
  beta[[7]] <- rep(0, p)
  trueindex[[7]] <- c(1,2,3,4)
  truecoef[[7]] <- c(3,-4,-4,5)
  truecoef[[7]] <- truecoef[[7]]/norm(truecoef[[7]], type = "2")
  beta[[7]][trueindex[[7]]] <- truecoef[[7]]
  #-------------------------------------------------
  beta[[8]] <- rep(0, p)
  trueindex[[8]] <- c(1,2,3,4)
  truecoef[[8]] <- c(3,-4,-4,5)
  truecoef[[8]] <- truecoef[[8]]/norm(truecoef[[8]], type = "2")
  beta[[8]][trueindex[[8]]] <- truecoef[[8]]
  #-------------------------------------------------
  beta[[9]] <- rep(0, p)
  trueindex[[9]] <- c(1,2,3,4)
  truecoef[[9]] <- c(3,-4,-4,5)
  truecoef[[9]] <- truecoef[[9]]/norm(truecoef[[9]], type = "2")
  beta[[9]][trueindex[[9]]] <- truecoef[[9]]
  #-------------------------------------------------
  beta[[10]] <- rep(0, p)
  trueindex[[10]] <- c(1,2,3,4)
  truecoef[[10]] <- c(3,-4,-4,5)
  truecoef[[10]] <- truecoef[[10]]/norm(truecoef[[10]], type = "2")
  beta[[10]][trueindex[[10]]] <- truecoef[[10]]
  #-------------------------------------------------
  
  err <- vector("list",nbex)
  err[[1]] <- rnorm(n)
  err[[2]] <- rnorm(n)
  err[[3]] <- rnorm(n)
  err[[4]] <- rnorm(n)
  err[[5]] <- rnorm(n)
  err[[6]] <- rnorm(n)
  err[[7]] <- rnorm(n)
  err[[8]] <- rnorm(n)
  err[[9]] <- rnorm(n)
  err[[10]] <- rnorm(n)
  
  f4 <- function(x){0.45*x^3}
  f5 <- function(x){1.2*nthroot(x,3)}
  f6 <- function(x){3*sin(10*x)/10 + 3*x}
  f7 <- function(x){2*(x+1/2)^2-1}
  f8 <- Vectorize(function(x){
    x <- 1.5*x
    if(x<=-3/2){return(4*x+10/2)}
    if(x<=-1/2){return(x+1/2)}
    if(x<=1/2){return(0)}
    if(x<=3/2){return(x-1/2)}
    return(4*x-10/2)},vectorize.args = "x")
  f9 <- function(x){x}
  f10 <- function(x){x}
  
  Y <- list(atan(2*t(X[[1]][[1]])%*%beta[[1]]) + 0.2*err[[1]],
            0.5*exp(t(X[[1]][[1]])%*%beta[[2]]) + 0.5*err[[2]],
            exp(atan(2*t(X[[1]][[1]])%*%beta[[3]])) + 0.2*err[[3]],
            f4(t(X[[1]][[1]])%*%beta[[4]]) + 0.2*err[[4]],
            f5(t(X[[1]][[1]])%*%beta[[5]]) + 0.5*err[[5]],
            f6(t(X[[1]][[1]])%*%beta[[6]]) + 0.2*err[[6]],
            f7(t(X[[1]][[1]])%*%beta[[7]]) + 0.2*err[[7]],
            f8(t(X[[1]][[1]])%*%beta[[8]]) + 0.2*err[[8]],
            f9(t(X[[1]][[1]])%*%beta[[9]]) + 0.2*err[[9]],
            f10(t(X[[1]][[1]])%*%beta[[10]]) + 0.2*err[[10]])
  
  return(list(X = X,Y = Y, beta = beta, index = trueindex, coef = truecoef))
}
tau <- function(Y, X, beta, comment = FALSE) {
  n <- length(Y)
  if(n != dim(X)[2] | dim(X)[1] != length(beta)) {return("Dimension Error")}
  SS <- 0
  
  Xb <- t(X)%*%beta
  for (ione in 1:n) {
    for (itwo in 1:n) {
      SS <- SS + sign(Y[itwo]-Y[ione])*sign(Xb[itwo]-Xb[ione])
    }
  }
  
  return( SS/(n*(n-1)) )
}
FalseCalc <- function(true,test){
  
  strue <- length(true)
  stest <- length(test)
  FN <- 0
  
  for (k in 1:strue) {
    FN <- FN + c(rep(0,length(test))[true[k] == test],1)[1]
  }
  TP <- strue - FN
  FP <- stest - TP
  
  return(list(FP = as.numeric(FP), FN = as.numeric(FN), TP = as.numeric(TP)))
}
tauS <- function(Y, X, beta, h, comment = FALSE) {
  n <- length(Y)
  if(n != dim(X)[2] | dim(X)[1] != length(beta)) {return("Dimension Error")}
  SS <- 0
  
  Xb <- t(X)%*%beta
  for (ione in 1:n) {
    for (itwo in 1:n) {
      SS <- SS + sign(Y[itwo]-Y[ione])*tanh((Xb[itwo]-Xb[ione])/h)
    }
  }

  return( SS/(n*(n-1)) )
}
optimtauS <- function(X, Y, index, p, h){
  size <- length(index)
  f <- function(x){
    bet <- rep(0, p)
    bet[index] <- x
    -tauS(Y,X,bet, h)
  }
  minv <- optim(rep(0,size),f,method = "L-BFGS-B",lower = rep(-1,size), upper = rep(1,size))$par
  minv <- minv/norm(minv, type = "2")
  bb <- rep(0,p)
  bb[index] <- minv
  return(bb)
}
stepalgo <- function(X, Y, beta, index, lambda, eps, h, comment = FALSE) {

  p <- dim(X)[1]
  n <- dim(X)[2]
  if(n != length(Y)){return("Dimension Error")}

  temp <- numeric(p-length(index))
  betaopt <- numeric(p-length(index))
  ind <- (1:p)[-index]
  for (var in 1:(p-length(index))) {
    betaopt[var] <- optimize(function(b) tauS(Y, X, beta + b*diag(p)[,ind[var]], h)-lambda*abs(b),
                             lower = -5, upper = 5, maximum = TRUE)$maximum
    temp[var] <- tauS(Y, X, beta+betaopt[var]*diag(p)[,ind[var]], h)
    if(comment){print(paste0("lambda=",lambda,"/var=",var,"/bopt=",round(betaopt[var],digits = 3),
                             "/tau(beta - bopt*e_j)=",round(temp[var],digits = 3)))}
  }
  wmax <- which.max(temp)
  maxj <- ind[wmax]
  if(comment){print(paste0("j_max = ",maxj, " / betamax = ",betaopt[wmax]))}
  be <- beta+betaopt[wmax]*diag(p)[,maxj]
  if(tauS(Y, X, be, h)-tauS(Y, X, beta, h)<eps) {
    if(comment){
      print(paste0("tauS(beta_new) = ",tauS(Y, X, be, h)," / tauS(beta) = ",tauS(Y, X, beta, h)))
      print(paste0("Succesful iteration: tauS(beta_new)-tauS(beta) = ",
                 tauS(Y, X, be, h)-tauS(Y, X, beta, h)," < eps"))
    }
    return(list(tauS(Y, X, be, h)-tauS(Y, X, beta, h),"STOP"))
  }
  be <- be/norm(be, type = "2")
  if(comment){
    print(paste0("tauS(beta_new) = ",tauS(Y, X, be, h)," / tauS(beta) = ",tauS(Y, X, beta, h)))
    print(paste0("Succesful iteration: tauS(beta_new)-tauS(beta) = ",
               tauS(Y, X, be, h)-tauS(Y, X, beta, h)," >= eps"))
  }
  return(list(be,maxj))
}
link_estimate <- function(X,Y,beta,Xtest,Ytest){
  
  p <- length(beta)
  n <- length(Y)
  index <- (1:p)[!(beta == 0)]
  
  Z <- t(X)%*%beta
  Zord <- Z[order(Z)]
  Yord <- Y[order(Z)]
  
  Ytop <- pava(Yord, stepfun = FALSE)
  
  g <- function(t,b){
    top <- 0
    bottom <- 0
    for (j in 1:n) {
      top <- top + Ytop[j]*kdgaussian((t-Zord[j])/b)
      bottom <- bottom + kdgaussian((t-Zord[j])/b)
    }
    return(top/bottom)
  }
  
  bot <- function(t,b){
    bottom <- 0
    for (j in 1:n) {
      bottom <- bottom + kdgaussian((t-Zord[j])/b)
    }
    return(bottom)
  }
  trS <- function(b){
    top <- 0
    for (j in 1:n) {
      top <- top + kdgaussian(0)/bot(Zord[j],b)
    }
    return(top)
  }
  link.gcv <- function(b){
    return( sum((Y-g(Z,b))^2)/(n*(1-trS(b)/n)^2) )
  }

  bopt <- optimize(link.gcv, interval = c(0.001,5))$minimum
  
  MSE <- mean((g(t(Xtest)%*%beta,bopt)-Ytest)^2)
  
  return(list(link = function(t) g(t,bopt), PE = MSE, Z = Z,Y = Y))
  
}
plotlink <- function(X,Y,beta,g,example,lambda,eps,sel,sim){
  
  f4 <- function(x){0.45*x^3}
  f5 <- function(x){1.2*nthroot(x,3)}
  f6 <- function(x){3*sin(10*x)/10 + 3*x}
  f7 <- function(x){2*(x+1/2)^2-1}
  f8 <- Vectorize(function(x){
    x <- 1.5*x
    if(x<=-3/2){return(4*x+10/2)}
    if(x<=-1/2){return(x+1/2)}
    if(x<=1/2){return(0)}
    if(x<=3/2){return(x-1/2)}
    return(4*x-10/2)},vectorize.args = "x")
  f9 <- function(x){x}
  f10 <- function(x){x}
  
  var <- (1:length(beta))[beta != 0]
  
  if(length(beta) == 13){
    plot(g, xlim = c(min(t(X)%*%beta)-0.1,max(t(X)%*%beta)+0.1), 
                              ylim = c(min(Y)-2,max(Y)+2), col = "blue", lwd = 2, 
                              xlab = TeX(r'($X^T \widehat{\beta}$)'), ylab = "")
    text(min(t(X)%*%beta)-0.1, max(Y)+1.8, paste0("Selected Var = ",str_c(var, collapse = ", ")),cex=1, pos=4) 
    }
  if(length(beta) != 13){
    plot(g, xlim = c(min(t(X)%*%beta)-0.1,max(t(X)%*%beta)+0.1), 
       ylim = c(min(Y)-2,max(Y)+2), col = "blue", lwd = 2, 
       xlab = TeX(r'($X^T \beta$)'), ylab = "")
    text(min(t(X)%*%beta)-0.1, max(Y)+1.9, paste0("True Var = ",str_c(var, collapse = ", ")),cex=1, pos=4) 
    text(min(t(X)%*%beta)-0.1, max(Y)+1.2, paste0("Selected Var = ",str_c(sel, collapse = ", ")),cex=1, pos=4) 
    }
  title(main = paste0("Example ",example," - lambda = ",lambda," / eps = ",eps,"  (Sim = ",sim,")" ),
        ylab = TeX(r'($g(X^T \beta)$)'),line=2.4, cex.lab=1)
  xseqbis <- seq(min(t(X)%*%beta)-1,max(t(X)%*%beta)+1,0.1)
  if(example == 1){lines(xseqbis,atan(2*xseqbis),col = "red", lty = 2)}
  if(example == 2){lines(xseqbis,0.5*exp(xseqbis),col = "red", lty = 2)}
  if(example == 3){lines(xseqbis,exp(atan(2*xseqbis)),col = "red", lty = 2)}
  if(example == 4){lines(xseqbis,0.45*xseqbis^3,col = "red", lty = 2)}
  if(example == 5){lines(xseqbis,f5(xseqbis),col = "red", lty = 2)}
  if(example == 6){lines(xseqbis,f6(xseqbis),col = "red", lty = 2)}
  if(example == 7){lines(xseqbis,f7(xseqbis),col = "red", lty = 2)}
  if(example == 8){lines(xseqbis,f8(xseqbis),col = "red", lty = 2)}
  if(example == 9){lines(xseqbis,f9(xseqbis),col = "red", lty = 2)}
  if(example == 10){lines(xseqbis,f10(xseqbis),col = "red", lty = 2)}
  points(t(X)%*%beta, Y)
  if(length(beta) != 13){legend(x= "bottomright",c(TeX(r'($\widehat{g}$)'),TeX(r'($g_\true$)'),"samples"),
         lty = c(1,2,NA), pch = c(NA,NA,1), lwd = c(2,1,NA), col = c("blue","red","black"))}
  
}
printResult <- function(ResSim, scale = "log"){
  if(dim(ResSim$Table)[1] == 1){
     print.table(ResultSim$Table) 
  }
  if(dim(ResSim$Table)[1] != 1){
    cat("  -----------------------------------------------------------\n")
    cat(c("  n =",ResSim$param["n"],", p =",ResSim$param["p"],", eps =",ResSim$param["eps"],
          ", h =",ResSim$param["h"],", sim =",ResSim$param["sim"],"\n"))
    exsize <- ResSim$param["exsize"]
    lsize <- ResSim$param["lsize"]
    cat("  -----------------------------------------------------------\n")
    par(mfrow = c(1,2))
    par(pty = "s")
    for (i in 1:exsize) {
      M <- cbind(ResSim$Table[ResSim$Table[,1] == unique(ResSim$Table[,1])[i],],rep("|",lsize))
      colnames(M)[dim(M)[2]] <- "|"
      rownames(M) <- rep("|",dim(M)[1])
      print(noquote(M))
      cat("  -----------------------------------------------------------\n")

      maxPE <- max(as.numeric(M[,"PE"]))
      minPE <- min(as.numeric(M[,"PE"]))
      maxFP <- max(as.numeric(M[,"FP"]))
      minFP <- min(as.numeric(M[,"FP"]))
      maxFN <- max(as.numeric(M[,"FN"]))
      minFN <- min(as.numeric(M[,"FN"]))
      maxTP <- max(as.numeric(M[,"TP"]))
      minTP <- min(as.numeric(M[,"TP"]))
      maxALL <- max(maxFP,maxFN,maxTP)
      minALL <- min(minFP,minFN,minTP)
      maxl <- max(as.numeric(M[,"lambda"]))
      minl <- min(as.numeric(M[,"lambda"]))
     
      ord <- order(M[,"lambda"])
      
      if(scale == "log"){plot(M[ord,"lambda"],M[ord,"PE"], main = paste0("Example ",i), pch = 16, type = "b", xlim = c(maxl,minl),ylim = c(minPE,maxPE), xlab = TeX(r'($\lambda)'), ylab = "PE", log = 'x', xaxt = "n")}
      if(scale == "lin"){plot(M[ord,"lambda"],M[ord,"PE"], main = paste0("Example ",i), pch = 16, type = "b", xlim = c(maxl,minl),ylim = c(minPE,maxPE), xlab = TeX(r'($\lambda)'), ylab = "PE", xaxt = "n")}
      axis(1, at = rev(M[ord,"lambda"]))
      if(scale == "log"){plot(M[ord,"lambda"],M[ord,"FP"], main = paste0("Example ",i), pch = 16, type = "b", xlim = c(maxl,minl), ylim = c(minALL,maxALL),xlab = TeX(r'($\lambda)'), ylab = "Average", col = "red", log = 'x', xaxt = "n")}
      if(scale == "lin"){plot(M[ord,"lambda"],M[ord,"FP"], main = paste0("Example ",i), pch = 16, type = "b", xlim = c(maxl,minl), ylim = c(minALL,maxALL),xlab = TeX(r'($\lambda)'), ylab = "Average", col = "red", xaxt = "n")}
      lines(M[ord,"lambda"],M[ord,"TP"], pch = 16, type = "b", col = "blue")
      lines(M[ord,"lambda"],M[ord,"FN"], pch = 16, type = "b", col = "purple")
      axis(1, at = rev(M[ord,"lambda"]))
      legend(x = "topright", legend = c("TP","FP","FN"), pch = rep(16,3), col = c("blue","red","purple"), cex = 0.7)
    }
    par(mfrow = c(1,1))
    par(pty = "m")
    
  }
  
}
tau_optim_plot <- function(n = 100, p = 150, h = 0.2, trueindex = c(1,2), truecoef = c(3,5), dimtooptim = c(1,2)){
  
  set.seed(15)
  data <- createdata(n,p,0.25)
  X <- data$X[[1]][[1]]
  beta <- rep(0, p)
  truecoef <- truecoef/norm(truecoef, type = "2")
  beta[trueindex] <- truecoef
  err <- rnorm(n)
  f <- function(x){3*sin(10*x)/10 + 3*x}
  Y <- f(t(X)%*%beta)+0.2*err
  
  fminus <- function(x){
    bet <- rep(0, p)
    bet[dimtooptim] <- x
    -tauS(Y,X,bet, h)
  }
  f <- function(x){
    bet <- rep(0, p)
    bet[dimtooptim] <- x
    tauS(Y,X,bet, h)
  }
  
  optimvec <- optimtauS(X,Y,dimtooptim,p,h)
  cat("beta_optim = (",optimvec[dimtooptim][1],",",optimvec[dimtooptim][2],") => f(beta_optim) = ",f(optimvec[dimtooptim]),"\n")
  cat("beta_true = (",truecoef[1],",",truecoef[2],") => f(beta_true) = ",f(beta[dimtooptim]),"\n")
  
  sx <- seq(-40,40,1)
  sy <- seq(-40,40,1)
  grid <- expand.grid(sx,sy)
  sz <- apply(grid,1,f)
  mz <- t(matrix(sz, nrow = length(sx)))
  
  #plotly -> https://plotly.com/r/3d-surface-plots/#surface-plot-with-contours
  require(plotly)
  fig <- plot_ly(x = sx, y = sy, z = mz) %>% 
    add_surface(contours = list(z = list(show=TRUE, start = f(optimvec[dimtooptim])-0.2,
                                         end = f(optimvec[dimtooptim])+0.2, size = 0.02, 
                                         color = "#b300ff",usecolormap=FALSE,
                                         highlightcolor="#ff0000",project=list(z=TRUE)))) %>% 
    layout(scene = list(xaxis = list(title ='beta_i'), yaxis = list(title = 'beta_j'),
                        zaxis = list(title = 'tau*(0,...,beta_i,...,beta_j,...,0)')))
  
  return(fig)
}

#' SKIMMS
#'
#' Smoothed Kendall’s association Iterative Maximizer Model Selector
#'
#' @param X The design matrix X 
#' @param Y The response vector Y
#' @param lambda The penalization parameter lambda > 0
#' @param eps The precision parameter epsilon > 0
#' @param h The smoothing parameter h > 0
#' @param comment Add comments?
#'
#' @return The beta obtained at the end of step 2, the variables chosen and the re-estimated beta
#' @export
#'
#' @examples
#' # Kendall Maximizer Single-Index estimation
#' SKIMMS(X,Y)
SKIMMS <- function(X,Y,lambda,eps,h, comment = FALSE){
  start_time <- Sys.time()
  p <- dim(X)[1]
  n <- dim(X)[2]
  if(n != length(Y)){return("Wrong dimensions")}
  T <- numeric(p)
  for (j in 1:p) {
    T[j] <- tau(Y, X, diag(p)[,j])
  }
  if(comment){
    print("|---------------------------------------------------------------------------|")
    print(paste0("Parameters: lambda = ", lambda," / eps = ",eps," / h = ",h))
  }
  
  index <- which.max(abs(T))
  betah <- sign(T[index])*diag(p)[,index]
  if(comment){
    print(paste0("Step 1: max(Tj) = ", index))
    print(paste0("--------------------------------"))
    print(paste0("Step 2: Iteration 1 starting... "))
  }
  
  templist <- stepalgo(X, Y, betah, index, lambda, eps, h, comment = comment)
  if(templist[[2]] == "STOP"){
    if(comment){
      print(paste0("--------------------------------"))
      print(paste0("Step 2: Ended"))
    }
    end_time <- Sys.time()
    return(list(beta = betah, var = index, beta_reEstimate = optimtauS(X, Y, index, p, h),
                run_time = (end_time - start_time)))
  }
  
  index <- rbind(index,templist[[2]])
  betah <- templist[[1]]
  
  if(comment){print(paste0("Index added -> ",templist[[2]]))}
  
  for (i in 1:round(1/eps)) {
    if(comment){
      print(paste0("--------------------------------"))
      print(paste0("Step 2: Iteration ",i+1," starting... "))
    }
    templist <- stepalgo(X, Y, betah, index, lambda, eps, h, comment = comment)
    if(templist[[2]] == "STOP"){
      if(comment){
        print(paste0("--------------------------------"))
        print(paste0("Step 2: Ended"))
      }
      end_time <- Sys.time()
      return(list(beta = betah, var = index, beta_reEstimate = optimtauS(X, Y, index, p, h),
                  run_time = (end_time - start_time)))
    }
    index <- rbind(index,templist[[2]])
    betah <- templist[[1]]
    
    if(comment){print(paste0("Index added -> ",templist[[2]]))}
    
  }
  
}
simulation <- function(sim = 1, n = 60,p = 80,rho = 0.25, lambda = 0.1, eps = 0.001, h = 0.2, example = 1, comment = FALSE, seed = 10){
  
  set.seed(seed)
  
  nrows <- length(lambda)*length(example); Times <- 1
  TP <- matrix(NA, nrow = nrows, ncol = sim); FP <- matrix(NA, nrow = nrows, ncol = sim); FN <- matrix(NA, nrow = nrows, ncol = sim); PE <- matrix(NA, nrow = nrows, ncol = sim)
  cuberoot = function(x){nthroot(x,3)}
  start_T <- proc.time()[[3]]
  print(noquote(paste0("Starting Simulation: sim = ",sim," / n = ",n," / p = ",p," / eps = ",eps)))
  
  for (k in 1:sim) {
    data <- createdata(n,p,rho)
    test <- createdata(400,p,rho)
    Xtest <- test$X[[1]][[1]]
    Ytest <- test$Y
    X <- data$X
    Y <- data$Y
    beta <- data$beta
    trueindex <- data$index
    truecoef <- data$coef
    
    for (j in 1:length(example)) {
      for (i in 1:length(lambda)) {
        start_T2 <- proc.time()[[3]]
        
        rem_iter <- length(example)*length(lambda)*sim - (length(lambda)*length(example)*(k-1)+length(lambda)*(j-1)+i)+1
        print(noquote(paste0("Sim Iteration: ",k,"/",sim," - Example ",example[j]," - lambda ", lambda[i]," - Tot Iteration ",length(lambda)*length(example)*(k-1)+length(lambda)*(j-1)+i,
                     "/",length(example)*length(lambda)*sim," - Estimated Time: ",floor(mean(Times)*rem_iter)," mins ",round((mean(Times)*rem_iter-floor(mean(Times)*rem_iter))*60)," secs")))
        
        idx <- i+length(lambda)*(j-1)
        
        Result <- SKIMMS(X[[1]][[1]], Y[[example[j]]], lambda[i], eps, h, comment = comment)
        run_time <- Result$run_time
      
        #print(paste0("Running time: ",round(run_time)," min ",round(60*(run_time-round(run_time)))," sec"))
        if(comment){print(noquote(paste0("Running time: ",run_time)))}
        FP[idx,k] <- FalseCalc(trueindex[[example[j]]],Result$var)$FP
        FN[idx,k] <- FalseCalc(trueindex[[example[j]]],Result$var)$FN
        TP[idx,k] <- FalseCalc(trueindex[[example[j]]],Result$var)$TP
        
        estim <- link_estimate(X[[1]][[1]],Y[[example[j]]],Result$beta_reEstimate,Xtest,Ytest[[example[j]]])
        PE[idx,k] <- round(estim$PE, digits = 2)
        
        if(k < 4){
          plotlink(X[[1]][[1]],Y[[example[j]]],beta[[example[j]]],estim$link,example[j],lambda[i],eps,Result$var,k)
        }
       
        end_T2 <- proc.time()[[3]]; IterTime <- (end_T2 - start_T2)/60; Times <- c(Times,IterTime)
        print(noquote(paste0("Var Selected: ",paste(Result$var,collapse = " ")," - Iteration lasted: ",
                             floor(IterTime)," mins ",round((IterTime-floor(IterTime))*60)," secs")))
      }
    }
  }
  vec <- c(if(c(rep("TRUE",length(example))[example == 1],"FALSE")[1]){rep(1,length(lambda))},if(c(rep("TRUE",length(example))[example == 2],"FALSE")[1]){rep(2,length(lambda))},
           if(c(rep("TRUE",length(example))[example == 3],"FALSE")[1]){rep(3,length(lambda))},if(c(rep("TRUE",length(example))[example == 4],"FALSE")[1]){rep(4,length(lambda))},
           if(c(rep("TRUE",length(example))[example == 5],"FALSE")[1]){rep(5,length(lambda))},if(c(rep("TRUE",length(example))[example == 6],"FALSE")[1]){rep(6,length(lambda))},
           if(c(rep("TRUE",length(example))[example == 7],"FALSE")[1]){rep(7,length(lambda))},if(c(rep("TRUE",length(example))[example == 8],"FALSE")[1]){rep(8,length(lambda))},
           if(c(rep("TRUE",length(example))[example == 9],"FALSE")[1]){rep(9,length(lambda))},if(c(rep("TRUE",length(example))[example == 10],"FALSE")[1]){rep(10,length(lambda))},
           rep(lambda,length(example)),rep("SKIMMS",nrows),apply(TP, 1, mean),
           apply(FP, 1, mean),apply(FN, 1, mean),apply(PE, 1, mean))
  
  Table <- matrix(vec,nrow = nrows); colnames(Table) <- c("Example","lambda","method","TP","FP","FN","PE")
  
  end_T <- proc.time()[[3]]; TotTime <- (end_T - start_T)/60
  print(noquote(paste0("Simulation Time Spent: ",floor(TotTime)," mins ",round((TotTime-floor(TotTime))*60)," secs")))
  
  param <- c(lsize = length(lambda), exsize = length(example),sim = sim, n = n,p = p,rho = rho, lambda = lambda, eps = eps, h = h,
             example = example, comment = comment)
  
  return(list(Table = Table, time = TotTime, param = param, PE = PE))
}

    ###########################
    ####### Simulations #######
    ###########################

ResultSim <- simulation(sim = 1, n = 100, p = 150, lambda = 0.01, eps = 0.01, example = 1, comment = FALSE, seed = 24)   

ResultSimAlt <- simulation(sim = 1, n = 100, p = 150, lambda = c(0.001,0.005,0.01,0.05,seq(0.1,0.5,0.1)), eps = 0.001, example = 1, comment = FALSE, seed = 24)
ResultSimAlt2 <- simulation(sim = 1, n = 100, p = 150, lambda = c(0.001,0.005,0.01,0.05,seq(0.1,0.5,0.1)), eps = 0.005, example = 1, comment = FALSE, seed = 24)
ResultSimAlt3 <- simulation(sim = 1, n = 100, p = 150, lambda = c(0.001,0.005,0.01,0.05,seq(0.1,0.5,0.1)), eps = 0.01, example = 1, comment = FALSE, seed = 24)
ResultSimAlt4 <- simulation(sim = 1, n = 100, p = 150, lambda = c(0.001,0.005,0.01,0.05,seq(0.1,0.5,0.1)), eps = 0.05, example = 1, comment = FALSE, seed = 24)

printResult(ResultSim)
printResult(ResultSimAlt)
printResult(ResultSimAlt2)
printResult(ResultSimAlt3)
printResult(ResultSimAlt4)

# saveRDS(ResultSimAlt, file = "ResultSimAlt.rds") # eps = 0.001 / lambda = c(0.001,0.005,0.01,0.05,seq(0.1,0.5,0.1))
# saveRDS(ResultSimAlt2, file = "ResultSimAlt2.rds") # eps = 0.005 / ||
# saveRDS(ResultSimAlt3, file = "ResultSimAlt3.rds") # eps = 0.01 /  ||
# saveRDS(ResultSimAlt4, file = "ResultSimAlt4.rds") # eps = 0.05 /  ||

printResult(MainSim) # Ex 1-3 # n = 100 , p = 150 , eps = 0.01 , h = 0.2 , sim = 20 , lambda = c(0.1,0.01,0.05,0.001)
printResult(MainSim2) # Ex 5-7 # n = 100 , p = 150 , eps = 0.01 , h = 0.2 , sim = 8 , lambda = c(0.5,0.1,0.01,0.05,0.001)
printResult(MainSim3) # Ex 8 # n = 100 , p = 150 , eps = 0.01 , h = 0.2 , sim = 8 , lambda = c(0.5,0.1,0.01,0.05,0.001)
printResult(MainSim4) # Ex 1-8 # n = 100 , p = 150 , eps = 0.001 , h = 0.2 , sim = 4 , lambda = c(0.5,0.1,0.01,0.05,0.001)

MainSim6 <- ResultSim

saveRDS(MainSim6, file = "MainSim6.rds")
   
   #############################
   #### Presentation Monday ####
   #############################

#1)
figure_step3 <- tau_optim_plot(truecoef = c(3,4), dimtooptim = c(1,2))
figure_step3_bis <- tau_optim_plot(truecoef = c(3,4), dimtooptim = c(2,3))
figure_step3
figure_step3_bis
# saveRDS(figure_step3, file = "FigureStep3.rds")

#2)
figureboptim
# saveRDS(figureboptim, file = "FigurebOptim.rds")


