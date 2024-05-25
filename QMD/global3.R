suppressPackageStartupMessages({

   library(moderndive)
   library(knitr)
   library(equatiomatic)
   library(data.table)
   library(dplyr)
   library(ggplot2)
   library(ggExtra)
   library(ggtext)
   library(sampler)
   library(kableExtra)
   library(ggstatsplot)
   library(formatR)
   library(moments)
   library(latex2exp)
   library(skedastic)
   library(cocor)
   library(olsrr)
   library(MASS)
   library(psych)
   library(Rfast)
   library(lme4)
   library(jtools)
   library(gap)
   library(lmPerm)
   library(performance)
   library(lmtest)
   library(nortest)
   library(boot)
   library(relaimpo)
   library(gtsummary)
   library(fastDummies)
   library(basicPlotteR)
   library(doBy)
   library(modelsummary)
   library(boot)
   library(foreign)
   library(xlsx)
   library(misty)
   library(ggside)
   library(varhandle)
   library(rcompanion)
   library(Hmisc)
   library(simstudy)
   library(QuantPsyc)
   library(infer)
   library(downloadthis)
   library(polycor)
   library(car)
   library(DescTools)
   library(survival)
   library(polycor)
   library(diagram)
   library(MASS)
   library(graphics)
   library(summarytools)
   library(xtable)
   library(stargazer)
   library(plot3D)
   library(ggplot2)
   library(ggpubr)
   library(relaimpo)
   library(GoodmanKruskal)
   library(DescTools)
   library(webexercises)

})


# From: https://r-coder.com/normal-distribution-r/
# mean: mean of the Normal variable
# sd: standard deviation of the Normal variable
# lb: lower bound of the area
# ub: upper bound of the area
# acolor: color of the area
# ...: additional arguments to be passed to lines function

normal_area <- function(mean = 0, sd = 1, lb, ub, acolor = "lightgray", ...) {
  x <- seq(mean - 3 * sd, mean + 3 * sd, length = 100) 
  
  if (missing(lb)) {
    lb <- min(x)
  }
  if (missing(ub)) {
    ub <- max(x)
  }
  
  x2 <- seq(lb, ub, length = 100)    
  plot(x, dnorm(x, mean, sd), type = "n", ylab = "")
  
  y <- dnorm(x2, mean, sd)
  polygon(c(lb, x2, ub), c(0, y, 0), col = acolor)
  lines(x, dnorm(x, mean, sd), type = "l", ...)
}

# Plot scattergram from fit <- lm()

ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}


# Calcul du nombre de concordances et de discordances
ConDisPairs <- function(x,y){
  n <- length(x)
  ix <- order(x)
  x <- x[ix]
  y <- y[ix]
  Nc <- sum(sapply(1:(n-1),function(i) sum(x[i]<x[(i+1):n] & y[i]<y[(i+1):n])))
  Nd <- sum(sapply(1:(n-1),function(i) sum(x[i]<x[(i+1):n] & y[i]>y[(i+1):n])))
  Nt <- sum(sapply(1:(n-1),function(i) sum(x[i]==x[(i+1):n] & y[i]==y[(i+1):n])))
  return(list(Nc=Nc, Nd=Nd, Nt=Nt))
}

#'Graph Sampling Distribution
#'
#'Graph Sampling Distribution
#'
#'@param xb vector of sampling statistics
#'@param nc scalar value: Confidence level (default = 0.95)
#'@param parm character value: name of the estimated parameter (default = 'Mean')
#'@examples
#'xb <- rnorm(1000,100,2)
#'GraphDist(xb,nc=0.99,parm='Medians')
#'@export

GraphDist <- function(xb, nc = 0.95, parm = "Moyennes") {
  intc <- quantile(xb,
                   c(0.5 * (1 - nc), nc + 0.5 * (1 - nc)))  # Calcul des percentiles
  gr1 <- hist(xb, breaks = "FD", plot = FALSE)
  cuts <- cut(gr1$breaks, c(-Inf, intc[1], intc[2], Inf))
  maintit <- paste("Distribution des", parm)
  plot(gr1,
       main = maintit,
       xlab = parm,
       col = c("red", "cadetblue1", "red")[cuts])
  txt1 <- paste("Moyenne    = ", round(mean(xb), 3))
  txt2 <- paste("Err.Std.   = ", round(sd(xb), 3))
  legend("topleft",
         legend = c(txt1, txt2),
         lty = 0:0,
         cex = 0.8,
         bg = "lightgreen")
  txt3 <- paste("IC(", nc * 100, "%) = [", round(intc[1], 3), ",", round(intc[2], 3), "]")
  legend("topright",
         legend = c(txt3),
         lty = 0:0, cex =
           0.8,
         bg = "cadetblue1")
}


#' Intervalle de confiance d'une corrélation - Simulation
#' Échantillons tirés d'une distribution normale
#' @param r = Corrélation au niveau de la population
#' @param n = Taille de l'échantillon
#' @param nrep = Nombre de ré-échantillonnages
#' @param nconf = Niveau de confiance
#' @export

R_CI <- function(r=0.36, n=50, nrep=100000, conf=95){

z   = .5*log((1+r)/(1-r))
se  = 1/sqrt(n-3)
zvec= rnorm(nrep)*se+z
ivec= (exp(2*zvec)-1)/(exp(2*zvec)+1)
low = (1-conf/100)/2
upp = ((1-conf/100)/2)+(conf/100)
LL  = quantile(ivec,low)
UL  = quantile(ivec,upp)
LL4 = format(LL,digits=3)
UL4 = format(UL,digits=3)
# LE=(exp(2*(.5*log((1+r)/(1-r))+se*qnorm(low,0,1)))-1)/(exp(2*(.5*log((1+r)/(1-r))+se*qnorm(low,0,1)))+1)
# UE=(exp(2*(.5*log((1+r)/(1-r))-se*qnorm(low,0,1)))-1)/(exp(2*(.5*log((1+r)/(1-r))-se*qnorm(low,0,1)))+1)
# LE4=format(LE,digits=3)
# UE4=format(UE,digits=3)

gr1 <- hist(ivec,
            breaks = "FD",
            plot = FALSE)
cuts <- cut(gr1$breaks, c(-Inf, LL, UL, Inf))
plot(gr1,
     main = "Distribution Empirique de r",
     col = c("red", "lightblue", "red")[cuts],
     freq=FALSE,
     xlab=paste("P[", LL4,"\u2264 r \u2264 ",UL4,"] = ",conf, "%"), 
     ylab = "Densité",
     cex=0.4)

#  (Empir) \n C[",LE4,"\u2264 \u03C1 \u2264 ",UE4, "] = ",conf, "% (Théor)
abline(v=mean(ivec), lty=2, lwd=2, col="#CC0000")
txt1 <- paste("Moyenne = ", round(mean(ivec), 3))
txt2 <- paste("Err.Std = ", round(sd(ivec), 3))
legend(ifelse(mean(ivec)>=0, "topleft", "topright"),
       legend = c(txt1, txt2),
       lty = 0:0,
       cex = 0.8,
       bg = "lightyellow")
}

#' Display matrix from R variable
#' From: https://stackoverflow.com/questions/45591286/for-r-markdown-how-do-i-display-a-matrix-from-r-variable
#' 
write_matex2 <- function(x) {
  begin <- "\\begin{bmatrix}"
  end <- "\\end{bmatrix}"
  X <-
    apply(x, 1, function(x) {
      paste(
        paste(x, collapse = "&"),
        "\\\\"
      )
    })
  paste(c(begin, X, end), collapse = "")
}



#' Fisher Transformation: r to z
#' Argument: Correlation
#' Output: Fisher's z
#' @param r = correlation
#' @export
r2z <- function(r){
  if((sum(r<(-1) | r>1))!=0) stop("Corrélation invalide!")
  z <- 0.5 * log((1 + r)/(1 - r))
  return(z)
}

#' Fisher Transformation: z to r
#' Argument: z (can be a vector)
#' Output: r
#' @param z = Fisher's z
#' @export
z2r <- function(z){
  if((sum(z<(-10) | z>10))!=0) stop("Z invalide!")
  r <- (exp(2 * z) - 1)/(1 + exp(2 * z))
  return(r)
}

#' Intervalle de confiance: corrélation transformée en Z (Fisher) - Simulation
#' Échantillons tirés d'une distribution normale
#' @param r = Corrélation au niveau de la population
#' @param n = Taille de l'échantillon
#' @param nrep = Nombre de réplications
#' @param nconf = Niveau de confiance
#' @export

Z_CI <- function(r=0.36, n=50, nrep=100000, conf=95){

  z=.5*log((1+r)/(1-r))
  se=1/sqrt(n-3)
  zvec=rnorm(nrep)*se+z
  low=(1-conf/100)/2
  upp=((1-conf/100)/2)+(conf/100)
  LL=quantile(zvec,low)
  UL=quantile(zvec,upp)
  LL4=format(LL,digits=4)
  UL4=format(UL,digits=4)

  gr1 <- hist(zvec,
              breaks = "FD",
              plot = FALSE)
  cuts <- cut(gr1$breaks, c(-Inf, LL, UL, Inf))
  plot(gr1,
       main = "Distribution Empirique de Zxy",
       col = c("red", "lightblue", "red")[cuts],
       freq=FALSE,
       xlab=paste("P[", LL4,"\u2264 \u005A \u2264 ",UL4,"] = ",conf, "%"), cex=0.4)
  abline(v=mean(zvec), lty=2, lwd=2, col="#CC0000")
  txt1 <- paste("Moyenne = ", round(mean(zvec), 3))
  txt2 <- paste("Err.Std = ", round(sd(zvec), 3))
  legend("topleft",
         legend = c(txt1, txt2),
         lty = 0:0,
         cex = 0.8,
         bg = "lightyellow")
  return(zst <- c(mean(zvec), sd(zvec)))
}


# Plot scattergram with marginals
#
# x <- rnorm(1000)
# y <- x*0.80+rnorm(1000)*sqrt(1-0.8**2)
# x <- x * 15 + 100
# y <- y * 6 + 20
# r <- cor(x,y)

PlotScat <- function(x, y){
  require(ggplot2)
  require(ggExtra)
  
df <- data.frame(x, y)
reg <- lm(y~x, data=df)
b0 <- reg$coefficients[1]
b1 <- reg$coefficients[2]
eq = paste0("Y' = ", round(reg$coeff[2],4), " x X ", round(reg$coeff[1],4))
p <- ggplot(df, aes(x, y)) + 
  geom_point() + 
  theme_classic() + 
  labs(x="X") +
  labs(y="Y") +
  geom_abline(intercept = reg$coeff[1], slope = reg$coeff[2], color="red", linewidth=1.3) + 
  ggtitle(eq) +
  geom_hline(yintercept = mean(y), linetype="dashed", color="blue", linewidth=1.2) +
  geom_vline(xintercept = mean(x), linetype="dashed", color="blue", linewidth=1.2) +
  geom_point(aes(x=mean(x), y=mean(y)), colour="yellow", size=5) +
  geom_text(x=2.6, y=-0.3, label="Moyenne Y") +
  geom_text(x=0, y=3.2, label="Moyenne X") +
  geom_text(x=-2, y=3, label=paste("R = ", round(r,3)), size=6)
ggExtra::ggMarginal(p, type = "densigram", col="black", fill="#7FFFD0")
}


# Fonction RegSim: Simulation de données bivariées
#
RegSim <- function(nrep=10000, n=50, xb=20, sx=5, yb=40, sy=10, r=0.6){
  
  b0 <- vector()
  b1 <- vector()
  rxy <- vector()
  
  for(i in 1:nrep) {
    x1 = rnorm(n)
    y1 = r * x1 + sqrt(1-r^2)*rnorm(n)
    x <- x1*sx + xb
    y <- y1*sy + yb
    reg <- lm(y~x)
    rxy[i] <- cor(x,y)
    b0[i] <- reg$coefficients[1]
    b1[i] <- reg$coefficients[2]
  }
  return(list(rxy=rxy, b0=b0, b1=b1))
}



# Tests d'hypothèse concernant des corrélations
CorTest <- function(r, n){
  if(length(r)!=length(n)){
    stop:"Les arguments doivent être de même longueur!"}
  if(length(r)==1){
    stop:"Les arguments doivent être des vecteurs de longueur supérieure à 1!"}
  z <- r2z(r)
  if(length(r)==2){
    Z <- (z[1]-z[2]) / sqrt(sum(1/(n-3)))
    p <- 2* (1 - pnorm(Z)) # Test bilatéral
  }else{
    chi2 <- sum((n-3) * z^2) - (sum((n-3) * z)^2 / sum(n-3))
    dl <- length(r) - 1      
    p <- pchisq(chi2, dl, lower.tail=FALSE)
  }
  if(length(r)==2){
    cat("Z = ", Z, "   p = ", p)
  }else{
    cat("Chi2 = ", chi2, "   p = ", p)
  }
}


#' @title PRESS
#' @author Thomas Hopper
#' @description Returns the PRESS statistic (predictive residual sum of squares).
#'              Useful for evaluating predictive power of regression models.
#' @param linear.model A linear regression model (class 'lm'). Required.
#' 
PRESS <- function(linear.model) {
  pr <- residuals(linear.model)/(1-lm.influence(linear.model)$hat)
  PRESS <- sum(pr^2)
  return(PRESS)
}


#' @title Predictive R-squared
#' @author Thomas Hopper
#' @description returns the predictive r-squared. Requires the function PRESS()
#' @param linear.model A linear regression model (class 'lm'). Required.
#'
pred_r2 <- function(linear.model) {
  lm.anova <- anova(linear.model)
  tss <- sum(lm.anova$'Sum Sq')
  r2pr <- 1-PRESS(linear.model)/(tss)
  return(r2pr)
}

ConfBand <- function(x, y, test=0, conf.lev=0.95){
  new <- data.frame(x=c(test))
  n <- length(x)
  cl <- (1-conf.lev)/2 + conf.lev
  reg <- lm(y ~ x)
  yp <- reg$fitted.values
  aov <- anova(reg)
  CMe <- aov$`Mean Sq`[2]
  IntBand <- predict(reg, newdata=new) + 
    c(-1, 1) * rep(qt(cl, n-2)) * 
    sqrt(CMe * (1/n + ((test-mean(x))^2)/ ((n-1)*var(x))))
  IntPred <- predict(reg, newdata=new) + 
    c(-1, 1) * rep(qt(cl, n-2)) * 
    sqrt(CMe * (1 + 1/n + ((test-mean(x))^2)/ ((n-1)*var(x))))
  return(list(IntBand, IntPred))
}

AdjustMeas <- function(VD, VI){
  if(!(length(VD)==length(VI))){stop:"Les vecteurs doivent être de mêmes longueurs!"}
  model <- lm(VD~VI)
  n <- length(VD)
  k <- 1
  res <-  model$residuals
  r2 <- cor(VD, VI)^2
  r2adj <- 1 - ((1-r2) * (n-1)) / (n - k - 1)
  press <- PRESS(model)
  r2pred <- pred_r2(model)
  mae <-  mean(abs(res))
  mse <-  mean((res)^2)
  rmse <- sqrt(mse)
  mape <- 100 * sum(abs(res / VD)) / n
  mpe <-  100 * sum(res / VD) / n
  out <- list(model, r2, r2adj, press, r2pred, mae, mse, rmse, mape, mpe)
  return(out)
}

# Functions to make SLR Confidence and Prediction Band Plots.
# Grant Brown, June 2013
# 171-162

# Internal function used by ConfidenceBandPlot and PredictionBandPlot
BasicRegressionPlot = function(X,Y, main = "", xlab = "",ylab = "", plt = TRUE, matchScale = FALSE)
{
  
  regression = lm(Y ~ X)
  if (plt)
  {
    if (matchScale)
    {
      lower = min(min(X), min(Y)) - 0.2*max(sd(X), sd(Y))
      upper = max(max(X), max(Y)) + 0.2*max(sd(X), sd(Y))
      plot(X,Y,main = main,xlab=xlab,ylab=ylab, xlim = c(lower, upper), ylim = c(lower, upper))
    }
    else
    {
      plot(X,Y,main = main,xlab=xlab,ylab=ylab, xlim = c(min(X)-sd(X), max(X) + sd(X)), ylim = c(min(Y)-sd(Y), max(Y) + sd(Y)))
    }
    abline(reg=regression)
  }
  return(regression)    
}
# Internal function used by ConfidenceBandPlot and PredictionBandPlot
BandPlot = function(X,Y,type = c("confidence", "prediction"),main="Regression + Confidence Band"
                    ,xlab="X",ylab="Y", col = "red", lty =2, lwd =1,
                    alpha=0.05, add = FALSE, matchScale = FALSE)
{
  
  
  if (! add)
  {
    regression = BasicRegressionPlot(X,Y,main,xlab,ylab,plt=TRUE,matchScale)
  }
  else
  {
    regression = BasicRegressionPlot(X,Y,plt=FALSE)
  }
  
  MSE = sum(regression$residuals^2)/regression$df.residual
  RSE = sqrt(MSE)
  
  n = length(Y)
  X0 = seq(min(X) - 10*sd(X), max(X) + 10*sd(X), length = 300)
  
  if (type == "confidence")
  {
    SE = RSE*sqrt(1/n + ((X0 - mean(X))^2)/(sum((X-mean(X))^2)))
  }
  else if (type == "prediction")
  {
    SE = RSE*sqrt(1 + 1/n + ((X0 - mean(X))^2)/(sum((X-mean(X))^2)))
  }
  else
  {
    print(paste("Unsupported type: ", type, sep = ""))
    print("Using Confidence Band")
  }
  
  LB = qt(alpha/2, df = n-2)*SE
  UB = qt(1-alpha/2,df=n-2)*SE
  Yhat = regression$coefficients[1] + regression$coefficients[2]*X0
  UB = UB + Yhat
  LB = LB + Yhat
  lines(X0, UB, lty=lty, col = col)
  lines(X0, LB, lty=lty, col = col)
  abline(h=0,lty=2,col="grey")
  abline(v=0,lty=2,col="grey")
}


# Make a Confidence Band Plot. The only required arguments are X and Y. 
ConfidenceBandPlot = function(X,Y,main="Regression + Confidence Band",xlab="X",
                              ylab="Y", col = "red", lty =2, lwd =1,matchScale = FALSE,
                              alpha=0.05, add = FALSE                             
)
{
  BandPlot(X,Y,type = "confidence", main = main, xlab =xlab,ylab=ylab,col=col,
           lty = lty, lwd = lwd, alpha = alpha, add = add, matchScale = matchScale)
}

# Make a Prediction Band Plot. The only required arguments are X and Y. 
PredictionBandPlot = function(X,Y,main="Regression + Prediction Band",xlab="X",
                              ylab="Y", col = "red", lty =2, lwd =1,matchScale = FALSE,
                              alpha=0.05, add = FALSE                             
)
{
  BandPlot(X,Y,type = "prediction", main = main, xlab =xlab,ylab=ylab,col=col,
           lty = lty, lwd = lwd, alpha =alpha, add = add, matchScale = matchScale)
}


MakeExamplePlots = function()
{
  # Simulate data for a regression
  X = seq(1,10, length = 100)
  Y = X + rnorm(length(X), 0, 5)
  # Open up a pdf device
  pdf(file = "Confidence and Prediction Bands.pdf")
  # Set up a grid of 2x2 plots
  par(mfrow = c(2,2))
  # Make a confidence band plot at alpha = 0.01
  ConfidenceBandPlot(X,Y, alpha = 0.01)
  
  # Make a prediction band plot at alpha = 0.01
  PredictionBandPlot(X,Y, alpha = 0.01)
  
  # Make a prediction + confidence band plot at alpha = 0.01
  ConfidenceBandPlot(X,Y, alpha = 0.01, main = "Regression + Confidence \n and Prediction Bands", lty = 1)
  PredictionBandPlot(X,Y, alpha = 0.01, add = TRUE, col = "blue")
  
  # Make a prediction + confidence band plot at alpha = 0.01 where the X and Y scale is identical.
  ConfidenceBandPlot(X,Y, alpha = 0.01, main = "Regression + Confidence \n and Prediction Bands", lty = 1, 
                     matchScale = TRUE)
  PredictionBandPlot(X,Y, alpha = 0.01, add = TRUE, col = "blue")
  # Write the plots to the pdf and close the file
  dev.off()
}


RegBandPlot <- function(x, y, type = c("confidence", "prediction", "both"), conflev = 0.95){
  options(warn=-1)
  if(length(x) != length(y)) {stop("X et Y sont de longueurs différentes!")
    return(NULL)}
  if(type == "confidence"){txt <- " de confiance"
  } else if (type == "prediction"){txt <- " de prédiction"
  } else {txt <- " de confiance et de prédiction"}
  
  df <- data.frame(x, y)
  model <- lm(y ~ x, data = df)
  
  predictions <- predict(model, interval = "prediction", level = conflev)
  confIntervl <- predict(model, interval = "confidence", level = conflev)
  
  df <- data.frame(cbind(df, predictions, confIntervl[, 2:3]))
  colnames(df) <- c("x", "y", "fit", "lwrPred", "uprPred", "lwrConf", "uprConf")
  
  gg <- ggplot(df, aes(x = x, y = y)) +
    geom_point(color = "blue") +
    geom_line(aes(y = predict(model, newdata = df)), 
              color = "red", 
              linetype = "dashed") +
    labs(x = "X", 
         y = "Y", 
         title = paste("Droite de régression avec bandes \n", txt)) +
    theme_classic()
  
  if(type == "confidence" | type == "both"){
    gg <- gg + geom_ribbon(data = df, 
                           aes(x = x, ymin = lwrConf, ymax = uprConf), 
                           fill = "blue", 
                           alpha = 0.2) }
  
  if(type == "prediction" | type == "both"){
    gg <- gg + geom_ribbon(data = df, 
                           aes(x = x, ymin = lwrPred, ymax = uprPred), 
                           fill = "green", 
                           alpha = 0.2) }
  
  print(gg)
#  out <- list(data = df, grf = gg)
#  return(out)
  
}

# Compute and print Confidence bands over scattergram and print table
RegInterval <- function(X, Y){
  # Créer un dataframe
  data <- data.frame(X = X, Y = Y)
  
  # Effectuer la régression linéaire
  reg <- lm(Y ~ X, data = data)
  
  # Calculer les prédictions et intervalles de confiance
  confid <- predict(reg, interval = "confidence")
  predictions <- suppressWarnings(predict(reg, interval = "predict"))
  
  # Créer un dataframe avec les résultats
  ConfTbl <- cbind(data, as.data.frame(predictions), confid[, 2:3])
  colnames(ConfTbl) <- c("X", "Y", "Y_Pred", "Pred_Inf", "Pred_Sup", "Conf_Inf", "Conf_Sup")
  
  # Créer le diagramme de dispersion avec les bandes de confiance et intervalles de prédiction
  p <- ggplot(ConfTbl, aes(x = X, y = Y)) + 
    geom_point() + 
    stat_smooth(method = lm) + 
    geom_line(aes(y = Pred_Inf), col = "coral2", linetype = "dashed") + 
    geom_line(aes(y = Pred_Sup), col = "coral2", linetype = "dashed") +
    ggtitle("Bande de Confiance (95%) et Intervalle de Prédiction")
  
  # Afficher le diagramme
  print(p)
  
  out = list(table = ConfTbl)
#   out = list(grf = p, table = ConfTbl)
  return(out)
  
}


#' Simpson's paradox dataset simulation
#'
#' @inheritParams bayestestR::simulate_correlation
#' @param n The number of observations for each group to be generated.
#' @param groups Number of groups.
#' @param difference Difference between groups.
#'
#' @return A dataset.
#'
#' @examples
#' data <- simulate_simpson(n = 100, groups = 5, r = 0.5)
#'
#' ggplot(data, aes(x = V1, y = V2)) +
#'   geom_point(aes(color = Group)) +
#'   geom_smooth(aes(color = Group), method = "lm") +
#'   geom_smooth(method = "lm")
#' @export
simulate_simpson <- function(n = 100, r = 0.5, groups = 5, difference = 1) {
  data <- data.frame()
  for (i in 1:groups) {
    dat <- bayestestR::simulate_correlation(n = n, r = r)
    dat$V1 <- dat$V1 + difference * i # (i * -sign(r))
    dat$V2 <- dat$V2 + difference * (i * -sign(r))
    dat$Group <- LETTERS[i]
    data <- rbind(data, dat)
  }
  return(data)
}

simulate_simpson2 <- function(n = 100, r = -0.5, groups = 5, difference = 1) {
  data <- data.frame()
  for (i in 1:groups) {
    dat <- bayestestR::simulate_correlation(n = n, r = runif(1, -0.7, 0.7))
    dat$V1 <- dat$V1 + difference * i # (i * -sign(r))
    dat$V2 <- dat$V2 + difference * (i * -sign(r))
    dat$Group <- LETTERS[i]
    data <- rbind(data, dat)
  }
  return(data)
}


# Calcul des valurs résiduelles 'studentisées' interne
ResStud <- function(model){
  if(!is(model, "lm")){stop: "Classe incorrecte: utilisez lm() au préalable."}
  res <- model$residuals
  CMerr <- anova(model)$`Mean Sq`
  res_stud <- res / (sqrt(CMerr[2] * (1 - hatvalues(model))))
  return(res_stud)
}

# Calcul d'une corrélation commune (groupes multiples)
CorComb <- function(r, n){
  if(length(r)!=length(n)){
    stop:"Les arguments doivent être de même longueur!"}
  if(length(r)==1){
    stop:"Les arguments doivent être des vecteurs de longueur supérieure à 1!"}
  zw <- sum(r2z(r) * (n - 3)) / sum(n - 3)
  rho = (exp(2 * zw)-1) / (exp(2 * zw)+1)
  cat("Rho = ", rho)
}

#' FROM: omyxjutsu library - https://rdrr.io/github/bryancquach/omixjutsu/
#' Histogram and jittered boxplots partitioned by categories
#'
#' Creates a histogram and jittered boxplot with color-coded, separate distributions for each value
#' of a user-specified categorical variable.
#'
#' Creates two plots that can be organized as a single column two row multi-figure. The top row
#' plot is multiple overlayed histograms from data partitioned by a variable. The bottom row plot
#' is multiple horizontal boxplots with jittered data points overlayed. The x-axis of both plots
#' are fixed to the same limits for comparability across plots.
#'
#' @param data A two-column data frame with numeric values for plotting in the first column and
#' the categorical variable as the second column.
#' @param binsize A numeric value for the histogram bin widths.
#' @param colors A string vector. The histogram bar colors and jitter colors for each categorical
#' variable value.
#' @param hist_alpha A numeric value for the alpha level for histogram bars.
#' @param box_fill A string. The fill color for the boxplot.
#' @param box_alpha A numeric value for the alpha level for the boxplot.
#' @param box_lwd A numeric value for the boxplot line width.
#' @param jitter_alpha A numeric for the alpha level for jittered data points.
#' @param jitter_size A numeric for the point size for jittered data points.
#' @param x_title A string denoting the x-axis title. Only added to the boxplot.
#' @param y_title A string denoting the y-axis title. Only added to the histogram.
#' @return A list of two ggplot objects `hist` and `boxplot`.
#' @seealso \code{\link{hist_boxplot}}
#' @export
hist_boxplot2 <- function(data,
                          binsize = diff(range(data[, 1], na.rm = T)) / 25,
                          colors = NULL,
                          hist_alpha = 0.75,
                          box_fill = "gray50",
                          box_alpha = 0.5,
                          box_lwd = 1,
                          jitter_alpha = 0.75,
                          jitter_size = 1.75,
                          x_title = "",
                          y_title = "") {
  colnames(data) <- c("value", "group")
  ggout_list <- list()
  data_min <- min(data$value, na.rm = T)
  data_max <- max(data$value, na.rm = T)
  axis_min <- switch(as.character(sign(data_min)),
                     "-1" = {
                       data_min * 1.1
                     },
                     "1" = {
                       data_min * 0.9
                     },
                     "0" = {
                       data_min - binsize
                     }
  )
  axis_max <- switch(as.character(sign(data_max)),
                     "-1" = {
                       data_max * 0.9
                     },
                     "1" = {
                       data_max * 1.1
                     },
                     "0" = {
                       data_max + binsize
                     }
  )
  ggout_list$hist <- ggplot(data, aes(x = value, fill = group)) +
    geom_histogram(
      position = "identity",
      binwidth = binsize,
      alpha = hist_alpha,
      color = "white",
    ) +
    xlim(axis_min, axis_max) +
    labs(x = x_title, y = y_title, fill = "") +
    theme(
      plot.margin = unit(c(0.5, 0.5, 0, 0.5), units = "cm"),
      plot.title = element_text(size = 18),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 18),
      axis.title.y = element_text(size = 18, vjust = 3),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16)
    )
  ggout_list$boxplot <- ggplot(data, aes(x = value, y = group, color = group)) +
    geom_boxplot(
      position = "identity",
      alpha = box_alpha,
      color = "black",
      fill = box_fill,
      outlier.alpha = 0,
      lwd = box_lwd
    ) +
    geom_jitter(
      shape = 16,
      size = jitter_size,
      alpha = jitter_alpha,
      height = 0.15
    ) +
    xlim(axis_min, axis_max) +
    labs(x = x_title, y = y_title, color = "") +
    theme(
      plot.margin = unit(c(0, 0.5, 0.5, 0.5), units = "cm"),
      plot.title = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 18),
      axis.title.x = element_text(size = 18, vjust = -1),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16)
    )
  if (is.null(colors)) {
    ggout_list$hist <- ggout_list$hist + scale_fill_brewer(palette = "Dark2")
    ggout_list$boxplot <- ggout_list$boxplot + scale_color_brewer(palette = "Dark2")
  } else {
    ggout_list$hist <- ggout_list$hist + scale_fill_manual(values = colors)
    ggout_list$boxplot <- ggout_list$boxplot + scale_color_manual(values = colors)
  }
  return(ggout_list)
}

hist_boxplot2 <- 
  function (x, freq = TRUE, density = FALSE, main = NULL, xlab = NULL, 
            ymax = NULL, ...) 
  {
    abs <- deparse(substitute(x))
    if (is.null(xlab)) {
      xlab = abs
    }
    ted <- hist(x, plot = FALSE)
    par(yaxs = "i")
    if (freq) {
      ylim = c(-(max(ted$counts)/5), max(ted$counts)) %>% round(1)
      boxwex = max(ted$counts)/6
    }
    else {
      if (density) {
        max = max(max(density(x)$y), max(ted$density)) %>% round(1)
      }
      else {
        max = max(ted$density) %>% round(1)
      }
      ylim = c(-max/5, max) %>% round(1)
      boxwex = max/6 %>% round(1)
    }
    if (!is.null(ymax)) 
      ylim <- c(-ymax/5, ymax)
    if (is.null(main)) {
      main = paste("Histogram of", abs)
    }
    hist(x, ylim = ylim, yaxt = "n", xlab = xlab, freq = freq, 
         main = main, ...)
    axis(2, at = seq(0, ylim[2], length = 5))
    usrs <- par()$usr
    par(new = TRUE, bty = "n", xaxs = "i", yaxs = "i", 
        ann = FALSE)
    boxplot(x, horizontal = TRUE, ylim = usrs[1:2], add = TRUE, 
            at = usrs[3]/2, pars = list(boxwex = 5), col = "lightblue", 
            boxwex = ylim[2]/6)
    if (density) {
      lines(density(x), lwd = 2, col = "red")
    }
  }

# Print upper or lower triangular square matrix
PrtCorMat <- function(R, tri = "lower"){
  require(kableExtra)
  if (!is.matrix(R)) 
    stop("L'argument R n'est pas une matrice")
  if(!(nrow(R) == ncol(R)))
    stop("La matrice R n'est pas carrée")
  
  if(!tri %in% c("lower", "upper")){tri = "lower"}
  
  cmdstr <- ifelse(tri == "lower", "upper.tri(R)", "lower.tri(R)")
  triangle <- R
  triangle[eval(parse(text = cmdstr))] <- ""
  triangle<-as.data.frame(triangle)
  triangle %>%
  kbl() %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), 
                full_width=FALSE) %>%
  column_spec(1, bold=TRUE, background="#89CFF0") %>%
  row_spec(0, bold=TRUE, background="#89CFF0") 
}


#PRESS - predicted residual sums of squares

Press <- function(linear.model) {
  #' Calcul de valeurs résiduelles prédictives
  pr <- residuals(linear.model)/(1-lm.influence(linear.model)$hat)
  Press <- sum(pr^2)
  return(Press)
}

pred_r_squared <- function(linear.model) {
  #' Calcul des sommes de carrés
  lm.anova <- anova(linear.model)
  #' Somme des carrés totale
  tss <- sum(lm.anova$'Sum Sq')
  # Calcul du R2 prédit
  pred.r.squared <- 1-Press(linear.model)/(tss)
  return(pred.r.squared)
}

rmspe <- function(DF){
  y <- DF[["Sales"]]
  y_obs <- y[DF[["Type"]] == "obser"]
  y_est <- y[DF[["Type"]] == "Forecast"]
  sqrt(mean(((y_obs - y_est)/y_obs)^2))
}