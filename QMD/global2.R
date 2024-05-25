suppressMessages({
	library(moments)
	library(formatR)
	library(kableExtra)
	library(sampler)
	library(dplyr)
	library(agrmt)
	library(latex2exp)
	library(ggplot2)
	library(gridExtra)
	library(fitdistrplus)
	library(MASS)
	library(qwraps2)
	library(data.table)
	library(interpretCI)
	library(glue)
	library(wBoot)
	library(reshape2)
	library(PASWR2)
	library(boot)
	library(rcompanion)
	library(pwr)
	library(moments)
})

#' From OpenIntroOrg
#' https://rdrr.io/github/OpenIntroOrg/oilabs-r-package/
#' Plot confidence intervals
#' 
#' @description This function takes vectors of upper and lower confidence 
#' interval bounds, as well as poplutation mean, and plots which confidence 
#' intervals intersect the mean
#' @param lo a vector of lower bounds for confidence intervals
#' @param hi a vector of upper bounds for confidence intervals
#' @param m the population mean
#' @export
#' 

plot_ci <- function(lo, hi, m) {
  par(mar=c(2, 1, 1, 1), mgp=c(2.7, 0.7, 0))
  k <- length(lo)
  ci.max <- max(rowSums(matrix(c(-1*lo,hi),ncol=2)))
  
  xR <- m + ci.max*c(-1, 1)
  yR <- c(0, 41*k/40)
  
  plot(xR, yR, type='n', xlab='', ylab='', axes=FALSE)
  abline(v=m, lty=2, col='#00000088')
  axis(1, at=m, paste("mu = ",round(m,4)), cex.axis=1.15)
  #axis(2)
  for(i in 1:k){
    x <- mean(c(hi[i],lo[i]))
    ci <- c(lo[i],hi[i])
    if((m < hi[i] & m > lo[i])==FALSE){
      col <- "#F05133"
      points(x, i, cex=1.4, col=col)
      #		  points(x, i, pch=20, cex=1.2, col=col)
      lines(ci, rep(i, 2), col=col, lwd=5)
    } else{
      col <- 1
      points(x, i, pch=20, cex=1.2, col=col)
      lines(ci, rep(i, 2), col=col)
    }
  }
}

#' Function normal_area(): Tracé d'une distribution normale avez zone colorée
#'
#' Arguments:
#' mean: Moyenne de la distribution
#' sd: Écart-type de la distribution
#' lb: Borne inférieure de la zone à colorer
#' ub: Borne supérieure de la zone à colorer
#' acolor: Couleur de la zone
#' ...: Arguments additionnels à passer à la fonction lines()
#'
#' @exemple
#' normal_area(mean=100, sd=15, lb=100, ub=125, acolor="red")
#' @export
#'


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

#' Function t_area(): draw a Student distribution with shaded area
#'
#' Arguments:
#' df: Degrees of freedom
#' lb: lower bound of the area
#' ub: upper bound of the area
#' acolor: color of the area
#' ...: additional arguments to be passed to lines function
#'
#' @example
#' t_area(df=25, lb=0, ub=1.25, acolor="red")
#' @export
#'

t_area <- function(df=25, lb, ub, acolor = "lightgray", ...) {
  x <- seq(-5, 5, length = 100)
  
  if (missing(lb)) {
    lb <- min(x)
  }
  if (missing(ub)) {
    ub <- max(x)
  }
  
  x2 <- seq(lb, ub, length = 100)
  plot(x, dt(x, df), type = "n", ylab = "")
  
  y <- dt(x2, df)
  polygon(c(lb, x2, ub), c(0, y, 0), col = acolor)
  lines(x, dt(x, df), type = "l", ...)
  lines(x,rep(0,length(x)))
}

#'Graph Sampling Distribution
#'
#'@param xb vector of sampling statistics
#'@param nc scalar value: Confidence level (default = 0.95)
#'@param parm character value: name of the estimated parameter (default = 'Mean')
#'@examples
#'xb <- rnorm(1000,100,2)
#'GraphDist(xb,nc=0.99,parm='Medians')
#'@export
GraphDist <- function(xb, nc = 0.95, parm = "Means") {
  intc <- quantile(xb,
                   c(0.5 * (1 - nc), nc + 0.5 * (1 - nc)))  # Calcul des percentiles
  gr1 <- hist(xb, breaks = "FD", plot = FALSE)
  cuts <- cut(gr1$breaks, c(-Inf, intc[1], intc[2], Inf))
  maintit <- paste("Distribution of", parm)
  plot(gr1,
       main = maintit,
       xlab = parm,
       col = c("red", "cadetblue1", "red")[cuts])
  txt1 <- paste("Mean    = ", round(mean(xb), 3))
  txt2 <- paste("StdErr = ", round(sd(xb), 3))
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


hist_boot <- function(data, shade = "over", threshold1 = NULL, threshold2 = NULL) {
  if (shade == "over") {
    ggplot(data, aes(x = X)) +
      geom_histogram(binwidth = 0.05, fill = "cadetblue1", color = "black") +
      geom_histogram(data = subset(data, X > threshold1), binwidth = 0.05, 
                     colour = "black", fill = "red") +
      theme_bw()
  } else if (shade == "below") {
    ggplot(data, aes(x = X)) +
      geom_histogram(binwidth = 0.05, fill = "cadetblue1", color = "black") +
      geom_histogram(data = subset(data, X < threshold1), binwidth = 0.05, 
                     colour = "black", fill = "red") +
      theme_bw()
  } else if (shade == "between") {
    ggplot(data, aes(x = X)) +
      geom_histogram(binwidth = 0.05, fill = "cadetblue1", color = "black") +
      geom_histogram(data = subset(data, X > threshold1 & X < threshold2), binwidth = 0.05, 
                     colour = "black", fill = "red") +
      theme_bw()
  } else {
    print("Invalid 'shade' argument. Choose from 'over', 'below', or 'between'.")
  }
}


plot_histogram <- function(x, shade = "over", threshold1 = NULL, threshold2 = NULL, obs = NULL) {
  data <- data.frame(X = x)
  
  n_intervals <- nclass.FD(x)
  breaks <- seq(min(x), max(x), length.out = n_intervals + 1)
  hist_data <- hist(x, 
                    breaks = breaks, 
                    plot = FALSE)
  MaxCnt <- max(hist_data$counts)
  
  p <- ggplot(data, aes(x = X)) +
    geom_histogram(bins = n_intervals, fill = "lightgrey", color = "black") +
    theme_bw() +
    labs(title = "Distribution de Ré-échantillonnage des moyennes") +
    theme(plot.title = element_text(hjust = 0.5))
  
  if (shade == "over") {
    p <- p +
      geom_histogram(data = subset(data, X > threshold2), bins = n_intervals, 
                     colour = "black", fill = "red") 
    
  } else if (shade == "below") {
    p <- p +
      geom_histogram(data = subset(data, X < threshold1), bins = n_intervals, 
                     colour = "black", fill = "red") 
    
  } else if (shade == "between") {
    p <- p +
      geom_histogram(data = subset(data, X > threshold2 | X < threshold1), bins = n_intervals, 
                     colour = "black", fill = "red") 
    
  } else {
    print("Invalid 'shade' argument. Choose from 'over', 'below', or 'between'.")
    return(NULL)
  }
  
  # Calculate mean and standard deviation
  mean_val <- mean(x)
  sd_val <- sd(x)
  qtl <- quantile(x, c(0.025, 0.05, 0.95, 0.975))
  IC <- if(shade == "over"){paste("$C[", round(qtl[3], 3), "\\leq \\mu \\leq \\infty] = 0.95$")
  } else if(shade == "below"){paste("$C[-\\infty \\leq \\mu \\leq ", round(qtl[2], 3), "] = 0.95$")
  } else {paste("$C[", round(qtl[1], 3), "\\leq \\mu \\leq ", round(qtl[4], 3), "] = 0.95$")}
  
  labxpos <- max(x) - 0.35 * diff(range(x)) 
  labypos <- 0.85 * MaxCnt 
  
  p <- p + xlab(expression(bar(X))) +
    annotate("text", x = mean_val, y = 0, label = TeX("$\\mu_{\\bar{X}}$"),
             vjust = 1, hjust = 0.5, color = "red", size = 5) +
    geom_segment(aes(x = mean_val, xend = mean_val, y = 0, yend = MaxCnt),
                 linetype = "solid", color = "red",linewidth = 0.8) +
    annotate("text", x = labxpos, y = labypos, parse = TRUE,
             label = as.character(TeX(paste("$\\mu_{\\bar{x}} =", round(mean_val, 2), "$"))), 
             hjust = 0, size = 4, color = "darkgreen", fontface = "bold") +
    annotate("text", x = labxpos, y = labypos - labypos/10, parse = TRUE,
             label = as.character(TeX(paste("$\\sigma_{\\bar{x}} =", round(sd_val, 2), "$"))), 
             hjust = 0, size = 4, color = "darkgreen", fontface = "bold") +
    annotate("text", x = labxpos, y = labypos - labypos/5, parse = TRUE,
             label = as.character(TeX(IC)), 
             hjust = 0, size = 4, color = "darkgreen", fontface = "bold") 
  if(!is.null(obs)){
    p <- p +
        geom_segment(aes(x = obs, xend = obs, y = 0, yend = MaxCnt * 0.1),
                     linetype = "solid", 
                     color = "darkblue",
                     linewidth = 1.3,
                    arrow = arrow(length = unit(0.4, "cm"), 
                                  ends = "first")) +
    
        annotate("text", x = obs, y = MaxCnt * 0.1, label = TeX("$\\bar{X}_{obs}$"),
                 vjust = -0.5, hjust = 0.5, color = "darkblue", size = 6) }
  return(p)
}


# Fonction boot.mean()
# Arguments: data = vecteur des données échantillonnales
#            num = Nombre de ré-échantillonnages à effectuer
# Sortie: std.err = Erreur-Standard de la distribution d'échantillonnage
#         resamples = liste des num échantillons
#         means = vecteur des moyennes observées pour chaque échantillon
boot.mean <- function(data, num) {
  resamples <- lapply(1:num, function(i) sample(data, replace=TRUE))
  r.mean <- sapply(resamples, mean) # Calcul des num moyennes échantillonnales
  std.err <- sd(r.mean)  # calcul de l'erreur-standard des moyennes
  list(std.err=std.err, resamples=resamples, means=r.mean)   
}


#' Calcul de puissance: Test d'hypothèse sur une moyenne unique
#' Arguments:
#' @param n scalaire: taille de l'échantillon (défaut: 30)
#' @param alpha scalaire: niveau de signification (défaut: 0.05)
#' @param H0 scalaire: valeur postulée sous H0 (défaut:0)
#' @param H1 scalaire: valeur proposée sous H1 (défaut: 2)
#' @param sig scalaire: Écart-Type de la population (défaut: 1)
#' @param lat scalaire [1, 2]: spécification d'un test uni- ou bilatéral
#' Sortie:  puissance
#' @export
PwrCalc <- function(n = 30,
                    alpha = .05,
                    H0 = 0,
                    H1 = 2,
                    sig = 1,
                    lat = 1) {
  if(!lat %in% c(1, 2)){stop: "Argument LAT doit être 1 ou 2 (Uni/Bilatéral"}
  # Déterminer la valeur critique (Z)
  z_alpha <- qnorm(p = alpha/lat,
                   mean = 0,
                   sd = 1,
                   lower.tail = FALSE)
  # Déterminer la moyenne critique
  y_bar_cv1 <- z_alpha*(sig/sqrt(n)) + H0
  if(lat==2){
    y_bar_cv2 <- -z_alpha*(sig/sqrt(n)) + H0}
  
  # Calcul de la puissance
  power <- pnorm(q = y_bar_cv1,
                 mean = H1,
                 sd = sig/sqrt(n),
                 lower.tail = FALSE)
  if(lat==2){
    power <- power + pnorm(q = y_bar_cv2,
                           mean = H1,
                           sd = sig/sqrt(n),
                           lower.tail = TRUE)
  }
  # Afficher la puissance
  power
}

#'
#'Sample from a uniform distribution, with random lower and upper bounds
#'@param n scalar value: Number of observations to draw (default = 50)
#'@param a scalar value: lower bound of distribution
#'@param b scalar value: upper bound of distribution
#'@return list of vector of observations, and lower/upper bounds (quantitative values)
#'@examples
#'x1 <- RandUnif() # use defaults
#'x2 <- RandUnif(25,c(1,5),c(10,15)) # draw 25 observations within
#'      lower limit set randomly from 1-5, and the upper limit set
#'      randomly from 10-15
#'@export
RandUnif <- function(n = 25, a = 10, b = 20) {
  if(a>=b){stop("Le paramètre a doit être inférieur au paramètre b!")}
  liminf <- runif(1, a, b) 
  limsup <- runif(1, liminf, b)  
  x <- runif(n, liminf, limsup)
  out <- list(x,liminf,limsup)
  return(out)
}

# Fonctions utilisées pour calculer le rapport de variation

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

RV = function(x){  
  mod = getmode(x)
  fm = sum(x==mod)
  rv = 1-fm/length(x)
  
}

# Calcul de la différence de RV (2 groupes)
rvDiff = function(x,y){  
  rv1 = RV(x)
  rv2 = RV(y)
  DiffRV = rv1 - rv2
  return(list(rv1, rv2, DiffRV))
}


#' Distribution de Student
t.Distr <- function(t, df, alpha, h1){
  x <- seq(-4,4,by=0.01)
  y <- dt(x, df)
  if(h1!=0){
    right <- qt(1-alpha, df)
    left <- qt(alpha, df)
  }else{
    right <- qt(1-alpha/2, df)
    left <- qt(alpha/2, df)
  }
  plot(x, y, type="l",
       xaxt="n",ylab="Probabilité",
       xlab=expression(paste('Distribution des  ',(bar(X) - bar(Y)),' supposant  ', H[0],'  vraie')),
       axes=FALSE,
       ylim=c(0,max(y)*1.05),
       xlim=c(min(x),max(x)),
       frame.plot=FALSE)
  axis(1,at=c(round(left, 3), 0, round(right, 3), round(t, 3)),
       pos = c(0,0),
       labels=TRUE)
  axis(2)
  xReject <- c(seq(-4,4,by=0.01))
  yReject <- dt(xReject, df)
  if(h1==-1){
    yReject[xReject > left ] <- 0}else{
      if(h1==1){yReject[xReject < right ] <- 0}else{
        if(h1==0){yReject[xReject > left & xReject < right] <- 0}
      }
    }
  polygon(c(xReject,xReject[length(xReject)],xReject[1]),
          c(0, yReject , 0), col='lightblue')
  if(h1!=0){
    p <- pt(abs(t), df, lower.tail=FALSE)
  }else{
    p <- 2 * pt(abs(t), df, lower.tail=FALSE)
  }
  legend("topleft", paste("p = ", round(p, 4)))
  arrows(t, .15, t, 0, length=0.1, col="red", lwd=2)
  text(t, 0.17, "t calculé", cex=1.2)
}

#' Calcul du nombre de concordances/discordances/égalités
#' dans un tableau de contingence (2 variables ordinales)
#' @param x Tableau de contingence bidimensionnel
#' x est un tableau (matrice, array, table ou data frame) bidimensionnel
#' @examples
#' x <- table(sample(1:3, 50, replace=TRUE), sample(1:5, 50, replace=TRUE) )
#' ConDiscord(x)
#' @export
#'
ConDiscord <- function(x){
  if(!is.array(x) &
     !is.matrix(x) &
     !is.data.frame(x) &
     !is.table(x)){stop:"L'argument présenté n'est pas du type requis!"}
  if(!length(dim(x)==2)){stop:"Tableau de contingence bidimensionnel requis!"}
  r_ndx = row(x)
  c_ndx = col(x)
  Nc <- sum(x * mapply(function(r, c){sum(x[(r_ndx > r) & (c_ndx > c)])}, r = r_ndx, c = c_ndx))
  Nd <- sum(x * mapply(function(r, c){sum(x[(r_ndx > r) & (c_ndx < c)])}, r = r_ndx, c = c_ndx) )
  Ntr <- sum(x * mapply(function(r, c){sum(x[(r_ndx == r)])}, r = r_ndx, c = c_ndx) )
  Ntc <- sum(x * mapply( function(r, c){sum(x[(c_ndx == c)])},r = r_ndx, c = c_ndx) )
  return(list(Nc=Nc, Nd=Nd, Ntr=Ntr, Ntc=Ntc))
}


#' Chi-squared test statistic for contingency tables
#'
#' Calculates the chi-squared test statistic for a two-way contingency table.
#' @param tab A \code{K x C} matrix (contingency table) of counts. See details.
#' @return The calculated value of the chi-squared statistic.
#' @details Suppose that \code{tab} consists of counts from \eqn{K} populations (rows) in \eqn{C} categories.  The chi-squared test statistic is computed as
#' \deqn{
#'   \sum_{i=1}^K \sum_{j=1}^C (E_{ij} - O_{ij})^2/E_{ij},
#' }{
#'   \sum_ij (E_ij - O_ij)^2/E_ij,
#' }
#' where \eqn{O_{ij}}{O_ij} is the observed number of counts in the \eqn{i}th row and \eqn{j}th column of \code{tab}, and \eqn{E_{ij}}{E_ij} is the expected number of counts under \eqn{H_0} that the populations have indentical proportions in each category:
#' \deqn{
#'   E_{ij} = \frac 1 N \sum_{i=1}^K O_{ij} \times \sum_{j=1}^C O_{ij}.
#' }{
#'   E_ij = \sum_i O_ij * \sum_j O_ij / N,
#' }
#' where \eqn{N} is the total number of counts in \code{tab}.
#' @examples
#' # simple contingency table
#' ctab <- rbind(pop1 = c(5, 3, 0, 3),
#'                 pop2 = c(4, 10, 2, 5))
#' colnames(ctab) <- LETTERS[1:4]
#' ctab
#' chi2.stat(ctab) # chi^2 test statistic
#' @export
chi2.stat <- function(tab) {
  N <- rowSums(tab) # number of fish in each lake
  p0 <- colSums(tab)/sum(N) # MLE of common probability vector
  tabE <- N %o% p0 # expected
  sum(((tab-tabE)^2/tabE)[tabE > 0]) # chi2
}

#' somersD: Calcul du coefficient D de Somers
#'
#' Arguments:
#' table: Tableau de contingence bidimensionnel
#' conflev: Niveau de confiance, 0.95 par défaut
#' digits: Nombre de décimales, 3 par défaut
#'
#' Sortie: liste comprenant:
#'           - Tableau contenant pour VD=R et VD=C,
#'               - D de Somers
#'               - Erreur Standard de D
#'               - Limite inférieure de l'IC
#'               - Limite supérieure de l'IC
#'           - Nombre de concordances
#'           - Nombre de discordances
#'
#' @examples
#' x <- matrix(c(12, 4, 3, 5, 10, 6, 3, 5, 14),nrow=3, byrow=TRUE)
#' somersD(x)
#' @export
#'

somersD = function(tbl, conflev=0.95, digits=3){
  if(!length(dim(tbl))==2){stop("Tableau de contingence bidimensionnel nécessaire")}
  N <- sum(tbl)
  nr <- nrow(tbl)
  nc <- ncol(tbl)
  r_ndx = row(tbl)
  c_ndx = col(tbl)
  C <- sum(tbl * mapply(function(r, c){sum(tbl[(r_ndx > r) & (c_ndx > c)])},
                        r = r_ndx, c = c_ndx))
  D <- sum(tbl * mapply( function(r, c){sum(tbl[(r_ndx > r) & (c_ndx < c)])},
                         r = r_ndx, c = c_ndx) )
  
  Dc.r <- 2 * (C - D) / (N^2 - sum(rowSums(tbl)^2))
  Dr.c <- 2 * (C - D) / (N^2 - sum(colSums(tbl)^2))
  
  # Asymptotic standard error: sqrt(sigma2)
  sigma2r <- (4 * (nr^2-1) * (nc + 1)) / (9 * N * nr^2 * (nc - 1))
  sigma2c <- (4 * (nc^2-1) * (nr + 1)) / (9 * N * nc^2 * (nr - 1))
  SEDc <- sqrt(sigma2c)
  SEDr <- sqrt(sigma2r)
  
  pr2 <- 1 - (1 - conflev)/2
  CIr <- qnorm(pr2) * SEDr * c(-1, 1) + Dr.c
  CIc <- qnorm(pr2) * SEDc * c(-1, 1) + Dc.r
  
  ResR <- c(Dr.c, SEDr, CIr)
  ResC <- c(Dc.r, SEDc, CIc)
  
  res <- data.frame(rbind(ResR, ResC))
  colnames(res) <- c("Somers'D", "StdErr.D", "Lower", "Upper")
  rownames(res) <- c("Rows", "Columns")
  
  
  return(list(ResTbl=res, concord=C, discord=D))
  
}

#' Display matrix from R variable
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

sim2x2 <- function(p1=.1, p2=.45, n1=200, N=10000){
  set.seed(NULL)
  or=rep(NA,N)
  lnor=rep(NA,N)
  a=rep(NA,N)
  b=rep(NA,N)
  c=rep(NA,N)
  d=rep(NA,N)
  rr=rep(NA,N)
  lrr=rep(NA,N)
  selnor=rep(NA,N)
  llor=rep(NA,N)
  ulor=rep(NA,N)
  selnrr=rep(NA,N)
  llrr=rep(NA,N)
  ulrr=rep(NA,N)
  #
  for(i in 1:N) {
#   n1=as.integer(runif(1,300,600))
    x1 = runif(n1,0,1)<p1
    x2 = runif(n1,0,1)<p2
    a[i]=sum((x1==1)&(x2==1))
    b[i]=sum((x1==0)&(x2==1))
    c[i]=sum((x1==1)&(x2==0))
    d[i]=sum((x1==0)&(x2==0))
# Relative Risk
    rr[i]=(a[i] / (a[i]+b[i])) / (c[i] / (c[i]+d[i]))
    lrr[i]=log(rr[i])
    selnrr[i]=sqrt((b[i]/a[i])/(a[i]+b[i])+(d[i]/c[i])/(c[i]+d[i]))
    llrr[i]=lrr[i]-1.96*selnrr[i]
    ulrr[i]=lrr[i]+1.96*selnrr[i]
# Odds Ratio
    or[i]=(b[i]*c[i])/(a[i]*d[i])
    lnor[i]=log(or[i])
    selnor[i]=sqrt(1/a[i]+1/b[i]+1/c[i]+1/d[i])
    llor[i]=lnor[i]-1.96*selnor[i]
    ulor[i]=lnor[i]+1.96*selnor[i]
  }
  dataout=data.frame(a,b,c,d,rr, lrr, selnrr,llrr, ulrr,or,lnor,selnor,llor,ulor)
}



#' Fonction DistHist(): Production d'histogrammes comportants des zones colorées
#' Arguments:
#'   x: Vecteur de données, de moyennes, ou de toute autre statistique
#'   lb: borne inférieure. Les valeurs inférieures à lb seront colorées
#'   ub: borne supérieure. Les valeurs supérieures à ub seront colorées
#'   stat: statistique considérée dans le cas d'une distribution d'échantillonnage. On peut spécifier une chaîne utilisant expression() ou TeX()
#'   title: Chaîne de caractères utilisée en tant que titre du graphique
#'   col: Couleur de la zone spécifiée par lb et ub
#'   xlim: Vecteur composé des limites de l'axe horizontal
#'   
#'   Note: En sauvegardant le graphique dans une variable, on peut ajouter des éléments au graphique
#'   Exemples d'utilisation: 
#'        DistHist(rnorm(1000, 100, 15), stat = TeX("$\\bar{X}$"), title = TeX("Distribution des $\\bar{X}$"))
#'
#'        data <- rnorm(100000, mean = 10, sd = 2)
#'        DistHist_ggplot(data, lb = 6, ub = 13, 
#'                          stat = expression(bar(X)), 
#'                          title = "Distribution d'échantillonnage des Moyennes", 
#'                          col = "lightblue")

DistHist2 <- function(x, 
                     lb = NULL, 
                     ub = NULL, 
                     stat = expression(bar(X)), 
                     title = "Distribution d'échantillonnage des Moyennes", 
                     col = "lightblue",
                     xlim = NULL){
  
  h <- hist(x,
            breaks = "FD",
            plot = FALSE)
  df <- data.frame(mids = h$mids, counts = h$counts)
  
  if (!is.null(lb) | !is.null(ub)) {
    if (!is.null(lb) & !is.null(ub)) {
      df$clr <- ifelse(df$mids < ub & df$mids > lb, "grey", col)
    } else if (!is.null(lb) & is.null(ub)) {
      df$clr <- ifelse(df$mids > lb, "grey", col)
    } else if (is.null(lb) & !is.null(ub)) {
      df$clr <- ifelse(df$mids < ub, "grey", col)
    }
  } else {
    df$clr <- col
  }
  
  p <- ggplot(df, aes(x = mids, y = counts, fill = clr)) +
    geom_bar(colour = "black", stat = "identity") +
    labs(title = title, x = stat, y = "Fréquence") +
    scale_fill_identity() +
    #   scale_fill_manual(values = c("grey" = "grey", "color" = col)) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 16), 
      axis.title.y = element_text(size = 16),  
      plot.title = element_text(size = 20, hjust = 0.5),
      axis.line = element_line(color = "black"),  
      axis.text = element_text(color = "black")) +
    guides(fill = "none") 
  
  # Applying xlim if provided
  if (!is.null(xlim)) {
    p <- p + coord_cartesian(xlim = xlim)
  }
  
  return(p)
}


#' Fonction DistHist(): Production d'histogrammes comportants des zones colorées
#' Arguments:
#'   x: Vecteur de données, de moyennes, ou de toute autre statistique
#'   lb: borne inférieure. Les valeurs inférieures à lb seront colorées
#'   ub: borne supérieure. Les valeurs supérieures à ub seront colorées
#'   stat: statistique considérée dans le cas d'une distribution d'échantillonnage. On peut spécifier une chaîne utilisant expression() ou TeX()
#'   title: Chaîne de caractères utilisée en tant que titre du graphique
#'   col: Couleur de la zone spécifiée par lb et ub
#'   xlim: Vecteur composé des limites de l'axe horizontal
#'   threshold: permet de colorer une zone intérieure plutôt que les extrémités
#'   
#'   Note: En sauvegardant le graphique dans une variable, on peut ajouter des éléments au graphique
#'   Exemples d'utilisation: 
#'        DistHist(rnorm(1000, 100, 15), stat = TeX("$\\bar{X}$"), title = TeX("Distribution des $\\bar{X}$"))
#'
#'        data <- rnorm(100000, mean = 10, sd = 2)
#'        DistHist(data, lb = 6, ub = 13, 
#'                 stat = expression(bar(X)), 
#'                 title = "Distribution d'échantillonnage des Moyennes", 
#'                 col = "lightblue")
DistHist <- function(x, 
                     lb = NULL, 
                     ub = NULL, 
                     stat = expression(bar(X)), 
                     title = "Distribution d'échantillonnage des Moyennes", 
                     col = "lightblue",
                     xlim = NULL,
                     threshold = NULL){
  
  h <- hist(x,
            breaks = "FD",
            plot = FALSE)
  num_bins <- nclass.FD(x)
  df <- data.frame(x = x)
  
  if (!is.null(threshold)) {
    df$clr <- ifelse(df$x > threshold, "grey", col)
  } else {
    df$clr <- col
  }
  
  if (!is.null(lb) & !is.null(ub)) {
    df$clr <- ifelse(df$x > lb & df$x < ub, "grey", df$clr)
  } else if (!is.null(lb) & is.null(ub)) {
    df$clr <- ifelse(df$x > lb, "grey", df$clr)
  } else if (is.null(lb) & !is.null(ub)) {
    df$clr <- ifelse(df$x < ub, "grey", df$clr)
  }
  
  p <- ggplot(df, aes(x = x, fill = clr)) +
    geom_histogram(binwidth = diff(range(x)) / num_bins, colour = "black", position = "identity") +
    labs(title = title, x = stat, y = "Frequency") +
    scale_fill_identity() +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 16), 
      axis.title.y = element_text(size = 16),  
      plot.title = element_text(size = 20, hjust = 0.5),
      axis.line = element_line(color = "black"),  
      axis.text = element_text(color = "black")) +
    guides(fill = "none") 
  
  if (!is.null(xlim)) {
    p <- p + coord_cartesian(xlim = xlim)
  }
  
  return(p)
}

# Ré-échantillonnage: moyenne unique
# Arguments: 
# - x : Observations
# - k : Nombre de réplications 
# - mu : Moyenne de la population (option)  
# Sortie:
# - means : moyennes observées pour chaque ré-échantillon
# - mean : moyenne des moyennes
# - se : écart-type des moyennes (erreur-standard)
# - meansCtr: moyennes centrées, si mu est spécifié

resamp.mean <- function(x, k = 1000, mu = NULL){
  m <- replicate(k,
                 mean(sample(x, 
                             length(x), 
                             replace = TRUE)))
  muxb <- mean(m)
  sexb <- sd(m)
  mC <- if(!is.null(mu)){m - muxb + mu}
  return(list("means" = m, "mean" = muxb, "se" = sexb, "meansCtr" = mC))
}


#' Show Chi2-value table
#'
#' Show Chi2-value table
#' @param DF Numeric degree of freedom
#' @param p Numeric probability
#' @param alternative Character One of c("two.sided","greater","less")
#' @importFrom flextable flextable highlight autofit
#' @importFrom purrr map_dfc
#' @return An object of class "flextable"
#' @export
#' @examples
#' show_Chi2_table()
show_Chi2_table=function(DF=20, p=0.05, alternative="greater"){
  DF=round(DF)
  if(alternative=="two.sided") p=p/2
  pp=c(0.999, 0.995, 0.99, 0.975, 0.95, 0.9, 0.1, 0.05, 0.025, 0.01, 0.005, 0.001)
  pp=unique(c(pp,p))
  pp=sort(pp,decreasing=TRUE)
  df=max(1,DF-2):(max(1,DF-2)+3)
  suppressMessages(res<-purrr::map_dfr(df,function(x){
    purrr::map_dfc(pp,~sprintf("%.3f",qchisq(.,x)))
  }))
  colnames(res)=pp
  res$df=df
  myrow=which(res$df==DF)
  res<-res %>% select(.data$df,everything())
  
  
  
  tbl <- flextable::flextable(res)%>%
    flextable::align(i=1,align="center",part="header") %>%
    flextable::align(align="center",part="body") %>%
    flextable::highlight(i=myrow,j=grep(p,colnames(res))) %>%
    flextable::autofit() 
  flextable::set_caption(tbl, paste("Critical values of chi-square, for df = ", DF, " |  Direction: ", alternative))
}

#' Adapted From: interpretCI package
#' Draw normal distribution curve
#' @param mean vector of means
#' @param sd vector of standard deviations
#' @param z vector of quantiles
#' @param p vector of probabilities
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of "two.sided" (default), "greater" or "less".
#' @importFrom stats dnorm qnorm
#' @return A ggplot
#' @export
#' @examples
#' draw_n()
#' draw_n(alternative="less")
#' draw_n(alternative="greater")
#' draw_n(z=-1.75)
#' draw_n(z=-1.75,alternative="greater")
#' draw_n(z=-1.75,alternative="less")
draw_n=function(mean=0,
                sd=1,
                z=NULL,
                p=0.05,
                alternative="two.sided"){
  # mean=0;sd=1;z=100;alternative="greater"
  if(!is.null(z)) {
    if(alternative=="two.sided"){
      p=pnorm(-abs(z))*2
    } else if(alternative!="less"){
      p=pnorm(z,lower.tail=FALSE)
    } else{
      p=pnorm(z)
    }
    showp=TRUE
  } else{
    showp=FALSE
  }
  x=c(seq(qnorm(p=0.0001,
                mean=mean,
                sd=sd),
          qnorm(p=0.9999,
                mean=mean,
                sd=sd),
          length=100))
  y=dnorm(x,mean=mean,sd=sd)
  data=data.frame(x=x,y=y)
  data=rbind(data,data[100,],data[100,])
  data$y[101:102]=0
  data$x[102]=data$x[1]
  data
  
  if(alternative=="two.sided"){
    if(near(p,0)){
      x1=-Inf
      x2=Inf
    } else if(near(p,1)){
      x1=-Inf
      x2=Inf
    } else{
      x1=seq(qnorm(p=0.0001,mean=mean,sd=sd),qnorm(p=p/2,mean=mean,sd=sd),length=50)
      x2=seq(qnorm(p=1-p/2,mean=mean,sd=sd),qnorm(p=0.9999,mean=mean,sd=sd),length=50)
    }
    y1=dnorm(x1,mean=mean,sd=sd)
    y2=dnorm(x2,mean=mean,sd=sd)
    data2=data.frame(x=x1,y=y1)
    if(nrow(data2)==50){
      data2=rbind(data2,data2[50,],data2[50,])
      data2$y[51:52]=0
      data2$x[52]=data2$x[1]
    }
    
    
    data3=data.frame(x=x2,y=y2)
    if(nrow(data3)==50){
      data3=rbind(data3,data3[50,],data3[50,])
      data3$y[51:52]=0
      data3$x[52]=data3$x[1]
    }
  } else if(alternative=="less"){
    if(near(p,1)){
      x1=seq(qnorm(p=0.0001,
                   mean=mean,
                   sd=sd),
             qnorm(p=0.9999,
                   mean=mean,
                   sd=sd),
             length=50)
      x1[50]=-Inf
    } else if(near(p,0)){
      x1=Inf
    }else{
      x1=seq(qnorm(p=0.0001,
                   mean=mean,
                   sd=sd),
             qnorm(p=p,
                   mean=mean,
                   sd=sd),
             length=50)
    }
    y1=dnorm(x1,mean=mean,sd=sd)
    data2=data.frame(x=x1,y=y1)
    if(nrow(data2)==50){
      data2=rbind(data2,data2[50,],data2[50,])
      data2$y[51:52]=0
      data2$x[52]=data2$x[1]
    }
  } else{
    if(near(p,1)){
      x2=seq(qnorm(p=0.0001,
                   mean=mean,
                   sd=sd),
             qnorm(p=0.9999,
                   mean=mean,
                   sd=sd),
             length=50)
      x2[1]=-Inf
    } else if(near(p,0)){
      x2=Inf
    } else{
      x2=seq(qnorm(p=1-p,
                   mean=mean,
                   sd=sd),
             qnorm(p=0.9999,
                   mean=mean,
                   sd=sd),
             length=50)
    }
    y2=dnorm(x2,mean=mean,sd=sd)
    data3=data.frame(x=x2,y=y2)
    if(nrow(data3)==50){
      data3=rbind(data3,data3[50,],data3[50,])
      data3$y[51:52]=0
      data3$x[52]=data3$x[1]
    }
  }
  
  if(alternative=="two.sided"){
    if(showp){
      labels=paste0("p=",rep(round(p/2,3),3))
    } else{
      labels=paste0("z=",c(round(qnorm(p=p/2,mean=mean,sd=sd),3),round(qnorm(p=1-p/2,mean=mean,sd=sd),3)))
    }
  } else{
    if(showp){
      labels=paste0("p=",rep(round(p,3),3))
    } else{
      labels=paste0("z=",c(round(qnorm(p=p,mean=mean,sd=sd),3),round(qnorm(p=1-p,mean=mean,sd=sd),3)))
    }
  }
  p1=ggplot(data,aes(x=x,y=y))+
    geom_line()
  
  if(alternative!="greater"){
    p1=p1+
      geom_polygon(data=data2,aes(x=x,y=y),fill="red",alpha=0.5)+
      annotate(geom="text",x=data2$x[50],y=data2$y[50],label=labels[1],hjust=ifelse(p>0.5,-0.3,1.3),vjust=-0.5)
  }
  if(alternative!="less"){
    p1=p1+ geom_polygon(data=data3,aes(x=x,y=y),fill="red",alpha=0.5)+
      annotate(geom="text",x=data3$x[1],y=data3$y[1],label=labels[2],hjust=ifelse(p>0.5,1.3,-0.3),vjust=-0.5)
  }
  
  p1 <- p1 + 
#      annotate(geom="text",x=0,y=0.1,label=paste0("N(",mean,",",sd,")"))+
#      annotate(geom="text",x=0,y=0.08,label=paste0("p=",round(p,3)))+
#      annotate(geom="text",x=0,y=0.06,label=paste0("alternative=",alternative))+
      labs(x="Score Z",y="Densité")+
      theme_classic()
  return(p1)
}


#' Adapted from interpretCI package
#' Draw t distribution curve
#' @param DF numeric degree of freedom
#' @param t numeric t value
#' @param p numeric p value
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of "two.sided" (default), "greater" or "less".
#' @importFrom dplyr near
#' @return A ggplot
#' @export
#' @examples
#' draw_t(DF=30)
#' draw_t(DF=20,t=2.5)
#' draw_t(DF=49,t=1.77)
#' draw_t(DF=49,p=0.005)
#' draw_t(DF=19,t=-0.894,alternative="less")
#' draw_t(DF=146,t=0.67,alternative="greater")
draw_t=function(DF=50,t=NULL,p=0.05,alternative="two.sided"){
  # DF=19;t=2;alternative="greater"
  if(!is.null(t)) {
    if(alternative=="two.sided") {
      p=pt(-abs(t),DF)
      p=p*2
    } else if(alternative=="less"){
      p=pt(t,DF)
    } else{
      p=pt(t,DF,lower.tail=FALSE)
      p
    }
    
    showp=TRUE
  } else{
    showp=FALSE
  }
  x=seq(qt(p=0.0001,df=DF),qt(p=0.9999,df=DF),length=100)
  if(near(p,1)){
    x[1]=-Inf
  } else if(near(p,0)){
    x[100]=Inf
  }
  
  y=dt(x,df=DF)
  data=data.frame(x=x,y=y)
  data=rbind(data,data[100,],data[100,])
  data$y[101:102]=0
  data$x[102]=data$x[1]
  
  
  if(alternative=="two.sided"){
    
    if(near(p,0)){
      x1=-Inf
      x2=Inf
    } else if(near(p,1)){
      x1=-Inf
      x2=Inf
    } else{
      x1=seq(qt(p=0.0001,df=DF),qt(p=p/2,df=DF),length=50)
      x2=seq(qt(p=1-p/2,df=DF),qt(p=0.9999,df=DF),length=50)
    }
    y1=dt(x1,df=DF)
    y2=dt(x2,df=DF)
    data2=data.frame(x=x1,y=y1)
    if(nrow(data2)==50){
      data2=rbind(data2,data2[50,],data2[50,])
      data2$y[51:52]=0
      data2$x[52]=data2$x[1]
    }
    data3=data.frame(x=x2,y=y2)
    if(nrow(data3)==50){
      data3=rbind(data3,data3[50,],data3[50,])
      data3$y[51:52]=0
      data3$x[52]=data3$x[1]
    }
    
  } else if(alternative=="less"){
    if(near(p,1)){
      x1=seq(qt(p=0.0001,df=DF),qt(p=0.9999,df=DF),length=50)
      x1[50]=-Inf
    } else if(near(p,0)){
      x1=Inf
    }else{
      x1=seq(min(qt(p=0.0001,df=DF),qt(p=p,df=DF)),qt(p=p,df=DF),length=50)
    }
    y1=dt(x1,df=DF)
    data2=data.frame(x=x1,y=y1)
    if(nrow(data2)==50){
      data2=rbind(data2,data2[50,],data2[50,])
      data2$y[51:52]=0
      data2$x[52]=data2$x[1]
    }
  } else{
    
    if(near(p,1)){
      x2=seq(qt(p=0.0001,df=DF),qt(p=0.9999,df=DF),length=50)
      x2[1]=-Inf
    } else if(near(p,0)){
      x2=Inf
    } else{
      x2=seq(qt(p=1-p,df=DF),max(qt(p=1-p,df=DF),qt(p=0.9999,df=DF)),length=50)
    }
    y2=dt(x2,df=DF)
    data3=data.frame(x=x2,y=y2)
    if(nrow(data3)==50){
      data3=rbind(data3,data3[50,],data3[50,])
      data3$y[51:52]=0
      data3$x[52]=data3$x[1]
    }
    
  }
  
  if(alternative=="two.sided"){
    if(showp){
      labels=paste0("p=",rep(round(p/2,3),3))
    } else{
      labels=paste0("t=",c(round(qt(p=p/2,DF),3),round(qt(p=1-p/2,DF),3)))
    }
  } else{
    if(showp){
      labels=paste0("p=",rep(round(p,3),3))
    } else{
      labels=paste0("t=",c(round(qt(p=p,DF),3),round(qt(p=1-p,DF),3)))
    }
  }
  p1=ggplot(data,aes(x=x,y=y))+
    geom_line()
  
  if(alternative!="greater"){
    no=ifelse(nrow(data2)>50,50,nrow(data2))
    hjust=ifelse(is.infinite(data2$x[no]),-1,1.3)
    p1=p1+
      geom_polygon(data=data2,aes(x=x,y=y),fill="red",alpha=0.5)+
      annotate(geom="text",x=data2$x[no],y=data2$y[no],label=labels[1],hjust=hjust,vjust=-0.5)
    
  }
  if(alternative!="less"){
    hjust=ifelse(is.infinite(data3$x[1]),1.3,-0.3)
    p1=p1+ geom_polygon(data=data3,aes(x=x,y=y),fill="red",alpha=0.5)+
      annotate(geom="text",x=data3$x[1],y=data3$y[1],label=labels[2],hjust=hjust,vjust=-0.5)
  }
  
  p1<-p1+ annotate(geom="text",x=0,y=0.1,label=paste0("t(dl=",DF,")"))+
#    annotate(geom="text",x=0,y=0.08,label=paste0("p=",round(p,3)))+
#    annotate(geom="text",x=0,y=0.06,label=paste0("alternative=",alternative))+
    labs(x="Score t",y="Densité")+
    theme_classic()
  p1
  
}

Curv_Chisq <- function(x, df){
  curve(dchisq(x, df = df), 
        from = 40, to = 150,
        main = paste('Distribution du Chi-Carré (dl = ', df, ')'),
        ylab = 'Densité',
        lwd = 2)
  x_lower95 <- seq(0, x)
  p_lower95 <- dchisq(x_lower95, df = df)
  polygon(c(x_lower95, rev(x_lower95)), c(p_lower95, rep(0, length(p_lower95))),
          col = adjustcolor('red', alpha=0.3), border = NA)
}


PwrPlot <- function(mu0=0,
                    mu1=1.9,
                    sigma=5,
                    nsig=0.05,
                    n=25,
                    tail=1){
  alph <- nsig
  nalph <- 1 - alph
  bet <- pnorm(qnorm(nalph, mu0, sigma/sqrt(n)), mu1, sigma/sqrt(n))
  power <- 1 - bet
  sig0 <- sigma/sqrt(n)
  liminf <- qnorm(0.001, mu0, sig0)
  limsup <- qnorm(0.999, mu1, sig0)
  ggplot(data.frame(x = seq(liminf, limsup, length.out=100)), aes(x)) +
    stat_function(fun = dnorm, args = c(mu0, sig0), geom = 'area',
                  xlim = c(qnorm(1-alph/tail, mu0, sig0), limsup),
                  aes(fill = 'red')) +
    stat_function(fun = dnorm, args = c(mu0, sig0)) +
    stat_function(fun = dnorm, args = c(mu0, sig0), geom = 'area',
                  xlim = c(liminf, qnorm(1-alph/tail, mu0, sig0)),
                  aes(fill = 'deepskyblue3')) +
    stat_function(fun = dnorm, args = c(mu0, sig0)) +
    stat_function(fun = dnorm, args = c(mu1, sig0)) +
    stat_function(fun = dnorm, args = c(mu1, sig0), geom = 'area',
                  xlim = c(qnorm(1-alph/tail, mu0, sig0), limsup),
                  aes(fill = 'green'), alpha=0.2) +
    stat_function(fun = dnorm, args = c(mu1, sig0), geom = 'area',
                  xlim = c(liminf, qnorm(1-alph/tail, mu0, sig0)),
                  aes(fill = 'chocolate3'), alpha=0.5) +
    scale_fill_identity(
      labels = c(paste("\u03B2=", round(bet, 3)),
                 paste("(1-\u03B1)=", nalph),
                 paste("(1-\u03B2)=", round(power, 3)),
                 paste("\u03B1=", alph)),
      guide = guide_legend(title="Probabilités", label.hjust = 0)) +
    theme(legend.title=element_text(size=15)) +
    annotate("text",
             size = 5,
             x = c(qnorm(0.05, mu0, sig0), qnorm(0.95, mu1, sig0)),
             y = c(dnorm(qnorm(0.25, mu0, sig0), mu0, sig0), 
                   dnorm(qnorm(0.75, mu1, sig0), mu1, sig0)),
             label = c(expression(H[0], H[1])),
             colour = rep("black", 2)) +
    labs(y="Densité", x="Test") +
    theme(legend.position="right",
          legend.text = element_text(margin = margin(r = 10, unit = "pt")))+
    theme(legend.text = element_text(size = 10))
  
}


#' Calcul de puissance: Test d'hypothèse sur une moyenne unique
#' Arguments:
#' @param n scalaire: taille de l'échantillon (défaut: 30)
#' @param alpha scalaire: niveau de signification (défaut: 0.05)
#' @param H0 scalaire: valeur postulée sous H0 (défaut:0)
#' @param H1 scalaire: valeur proposée sous H1 (défaut: 2)
#' @param sig scalaire: Écart-Type de la population (défaut: 1)
#' @param lat scalaire [1, 2]: spécification d'un test uni- ou bilatéral
#' Sortie:  puissance
#' @export
PwrCalc <- function(n = 30,
                    alpha = .05,
                    H0 = 0,
                    H1 = 2,
                    sig = 1,
                    lat = 1) {
  if(!lat %in% c(1, 2)){stop: "Argument LAT doit être 1 ou 2 (Uni/Bilatéral"}
  # Déterminer la valeur critique (Z)
  z_alpha <- qnorm(p = alpha/lat,
                   mean = 0,
                   sd = 1,
                   lower.tail = FALSE)
  # Déterminer la moyenne critique
  y_bar_cv1 <- z_alpha*(sig/sqrt(n)) + H0
  if(lat==2){
    y_bar_cv2 <- -z_alpha*(sig/sqrt(n)) + H0}
  
  # Calcul de la puissance
  power <- pnorm(q = y_bar_cv1,
                 mean = H1,
                 sd = sig/sqrt(n),
                 lower.tail = FALSE)
  if(lat==2){
    power <- power + pnorm(q = y_bar_cv2,
                           mean = H1,
                           sd = sig/sqrt(n),
                           lower.tail = TRUE)
  }
  # Afficher la puissance
  power
}

# Adapté de wBoot::print.boot.two()
print.boot.2 <- function(x, ...)
{
  test <- !is.null(x$Null)
  
  obs_diff <- x$Observed
  xbd <- x$Boot.values
  qtl <- c(x$Confidence.limits)
  h <- hist(xbd, 
            breaks = "FD", 
            plot = FALSE)
  breaks <- h$breaks
  counts <- h$counts
  maxf <- max(h$counts)
  txt <- paste("$C\\left [", round(qtl[1], 3), "\\leq \\Delta \\leq ", round(qtl[2], 3), "\\right ]=95\\%$")
  txt1 <- paste("$\\mu_{\\Delta} =", round(x$Mean, 3))
  txt2 <- paste("$\\sigma_{\\Delta}=", round(x$SE, 3), "$")
  left_tail <- which(breaks >= min(xbd) & breaks < qtl[1])
  right_tail <- which(breaks > qtl[2] & breaks <= max(xbd))
  
  cols <- rep("lightblue", length(breaks))
  cols[left_tail] <- "red"
  cols[right_tail] <- "red"
  
  hist(xbd, 
       breaks = breaks, 
       col = cols, 
       border = "white", 
       main = "Distribution des Différences",
       xlab = NULL,
       sub = TeX(txt))
  
  points(x$Observed, 0, pch = 17, col = "orange", cex = 1.2)
  points(x$Null, 0, pch = 19, col = "darkgreen", cex = 1.2)
  abline(v = x$Mean, lty = 2, lwd = 1.5, col = "darkblue")
  legend("topright", legend = c("Int. Conf", "Zone de Rejet", "Diff. observée", "Moyenne", expression(mu[0])),
         bty = "n",
         lwd = 2, cex = 1.2, col = c("lightblue", "red", "orange", "darkblue", "darkgreen"), 
         lty = c(NA, NA, NA, 2, NA), pch = c(15, 15, 17, NA, 19))
  legend("topleft", legend = c(TeX(txt1), TeX(txt2)))
  
}


F.Distr <- function(f, df1, df2, alpha=.05){
  x <- seq(0,4,by=0.01)
  y <- df(x, df1, df2)
  right <- qf(1-alpha, df1, df2)
  plot(x, y, type="l",
       xaxt="n",ylab="Probabilité",
       xlab="Rapport F",
       axes=FALSE,
       ylim=c(0,max(y)*1.05),
       xlim=c(min(x),max(x)),
       frame.plot=FALSE)
  axis(1,at=c(0, round(right, 3), round(f, 3)),
       pos = c(0,0),
       labels=TRUE)
  axis(2)
  xReject <- c(seq(0,4,by=0.01))
  yReject <- df(xReject, df1, df2)
  yReject[xReject < right ] <- 0
  
  polygon(c(xReject, xReject[length(xReject)], xReject[1]),
          c(0, yReject , 0), col='lightblue')
  p <- pf(f, df1, df2, lower.tail=FALSE)
  legend("topleft", paste("p = ", round(p, 4)))
  arrows(f, .15, f, 0, length=0.1, col="red", lwd=2)
  text(f, 0.25, "F \n calculé", cex=1.2)
  arrows(right, .35, right, 0, length=0.1, col="red", lwd=2)
  text(right, 0.45, "F \n critique", cex=1.2)
}
