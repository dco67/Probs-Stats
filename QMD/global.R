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

#'
#'Generate a frequency table
#'@param x a numeric vector
#'@param bins number of intervals (bins)
#'@return data frame:
#'      Lower & Upper limits of intervals,
#'      Absolute frequencies,
#'      relative frequencies,
#'      cumulative frequencies,
#'      Cumulative/Relative frequencies
#'@examples
#'TblFreq(rnorm(200,50,8),n_classes(x))
#'@export
TblFreq <- function(x, bins = 20) {
  hg <- hist(x, bins)  # Obtenir les paramètres de l'histogramme
  FreqAbs <- hg$counts  # Fréquences absolues
  FreqRel <- prop.table(FreqAbs)  # Fréquences relatives
  FreqCum <- cumsum(FreqAbs)  # Fréquences cumulatives
  FreqRelCum <- cumsum(FreqRel)  # Fréquences Relatives et Cumulatives
  lim <- hg$breaks[1:(length(hg$breaks))] - 0.5  # Limites exactes des intervalles
  LimInf <- lim[1:(length(lim) - 1)]
  LimSup <- lim[2:(length(lim))]
  output <- data.frame(LimInf, LimSup, FreqAbs, FreqRel, FreqCum, FreqRelCum)
  return(output)
}

#'
#'Generate a frequency table with various graphs
#'@param x a numeric vector
#'@param bins number of intervals (bins)
#'@return data frame:
#'      Lower & Upper limits of intervals,
#'      Absolute frequencies,
#'      relative frequencies,
#'      cumulative frequencies,
#'      Cumulative/Relative frequencies
#'      graphs (histogram, boxplot, cumulative, ogive, frequency polygon...)
#'@examples
#'TblFreq2(rnorm(200,50,8),n_classes(x))
#'@export
TblFreq2 <- function(x, bins = 20) {
  # Fonction TblFreq2: x = vecteur de données; bins = nombre d'intervalles
  hg <- hist(x, bins, main = "Histogramme", xlab = "X")  # Obtenir les paramètres de l'histogramme
  MidPt <- hg$mids
  FreqAbs <- hg$counts  # Fréquences absolues
  FreqRel <- prop.table(FreqAbs)  # Fréquences relatives
  FreqCum <- cumsum(FreqAbs)  # Fréquences cumulatives
  FreqRelCum <- cumsum(FreqRel)  # Fréquences Relatives et Cumulatives
  lim <- hg$breaks[1:(length(hg$breaks))] - 0.5  # Limites exactes des intervalles
  LimInf <- lim[1:(length(lim) - 1)]
  LimSup <- lim[2:(length(lim))]
  hist(x, bins, main = "Histogramme et Polygone de Fréquences", xlab = "X")
  par(new = TRUE)
  # Polygone de fréquences
  lines(MidPt, FreqAbs, col = "red", lwd = 3)
  par(new = FALSE)
  plot(MidPt, FreqAbs, type = "l", main = "Polygone de Fréquences", xlab = "X", ylab = "Fréquence",
       col = "red", lwd = 3)
  boxplot(x, main = "Diagramme en Boîte", xlab = "X", range = 1.5, horizontal = TRUE, boxwex = 0.25)
  # Distribution Cumulative
  barplot(FreqCum, main = "Distribution de Fréquences Cumulatives", xlab = "X", ylab = "Fréquences Cumulatives",
          space = FALSE, names.arg = MidPt)
  plot(MidPt, FreqRelCum, type = "l", main = "Distribution des fréquences Relatives Cumulatives",
       xlab = "X", ylab = "Proportion", col = "red", lwd = 3)
  output <- data.frame(LimInf, LimSup, FreqAbs, FreqRel, FreqCum, FreqRelCum)
  return(output)
}


if (!require("pacman")) install.packages("pacman")
pacman::p_load(formatR, 
               randtests, 
               VennDiagram, 
               ggvenn, 
               arrangements, 
               tidyr, 
               openintro, 
               Rfast, 
               extraDistr, 
               kableExtra, 
               dplyr, 
               DT,
               gtools,
               knitr,
               bookdown,
               tidyverse,
               reticulate,
               ggplot2,
               scales,
               kableExtra,
               moments,
               psych,
               DescTools,
               wakefield,
               agrmt, 
               mnonr, 
               multimode, 
               modeest,
               latex2exp)

blank_theme <- theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold", hjust = 0.5)
  )

options(scipen = 999)

#' Permutations of n objects
#'
#' Generate a full list of permutations of n objects
#' @param x A vector containing n items (numeric/alphanumeric)
#' @return Array of permuted objects, one permutation per line
#' @examples
#' p1 <- permute(1:4)
#' p2 <- permute(letters[1:5])
#' @export
permute <- function(x) {
  if (length(x) == 1) {
    return(x)
  } else {
    res <- matrix(nrow = 0, ncol = length(x))
    for (i in seq_along(x)) {
      res <- rbind(res, cbind(x[i], Recall(x[-i])))
    }
    return(res)
  }
}

  #' Dice Roll
  #'
  #'@param n number of Dice to roll (integer value, default = 1)
  #'@param p vector of probabilities of a 'success' (6 probs), default = 6X(1/6)
  #'       if it is a single value, the probabilities are equals
  #'@return Value of each die
  #'@export
  DiceRoll <- function(n = 1, p = rep(1, 6)) {
    if (!testInt(n)) {
      stop("n must be an integer!")
    }
    if (length(p) == 1) {
      p <- rep(p, 6)
    }
    if (!is.numeric(p)) {
      stop("p must be between numeric!")
    }
    x <- sample(1:6, n, p, replace = TRUE)
    return(x)
  }


    #' Urn
    #'
    #'@param n number of balls to draw from urn (integer value, default = 20)
    #'@param r number of RED balls in urn (default = 5)
    #'@param b number of BLUE balls in urn (default = 7)
    #'@param g number of GREEN balls in urn (default = 10)
    #'@param repl: with replacement (TRUE/FALSE)
    #'@return vector of 3 elements: numbers of balls drawn from each color
    #'@export
    polyUrn <- function(n = 20, r = 5, b = 7, g = 10, repl = FALSE) {
      if (!testInt(n)) {
        stop("n must be an integer!")
      }
      if (!testInt(r)) {
        stop("r must be an integer!")
      }
      if (!testInt(b)) {
        stop("b must be an integer!")
      }
      if (!testInt(g)) {
        stop("g must be an integer!")
      }
      if (n > (r + b + g)) {
        stop("Trying to draw more balls that urn content!")
      }
      urn <- c(rep("R", r), rep("B", b), rep("G", g))
      draw <- sample(urn, size = n, replace = repl)
      x <- c(sum(draw == "R"), sum(draw == "B"), sum(draw == "G"))
      return(x)
    }


      #' Tossing a coin
      #'
      #'@param n number of coins to toss (integer value, default = 1)
      #'@param p probability of a 'success' (Tail/Head), default = 0.5
      #'@return either 'Pile' or 'Face'
      #'@export
      CoinToss <- function(n = 1, p = 0.5) {
        if (!testInt(n)) {
          stop("n must be an integer!")
        }
        if (p < 0 | p > 1) {
          stop("p must be between 0 and 1!")
        }
        coin <- c("Pile", "Face")
        x <- sample(coin, n, c(p, 1 - p), replace = TRUE)
        return(x)
      }

      
      #'
      #' Verify that x is an integer value
      #'
      #'@param x integer value to check: Must be numeric
      #'@return TRUE if integer, FALSE if not
      #'@export
      testInt <- function(x) {
        if (!is.numeric(x)) {
          stop("Argument must be numeric")
        }
        t <- all.equal(x, as.integer(x))
        if (t == TRUE) {
          return(TRUE)
        } else {
          return(FALSE)
        }
      }
      
      
      # Calcul du nombre d'intervalles de classe pour n données quantitatives
      n_classes <- function(x){
        nc <- ceiling((length(x)^(1/3)*(max(x)-min(x)))/(2*IQR(x)))
        return(nc)
      }
      
      
      # Calcul du mode à partir de la courbe de densité
      DensMode <- function(x){
        d<-density(x)
        m = d$x[which.max(density(x)$y)]
        return(m)
      }
      
      # From: https://r-coder.com/exponential-distribution-r/
      # rate: lambda parameter
      # lb: lower bound of the area
      # ub: upper bound of the area
      # acolor: color of the area
      # ...: additional arguments to be passed to the lines function
      exp_area <- function(rate = 1, lb, ub, acolor = "lightgray", ...) {
        x <- seq(0, 12/rate, 0.01) 
        
        if (missing(lb)) {
          lb <- min(x)
        }
        if (missing(ub)) {
          ub <- max(x)
        }
        
        x2 <- seq(lb, ub, length = 100)    
        plot(x, dexp(x, rate = rate), type = "n", ylab = "Densité", xlim=c(0,qexp(0.999,rate)))
        
        y <- dexp(x2, rate = rate)
        polygon(c(lb, x2, ub), c(0, y, 0), col = acolor)
        lines(x, dexp(x, rate = rate), type = "l", ...)
      }