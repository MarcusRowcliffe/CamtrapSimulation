  library(circular)

## rautonorm
##Generates a set of autocorrelated random normal variates
#
#INPUT
# n: number of variates to generate
# mean, sd: mean and standard deviation of the normal distribution
# r: the autocorrelation coefficient (between 0 and 1)
rautonorm <- function(n,mean=0,sd=1,r){
  ranfunc <- function(i,z,r) sqrt(1-r^2) * sum(z[2:(i+1)]*r^(i-(1:i))) + z[1]*r^i
  z <- rnorm(n)
  mean + sd*c(z[1], sapply(1:(n-1), ranfunc, z, r))
}

## pathgen
##Generates a path of x,y positions using a correlated random walk
#
#INPUT
# n: number of steps
# pTurn: probability of turning at each step
# kTurn: mean vonMises concentration parameter (kappa) for turn angle (higher=more concentrated)
# logspeed: mean log speed
# speedSD: standard deviation of log speed
# speedCor: autocorrelation in speeed
# kCor: whether to correlate kappa with speed
# xlim, ylim: x and y axis limits within which to pick the starting point
# wrap: whether to wrap the path
#
#OUTPUT
#A list with elements:
# path: a dataframe with columns x and y (path co-ordinates) and, if wrap=TRUE, breaks indicating where wrap breaks occur
# turn, absturn: radian (absolute) turn angles for each step (turn ranging 0 to 2pi; absturn ranging 0 to pi)
# speed: step speeds
pathgen <- function(n, kTurn=0, logspeed=0, speedSD=0, speedCor=0, kCor=TRUE, pTurn=1, xlim=c(0,0), ylim=xlim, wrap=FALSE){
  spds <- exp(rautonorm(n, logspeed, speedSD, speedCor))
  tTurn <- rbinom(n,1,pTurn)
  if(kCor==TRUE){
    kappas <- kTurn * spds / mean(spds)
    deviates <- sapply(kappas, function(x) as.numeric(rvonmises(1,circular(0),x)))
  } else
    deviates <- as.numeric(rvonmises(n, circular(0), kTurn))
  deviates[tTurn==0] <- 0
  angles <- runif(1)*2*pi + cumsum(deviates)
  x <- c(0, cumsum(spds*sin(angles))) + runif(1,xlim[1],xlim[2])
  y <- c(0, cumsum(spds*cos(angles))) + runif(1,ylim[1],ylim[2])
  absdevs <- deviates
  i <- absdevs>pi
  absdevs[i] <- 2*pi-absdevs[i]
  absdevs <- abs(absdevs)
  res <- list(path=data.frame(x,y), turn=deviates, absturn=absdevs, speed=spds)
  if(wrap) res <- wrap(res, xlim, ylim)
  res
}

# wrap
#Takes a path object created with pathgen and wraps the co-ordinates within given limits
#
#INPUT
# pth: a two column array of x,y positions defining the path
# xlim, ylim: the x,y limits within which to wrap the path
#
#OUTPUT
# A path object with x,y co-ordinates wrapped and breaks column added indicating wrap breaks
wrap <- function(path, xlim, ylim=xlim){
  pth <- path$path
  n <- nrow(pth)
  brkpnts <- vector()
  repeat{
    xout <- which(pth$x<xlim[1] | pth$x>xlim[2])[1]
    yout <- which(pth$y<ylim[1] | pth$y>ylim[2])[1]
    if(is.na(xout) & is.na(yout)) break else {
      if(!is.na(xout)){
        brkpnts <- c(brkpnts, xout)
        addn <- if(pth$x[xout]<xlim[1]) diff(xlim) else -diff(xlim) 
        pth$x[xout:n] <- pth$x[xout:n]+addn
      }
      if(!is.na(yout)){
        brkpnts <- c(brkpnts,yout)
        addn <- if(pth$y[yout]<ylim[1]) diff(ylim) else -diff(ylim) 
        pth$y[yout:n] <- pth$y[yout:n]+addn
      }
    }
  }
  brkn <- diff(c(1, sort(brkpnts), n+1))
  breaks <- rep(1:length(brkn), brkn)
  path$path <- cbind(pth, breaks)
  path
}

#plot_wrap
#Plots a wrapped path
#
#INPUT
# path: a wrapped path object created by pathgen
# type: l(ine), p(oint) or b(oth)
# add: add to existing plot or create new one
# axisargs, lineargs, pointargs: lists of arguments to control axis lines or point characteristics
plot_wrap <- function(path, type=c("l","p","b"), add=FALSE, axisargs=list(), lineargs=list(), pointargs=list()){
  type <- match.arg(type)
  if(!"xlab" %in% names(axisargs)) axisargs <- c(xlab="", axisargs)
  if(!"ylab" %in% names(axisargs)) axisargs <- c(ylab="", axisargs)
  if(!add) do.call("plot", c(list(path$path[,1:2], type="n"), axisargs, asp=1))
  for(i in unique(path$path$breaks)) {
    j <- path$path$breaks==i
    xy <- subset(path$path, j)
    pargs <- lapply(pointargs, function(x) if(length(x)==nrow(path$path)) x[j] else x)
    if(type %in% c("l", "b")) do.call("lines", c(list(xy$x, xy$y), lineargs))
    if(type %in% c("p", "b")) do.call("points", c(list(xy$x, xy$y), pargs))
  }
}

#is_in_dz
#Defines whether points are with detection zones
#
#INPUT
# point: a two colummn x,y array of point positions
# dzone: four column array of parameters defining a sector-shaped detection zone
#        required column headings:
#           x,y: x,y coordinates of camera
#           r, th: detection zone radius and angle
#           dir: radian direction in which the camera is facing
#
#OUTPUT
#A logical array defining whether each point (rows) is in each detection zone (columns)
is_in_dz <- function(point, dzone){
  ij <- expand.grid(1:nrow(point), 1:nrow(dzone)) #expanding rows for each point and dzone
  pt <- point[ij$Var1, ]
  dz <- dzone[ij$Var2, ]
  dist <- sqrt((pt[, 1]-dz$x)^2 + (pt[, 2]-dz$y)^2) #distance from camera to point
  bear <- atan((pt[, 1]-dz$x) / (pt[, 2]-dz$y)) + #bearing from camera point
    ifelse(pt[, 2]<dz$y, pi, ifelse(pt[, 1]< dz$x, 2*pi, 0))
  beardif <- (bear-dz$dir) %% (2*pi) #abs angle between bear and dzone centre line
  beardif <- ifelse(beardif>pi, 2*pi-beardif, beardif)
  res <- ifelse(beardif < dz$th/2 & dist<dz$r, TRUE, FALSE)
  matrix(res, nrow=nrow(point))
}

# plot_dzone
#Convenience function for plotting detection zones (adds to existing plot)
#INPUT
# dzone: a four column array of parameters as defined above
plot_dzone <- function(dzone, ...){
  for(i in 1:nrow(dzone)){
    sq <- with(dzone[i, ], seq(dir-th/2, dir+th/2, len=50))
    poly <- with(dzone[i, ], cbind(x + c(0, r*sin(sq)), y + c(0, r*cos(sq))))
    polygon(poly, ...)
  }
}

## sequence_data
#1. Takes dataframes defining a path and a detection zone (as defined above)
#2. Filters the path points falling within the detection zone
#3. Assigns each contiguous sequence of points a unique sequence identifier
#4. calculates the distances between points with sequences
#
#INPUT
# pth: a path object
# a detection zone array
#OUTPUT#
# A data frame with columns:
# x,y: x,y co-ordinates of sequence points in detection zones
# sequenceID: integer sequence identifier
# distance: distance traveled for each step between points
sequence_data <- function(pth, dzone){
  pth <- pth$path[, c("x","y")]
  isin <- is_in_dz(pth, dzone)
  isin[1,] <- FALSE
  isin <- as.vector(isin)
  pth <- pth[rep(1:nrow(pth), nrow(dzone)), ]
  newseq <- tail(isin, -1) > head(isin, -1)
  seqid <- c(0, cumsum(newseq))[isin]
  xy <- pth[isin, ]
  dist <- sqrt(diff(xy$x)^2 + diff(xy$y)^2)
  newseq <- tail(seqid, -1) > head(seqid, -1)
  dist[newseq] <- NA
  dist <- c(NA, dist)
  data.frame(xy, sequenceID=seqid, distance=dist)
}

## calc_speed
#Summarises speeds for a dataframe of position sequences
#
#INPUT
# dat: a dataframe of position observations created by sequence_data (above)
#
#OUTPUT
#A dataframe with a row per sequence and columns:
# sequenceID: integer sequence identifier
# distance: total distance travelled during the sequence
# points: number of points in the sequence
# speed: overall sequence speed
calc_speed <- function(dat){
  dist <- with(dat, tapply(distance, sequenceID, sum, na.rm=TRUE))
  points <- with(dat, tapply(distance, sequenceID, length))
  speed <- dist/(points-1)
  data.frame(sequenceID=unique(dat$sequenceID), distance=dist, points=points, speed=speed)
}

