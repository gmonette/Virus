\section{Crescent-shaped distribution}

To better display how HMC proposals follow the coutours of the target density, we need a bivariate distribution with unusual contours.  In this appendix we discuss the crescent-shaped distributions used in this paper.

The distributions are obtained with the composition of two transformations of a bivariate normal distribution.
The first transformation turns the normal into a crescent and the second rotates it.
Both transformations have unit Jacobian determinants so that the contours of
the transformed distribution are the transformation of the original elliptical normal contours.

We start with a normal variate $\left( {\begin{array}{*{20}{c}} x \\ y \end{array}} \right)$, centered at the origin with diagonal variance $\Sigma = $\left( {\begin{array}{*{20}{cc}} \sigma_x^2 & 0 \\ 0 & y \end{array}} \right) . 
the transformation:
  #' \[\left( {\begin{array}{*{20}{c}}
#' {x'} \\ 
#' {y'} 
#' \end{array}} \right) = \left( {\begin{array}{*{20}{c}}
#' x \\ 
#' {y + a{x^2}} 
#' \end{array}} \right)\]
#' and then rotating the result counterclockwise through an angle $\phi$.
#' 
#' Both of these transformations have unit Jacobian determinants so that the
#' density contours are transformed though the same transformation.
#'
#' The following functions implement these transformations.
#'
#' Note that $Rot^{-1}(.,\phi) = Rot(.,-\phi)$ and $Tr^{-1}(.,a) = Tr(.,-a)$
#'
#'  
#' ## Rotation in $R^2$
#' 
Rot <- function(x, phi) {
  # Counterclockwise rotation in R2 by phi radians 
  as.vector(cbind(c(cos(phi), sin(phi)), c(-sin(phi), cos(phi))) %*% as.vector(x))
}
Rot_jac <- function(x, phi)  cbind(c(cos(phi), sin(phi)), c(-sin(phi), cos(phi)))
Rot_inv <- function(y, phi) Rot(y, -phi)
Rot_inv_jac <- function(y, phi) Rot_jac(y, -phi)
#'
#' ## Quadratic transformation  
#' 
#' A transformation with a unit Jacobian determinant
#' that turns a horizontally stretched ellipse into a 'smile'
#' 
Tr <- function(x, a) c(x[1], x[2] + a * x[1]^2)
Tr_jac <- function(x, a) cbind(c(1, 2*a*x[1]), c(0,1))
Tr_inv <- function(y, a) Tr(y, -a) 
Tr_inv_jac <- function(y, a) Tr_jac(y, -a) 
#'
#' To use HMC we need the negative of the log density function
#' and its gradient. To visualize the trajectories,
#' We also wish to have a function to generate
#' contours of the density function.
#' 
#' Note that $a = 0$ produces a normal densities centered at 0.
#' 
#' 




f_gen <- function(Sigma) {
  # negative log-density of normal distribution up to additive constant
  function(x) {
    crossprod(x, solve(Sigma, x))/2
  }
}
f_grad_gen <- function(Sigma) {
  # gradient for f_gen
  function(x) {
    solve(Sigma, x)
  }
}
f_contour_gen <- function(Sigma) {
  # contour for f_gen
  function(deviance = 1, p, n = 100) {
    if(!missing(p)) deviance <- qchisq(p, df = nrow(Sigma))
    spida2::ell(c(0,0), radius = sqrt(deviance), shape = Sigma, n = n)
  }
}


#' ## Parametric family of banana distributions
#' 

f_gen <- function(phi, a, Sigma) function(x) {
  # negative of log of density
  # using the fact that Tr and Rot have unit jacobien determinants
  - dmvnorm(Tr_inv(Rot_inv(x, phi), a), sigma = Sigma, log = TRUE)
}
f_grad_gen <- function(phi, a, Sigma) function(x) {
  # gradient of f
  x1 <- Rot_inv(x, phi)
  x2 <- Tr_inv(x1, a)
  c(t(Rot_inv_jac(x, phi)) %*% t(Tr_inv_jac(x1, a)) %*% solve(Sigma, x2))
}
f_contour_gen <- function(phi, a, Sigma) function(deviance = 1, n = 100) {
  # With unit jacobian determinants, the contour of the density of
  # the transformation is the transformation of the contour of the density
  e <- spida2::ell(c(0,0), radius = sqrt(deviance), shape = Sigma, n = n)
  e2 <- apply(e, 1, Tr, a)
  t(apply(e2, 2, Rot, phi))
}

#'
#' ### Generating f, fg and fc
#'

phi <- pi/3 
a <- 0.03
Sigma <- cbind( c(36,0),c(0,1))
f <- f_gen(phi, a, Sigma)
fg <- f_grad_gen(phi, a, Sigma)
fc <- f_contour_gen(phi, a, Sigma)



hamTraj <- function(x0, m0, neglogf, grad, stepsize=1, steps=20, verbose = FALSE){
  # Hamiltonian trajectory with leapfrog
  xmat <- matrix(NA,steps+1,length(x0))
  mmat <- matrix(NA,steps+2,length(x0))
  xmat[1, ] <- x <- x0
  mmat[1, ] <- mom <- m0
  mmat[2, ] <- mom <- mom - 0.5*stepsize*grad(x0)  # half-step for momentum
  for (step in 1:(steps - 1)){
    xmat[step+1, ] <- x <- x + stepsize*mom
    mmat[step+2, ] <- mom <- mom - stepsize*grad(x)
  }
  xmat[steps + 1, ] <- x <- x + stepsize*mom # last step for position
  mmat[steps + 2, ] <- mom <- mom  - 0.5 * stepsize * grad(x)  # last half-step for momentum
  a <- exp(neglogf(x0) - neglogf(x) + 0.5*sum(m0^2) - 0.5*sum(mom^2))
  p <- runif(1)
  move <- p < a
  list(x = xmat, m = mmat, a = a , p = p, move = move)
}

#'
#' Trajectory using a stepsize of 0.01
#'
#+ fig.keep='all'
# initialize

par(bg='white')
plot(fc(c(.001,seq(1,6)^2)),type = 'l', lty = 2, ylim = c(-15,50), xlim=c(-50,15), asp = 1, xlab = '',ylab = '')

stepsize <- 0.01
steps <- 200
pch <- 16
cex = 1
col  <- 'red'

x <- c(0,30)
m <- c(0,0)

points( rbind(x), pch = 16, col = 'red')
ac <- 0
for(i in 1:10) {
  ret <- hamTraj(x,m, f, fg, stepsize = stepsize, steps = steps)
  lines(ret$x,col=col)
  ac <- ac + (ret$a-1)^2
  points(tail(ret$x,1), pch = pch, col = col)
  x <- c(tail(ret$x,1))
  m <- c(tail(ret$m,1))
}
for(i in 11:50) {
  ret <- hamTraj(x,m, f, fg, stepsize = stepsize, steps = steps)
  lines(ret$x,col=col)
  ac <- ac + (ret$a-1)^2
  points(tail(ret$x,1), pch = pch, col = col)
  x <- c(tail(ret$x,1))
  m <- c(tail(ret$m,1))
}
for(i in 51:100) {
  ret <- hamTraj(x,m, f, fg, stepsize = stepsize, steps = steps)
  lines(ret$x,col=col)
  ac <- ac + (ret$a-1)^2
  points(tail(ret$x,1), pch = pch, col = col)
  x <- c(tail(ret$x,1))
  m <- c(tail(ret$m,1))
}
cat("Cumulative squared energy disparity:", ac,'\n' )

#' 
#' Note that the particle seems to escape at the end, or, at least, achieves a higher
#' energy than its initial energy.
#' 
#' Repeat with smaller stepsize
#'
#+ fig.keep='all',echo=FALSE
{
  plot(fc(c(.001,seq(1,6)^2)),type = 'l', lty = 2, ylim = c(-15,50), xlim=c(-50,15),asp=1,bg='white', xlab = '',ylab = '')
  
  stepsize <- .001
  steps <- 2000
  pch <- 16
  cex = 1
  col  <- 'blue'
  
  x <- c(0,30)
  m <- c(0,0)
  
  ac <- 0
  for(imax in c(5,10,40,100)){
    x <- c(0,30)
    m <- c(0,0)
    plot(fc(c(.001,seq(1,6)^2)),type = 'l', lty = 2, ylim = c(-15,50), xlim=c(-50,15),asp = 1, xlab = '',ylab = '')
    for(i in 1:imax) {
      ret <- hamTraj(x,m, f, fg, stepsize = stepsize, steps = steps)
      lines(ret$x,col=col)
      ac <- ac + (ret$a-1)^2
      points(tail(ret$x,1), pch = pch, col = col)
      x <- c(tail(ret$x,1))
      m <- c(tail(ret$m,1))
    }
  }
  cat("Cumulative squared energy disparity:", ac,'\n' )
}
#' 
#' There was no escape this time.
#' 
#' Repeating with larger stepsize
#'
#+ fig.keep='all'
{
  plot(fc(c(.001,seq(1,6)^2)),type = 'l', lty = 2, ylim = c(-15,50), xlim=c(-50,15), asp = 1, xlab = '',ylab = '')
  stepsize <- .1
  steps <- 20
  pch <- 21
  cex = 1.3
  col  <- 'black'
  
  for(imax in c(5,10,40,100)){
    x <- c(0,30)
    m <- c(0,0)
    ac <- 0
    plot(fc(c(.001,seq(1,6)^2)),type = 'l', lty = 2, ylim = c(-15,50), xlim=c(-50,15), asp = 1, xlab = '',ylab = '')
    for(i in 1:imax) {
      ret <- hamTraj(x,m, f, fg, stepsize = stepsize, steps = steps)
      lines(ret$x,col=col)
      ac <- ac + (ret$a-1)^2
      points(tail(ret$x,1), pch = pch, col = col)
      x <- c(tail(ret$x,1))
      m <- c(tail(ret$m,1))
    }
  }
  cat("Cumulative squared energy disparity:", ac,'\n' )
}
#'
#'
#' ## Exploring different trajectories
#'
#' ### Dropping
#' 
#+ fig.keep='all' 
{
  stepsize <- .001
  steps <- 2000
  pch <- 16
  cex = 1
  col  <- 'red'
  
  for( imax in c(5,10,40,100)){
    plot(fc(c(.001,seq(1,6)^2)),type = 'l', lty = 2, ylim = c(-15,50), xlim=c(-50,15),asp=1,, xlab = '',ylab = '')
    x <- c(0,29)
    m <- c(0,0)
    for(i in 1:imax) {
      points(rbind(x), col = col, pch = pch)
      ret <- hamTraj(x,m, f, fg, stepsize = stepsize, steps = steps)
      lines(ret$x,col=col)
      points(tail(ret$x,1), pch = pch, col = col)
      x <- c(tail(ret$x,1))
      m <- c(tail(ret$m,1))
    }
  }
}
#'
#' Stepsizes of 0.01 and 0.001 yielded the same trajectory (visually) but
#' using 0.1 produced divergence towards the end of the trajectory.
#' 
#' ### Banana-shaped orbit
#' 
#+ fig.keep='all' 
{
  stepsize <- .001
  steps <- 2000
  pch <- 16
  cex = 1
  col  <- 'red'
  
  for( imax in c(5,10,40,100)){
    plot(fc(c(.001,seq(1,6)^2)),type = 'l', lty = 2, ylim = c(-15,50), xlim=c(-50,15),asp=1, xlab = '',ylab = '')
    x <- c(0,29)
    m <- 1.5* c(1,-1)*rev(fg(x))
    for(i in 1:imax) {
      points(rbind(x), col = col, pch = pch)
      ret <- hamTraj(x,m, f, fg, stepsize = stepsize, steps = steps)
      lines(ret$x,col=col)
      points(tail(ret$x,1), pch = pch, col = col)
      x <- c(tail(ret$x,1))
      m <- c(tail(ret$m,1))
    }
  }
}
#'
#'
#'
#' ## Trajectories in elliptical bowl
#'

Sigma <- cbind( c(36,0), c(0,1))
f <- f_gen(0, 0, Sigma)
fg <- f_grad_gen(0, 0, Sigma)
fc <- f_contour_gen(0, 0, Sigma) 

#'
#' ### Dropping
#'
#+ fig.height=5,fig.width=8,fig.keep='all'
# windows( 16, 10)
{
  stepsize <- .001
  steps <- 2000
  pch <- 16
  cex = 1
  col  <- 'red'
  
  for( imax in c(5,10,40,100)){
    plot(fc(c(.001,seq(1,6)^2)),type = 'l', lty = 2, ylim = c(-10,10),asp=1, xlab = '',ylab = '')
    x <- c(-10,8)
    m <- c(0,0)
    for(i in 1:imax) {
      points(rbind(x), col = col, pch = pch)
      ret <- hamTraj(x,m, f, fg, stepsize = stepsize, steps = steps)
      lines(ret$x,col=col)
      points(tail(ret$x,1), pch = pch, col = col)
      x <- c(tail(ret$x,1))
      m <- c(tail(ret$m,1))
    }
  }
}

#'
#' ### Orbit?
#'
#+ fig.height=5,fig.width=8,fig.keep='all'
# windows( 16, 10)
{
  stepsize <- .001
  steps <- 2000
  pch <- 16
  cex = 1
  col  <- 'red'
  
  for( imax in c(5,10,40,100)){
    plot(fc(c(.001,seq(1,6)^2)),type = 'l', lty = 2, ylim = c(-10,10),asp=1, xlab = '',ylab = '')
    x <- c(10,8)
    m <- .5* c(1,-1)*rev(fg(x))
    for(i in 1:imax) {
      points(rbind(x), col = col, pch = pch)
      ret <- hamTraj(x,m, f, fg, stepsize = stepsize, steps = steps)
      lines(ret$x,col=col)
      points(tail(ret$x,1), pch = pch, col = col)
      x <- c(tail(ret$x,1))
      m <- c(tail(ret$m,1))
    }
  }
}

#'
#' ## Trajectories in spherical bowl
#'

Sigma <- cbind( c(1,0), c(0,1))
f <- f_gen(0, 0, Sigma)
fg <- f_grad_gen(0, 0, Sigma)
fc <- f_contour_gen(0, 0, Sigma) 

#'
#' ### Dropping
#'
#+ fig.height=7,fig.width=8,fig.keep='all'
# windows( 16, 10)
{
  stepsize <- .001
  steps <- 2000
  pch <- 16
  cex = 1
  col  <- 'red'
  
  for( imax in c(5,10,40,100)){
    plot(fc(c(.001,seq(1,6)^2)),type = 'l', lty = 2, ylim = c(-10,10),asp=1, xlab = '',ylab = '')
    x <- c(-6,5)
    m <- c(0,0)
    for(i in 1:imax) {
      points(rbind(x), col = col, pch = pch)
      ret <- hamTraj(x,m, f, fg, stepsize = stepsize, steps = steps)
      lines(ret$x,col=col)
      points(tail(ret$x,1), pch = pch, col = col)
      x <- c(tail(ret$x,1))
      m <- c(tail(ret$m,1))
    }
  }
}

#'
#' ### Orbit?
#'
#+ fig.height=5,fig.width=8,fig.keep='all'
# windows( 16, 10)
{
  stepsize <- .001
  steps <- 2000
  pch <- 16
  cex = 1
  col  <- 'red'
  
  for( imax in c(5,10,40,100)){
    plot(fc(c(.001,seq(1,6)^2)),type = 'l', lty = 2, ylim = c(-10,10),asp=1, xlab = '',ylab = '')
    x <- c(-6,5)
    m <- 0.7* c(1,-1)*rev(fg(x))
    for(i in 1:imax) {
      points(rbind(x), col = col, pch = pch)
      ret <- hamTraj(x,m, f, fg, stepsize = stepsize, steps = steps)
      lines(ret$x,col=col)
      points(tail(ret$x,1), pch = pch, col = col)
      x <- c(tail(ret$x,1))
      m <- c(tail(ret$m,1))
    }
  }
}
