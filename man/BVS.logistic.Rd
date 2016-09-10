\name{BVS.logistic}
\alias{BVS.logistic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian Variable Selection Algorithm for Logistic Regression
}
\description{
Fit a Bayesian linear model via VB.
}
\usage{
BVS.logistic(x, y, m = NULL, v1 = 1, a0 = 1, b0 = 1, iter.max = 10^3, thres = 10^-3, normalize = TRUE, trace = FALSE, param.ini = list(alpha = NULL, mu = NULL, phi = NULL, sigma2 = NULL, theta = NULL, w.E = NULL), keep.data = FALSE, phi.c = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Design matrix, of dimension \eqn{n \times p}. Each row is an observation vector, and each column is a feature vector.
}
  \item{y}{
Response variable of length \eqn{n}. \code{y} could be either binary taking 0 or 1, or positive integers (count of successes).
}
  \item{m}{
If \code{y} measures the count of successes, then \code{m} is a vector of number of random trials.
}
  \item{v1}{
The parameter related to the slab component of the prior distribution of \eqn{\beta_j}'s. Default is 1. 
}
  \item{a0}{
The shape parameter \eqn{a_0} of the prior of \eqn{\theta}.
}
  \item{b0}{
The shape parameter \eqn{b_0} of the prior of \eqn{\theta}.
}
  \item{iter.max}{
The maximum number of iterations. 
}
  \item{thres}{
The convergence threshold.
}
  \item{normalize}{
Normalize the data or not. After normalization, each column of the design matrix will have mean 0, and variance equal to \eqn{n}. Default is \code{TRUE}.
}
  \item{trace}{
Record the intermediate output or not. 
}
  \item{param.ini}{
User-specified initial values in a list. 
}
  \item{keep.data}{
Keep the data in the final output or not. Default is \code{FALSE}.
}
  \item{phi.c}{
Threshold value that controls whether to update \code{phi} or not: if \code{phi[j]} is either greater than \code{1/(1+exp(-phi.c)} or less than \code{1/(1+exp(phi.c)}, \code{phi[j]} will not be updated in the following iterations.
}
}
\details{
Performs VB algorithm for Bayesian logistic regression with a spike-and-slab prior on \eqn{\boldsymbol{\beta}}.
}
\value{
The output is a list of parameters of variational distributions and point estimates, including: 
\item{mu}{Vector of means of the slab component of the posterior of \eqn{\beta_j}'s.}
\item{sigma_2}{Vector of variances of the slab component of the posterior of \eqn{\beta_j}'s.}
\item{phi}{Vector of posterior inclusion probabilities of \eqn{\beta_j}'s.}
\item{w.E}{Posterior mean of the Polya-Gamma latent variables.}
\item{alpha}{Point estimate of the intercept.}
\item{theta}{Point estimates of \eqn{\theta}.}
\item{post.save}{A list contains intermediate estimates of \code{mu}, \code{sigma2}, \code{phi}, \code{w.E}, \code{alpha}, and \code{theta}.}
\item{x.mean}{The column means of the design matrix.}
\item{x.sd}{The column sample standard deviation of the design matrix.}
\item{remove.var}{The removed variables' id, i.e., if a variable is a constant, it will be removed.}
\item{x}{The design matrix if \code{keep.data=TRUE}. If \code{Normalize=TRUE}, \code{x} is the matrix after normalization.}
\item{y}{The response variable subtracted by its mean if \code{keep.data=TRUE}.}
}

\references{
%% ~put references to the literature/web site here ~
}
\author{
Xichen Huang
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Generate data
n = 100; p = 20; 
v1 = 1
beta.0 = c(3,2,1)
beta.0 = c(beta.0, rep(0, p-length(beta.0)))
set.seed(100)
x = matrix(rnorm(n*p), n, p)
y = x\%*\%beta.0 + rnorm(n, sd=1)
y = ifelse(y>0, 1, 0)
# Model fit
bvs.fit = BVS.logistic(x, y, v1=1)
# cbind(beta.0, with(bvs.fit, cbind(mu, sigma2, phi)))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
