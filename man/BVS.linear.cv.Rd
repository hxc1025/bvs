\name{BVS.linear.cv}
\alias{BVS.linear.cv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cross-Validation for BVS.linear.1 and BVS.linear.2
}
\description{
Given a sequence of \eqn{v_1}, perform a cross-validation on the training data, and select the best \eqn{v_1} that returns the minimal cross-validation error, in terms of the \eqn{L-2} loss. 
}
\usage{
BVS.linear.cv(x, y, v1.all = NULL, v1.length = 100, FUN = BVS.linear.1, fold = 5, L = 100, a0 = 2, b0 = 2, lambda = 1, nu = 1, iter.max = 100, thres = 10^-3, normalize = TRUE, trace = FALSE, param.ini = list(mu = NULL, phi = NULL, sigma_j_sq = NULL, sigma2 = NULL, theta = NULL), my.seed = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Design matrix, of dimension \eqn{n \times p}. Each row is an observation vector, and each column is a feature vector.
}
  \item{y}{
Response variable of length \eqn{n}.
}
  \item{v1.all}{
A sequence of \eqn{v_1}.
}
  \item{v1.length}{
Not available now.
}
  \item{FUN}{
Which algorithm to use: Algorithm 1 (component-wise) vs Algorithm 2 (Batch-wise. Use \code{BVS.linear.1} or \code{BVS.linear.2}.)
}
  \item{fold}{
Number of folds of data in cross-validation. 
}
  \item{L}{
Not available now.
}
  \item{a0}{
The shape parameter \eqn{a_0} of the prior of \eqn{\theta}.
}
  \item{b0}{
The shape parameter \eqn{b_0} of the prior of \eqn{\theta}.
}
  \item{lambda}{
The parameter \eqn{\lambda} of the prior of \eqn{\sigma^2}.
}
  \item{nu}{
The parameter \eqn{\nu} of the prior of \eqn{\sigma^2}.
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
  \item{my.seed}{
Random seed for sampling cross-validation data.}
}
\details{
This is a cross-validation function to choose \eqn{v_1} for \code{BVS.linear.1} and \code{BVS.linear.2}.  
}
\value{
The output is a list that contains: 
  \item{v1}{All the \eqn{v_1} used in CV.}
  \item{cvm}{Average CV error.}
  \item{cvsd}{Standard deviation of CV error.}
  \item{cvup}{Average CV error plus one s.d..}
  \item{cvlo}{Average CV error minus one s.d..}
  \item{v1.min}{The \eqn{v_1} that returns the minimum CV error.}
  \item{v1.1se}{The largest \eqn{v_1} that has CV error within one s.d. of the minimum CV error.}
  \item{x.mean}{The column means of the design matrix.}
  \item{x.sd}{The column sample standard deviation of the design matrix.}
  \item{remove.var}{The removed variables' id, i.e., if a variable is a constant, it will be removed.}
  \item{y.mean}{Sample mean of y.}
}

\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{BVS.linear.1}, \code{BVS.linear.2}
}
\examples{
# Generate data
n = 40; sigma = 3
beta.0 = c(3, 1.5, 0, 0, 2, 0, 0, 0)
p = length(beta.0)
rho = 0.5
Sigma = rho^abs(outer(1:p, 1:p,"-"))
set.seed(100)
x = rMVN(n, rep(0,p), Sigma)
y = x\%*\%beta.0 + rnorm(n, 0, sigma)
# CV
set.seed(200)
bvs.cv = BVS.linear.cv(x, y, v1.all = 10^seq(-2,2,0.1), FUN=BVS.linear.2)
bvs.fit = BVS.linear.2(x, y, v1=bvs.cv$v1.min)
# result
# with(bvs.fit, cbind(mu, sigma_j_sq, phi))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
