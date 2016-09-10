\name{BVS.logistic.cv}
\alias{BVS.logistic.cv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cross-Validation for BVS.logistic
}
\description{
Given a sequence of \eqn{v_1}, perform a cross-validation on the training data, and select the best \eqn{v_1} that returns the minimal cross-validation error. 
}
\usage{
BVS.logistic.cv(x, y, m = NULL, v1.list = NULL, fold = 10, a0 = 1, b0 = 1, iter.max = 10^3, thres = 10^-3, normalize = TRUE, param.ini = list(alpha = NULL, mu = NULL, phi = NULL, sigma2 = NULL, theta = NULL, w.E = NULL), my.seed = NULL, LogLoss = FALSE, type = 1)
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
  \item{v1.list}{
A sequence of \eqn{v_1}.
}
  \item{fold}{
Number of folds of data in cross-validation. 
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
  \item{param.ini}{
User-specified initial values in a list. 
}
  \item{my.seed}{
Random seed for sampling cross-validation data.
}
  \item{LogLoss}{
Use log-loss or classification-error. 
}
  \item{type}{
Type of prediction. \code{type=1} uses \code{mu*(phi>0.5)}, while \code{type=2} uses \code{mu*phi}. 
}
}
\details{
This is a cross-validation function to choose \eqn{v_1} for \code{BVS.logistic}. 
}
\value{
The output is a list that contains: 
  \item{v1.best}{The \code{v1} that produces the lowest CV error.}
  \item{cv.err}{A matrix containing the CV error.}
  \item{BVS}{A fitted model output of BVS.logistic. Can be directly used for prediction using \code{BVS.logistic}.}
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
# Cross-validation
bvs.cv = BVS.logistic.cv(x, y, v1.list=10^seq(-2,2,0.1), my.seed=200)
% cbind(beta.0, with(bvs.cv$BVS, cbind(mu, sigma2, phi)))
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
