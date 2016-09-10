\name{predict.BVS.logistic}
\alias{predict.BVS.logistic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prediction for Bayesian Logistic Regression
}
\description{
To predict the response given new \eqn{X}. 
}
\usage{
predict.BVS.logistic(BVS, x.new, m.new = NULL, type = 1, phi.thres = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BVS}{
The output of \code{BVS.logistic}.
}
  \item{x.new}{
New \eqn{X} matrix.
}
  \item{m.new}{
New number of random trials \eqn{m}. 
}
  \item{type}{
\code{type=1}: use sparse estimated coefficients, \code{beta=mu*[phi>phi.thres]}. \code{type=2}: use dense estimated coefficients (posterior mean), \code{beta=mu*phi}.
}
  \item{phi.thres}{
Cut-off value of \code{phi}, so that only those with \code{phi>phi.thres} will be selected when \code{type=1}. 
}
}
\details{
Perform prediction from the model output of \code{BVS.logistic} and \code{BVS.logistic.cv}. 
}
\value{
A list of:
\item{y}{Class label or number of success.}
\item{p}{Probability of success.}
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
# Generate training data
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
# Generate test data
n.test = 100
x.test = matrix(rnorm(n.test*p), n.test, p)
y.test = x.test\%*\%beta.0 + rnorm(n.test, sd=1)
y.test = ifelse(y.test>0, 1, 0)
py.hat = predict.BVS.logistic(bvs.fit, x.test)
y.hat = py.hat$y
# table(y.test, y.hat)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
