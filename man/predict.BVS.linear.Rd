\name{predict.BVS.linear}
\alias{predict.BVS.linear}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prediction for Bayesian Linear Regression
}
\description{
To predict the response given new \eqn{X}. 
}
\usage{
predict.BVS.linear(BVS, x.new, type = 1, phi.thres = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BVS}{
The output of \code{BVS.linear.1} or \code{BVS.linear.2}.
}
  \item{x.new}{
New \eqn{X} matrix.
}
  \item{type}{
Takes value of 1, 2, OR 3. \code{type=1}: use sparse estimated coefficients, \code{beta=mu*[phi>phi.thres]}. \code{type=2}: use dense estimated coefficients (posterior mean), \code{beta=mu*phi}. \code{type=3}: refit a logistic regression using \code{glm} and predict based on the refitted model. 
}
  \item{phi.thres}{
Cut-off value of \code{phi}, so that only those with \code{phi>phi.thres} will be selected when \code{type=1} or \code{type=3}. 
}
}
\details{
Perform prediction from the model output of \code{BVS.linear.1} and \code{BVS.linear.2}. 
}
\value{
A vector of predicted response.
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
# training data
n = 100
p = 10
b0 = c(1,2,3, rep(0, 7))
x = matrix(rnorm(n*p), n, p)
y = 2 + x\%*\%b0 + rnorm(n)
# fit model
m2 = BVS.linear.2(x, y)
# test data
x.test = matrix(rnorm(10*p), 10, p)
y.test = 2 + x.test\%*\%b0 + rnorm(10)
# prediction
y.hat = predict.BVS.linear(m2, x.test, type = 1)
cbind(y.test, y.hat)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
