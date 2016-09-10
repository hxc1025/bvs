\name{BVS.linear.2}
\alias{BVS.linear.2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian Variable Selection Algorithm for Linear Regression (Algorithm 2: Batch-wise)
}
\description{
Fit a Bayesian linear model via VB. The updating scheme of \eqn{\mu}, \eqn{\sigma_j^2}, and \eqn{\phi} is batch-wise.
}
\usage{
BVS.linear.2(x, y, v1 = 1, a0 = 2, b0 = 2, lambda = 1, nu = 1, iter.max = 100, thres = 10^-3, normalize = TRUE, trace = FALSE, param.ini = list(mu = NULL, phi = NULL, sigma_j_sq = NULL, sigma2 = NULL, theta = NULL), use.lambda_n1=FALSE, keep.data=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Design matrix, of dimension \eqn{n \times p}. Each row is an observation vector, and each column is a feature vector.
}
  \item{y}{
Response variable of length \eqn{n}.
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
  \item{use.lambda_n1}{
Whether use the minimum none-zero eigenvalue of \eqn{\mathbf{X}^T \mathbf{X}} or not. Default is \code{FALSE}. 
}
  \item{keep.data}{
Keep the data in the final output or not. Default is \code{FALSE}.
}
}
\details{
The linear regression model is
\deqn{
  y_i = \x_i^T \boldsymbol{\beta} + \epsilon_i, \, i=1,\ldots, n,
}
where \eqn{\epsilon \sim \N(0, \sigma^{2})}, and \eqn{\boldsymbol{\beta} = (\beta_{1}, \ldots, \beta_{p})^{T} } is the unknown coefficient vector of length \eqn{p}.

The hierarchical prior is specified as follows:
\deqn{
 \beta_{j} \vert \gamma_{j} &\sim &\gamma_{j}\N(0, v_{1}\sigma^2) + (1-\gamma_{j})\delta_{0}, \\
  \gamma_{j} &\overset{i.i.d.}{\sim} & \operatorname{Bern}(\theta), \nonumber  \\
  \sigma^2 &\sim& \IG\left(\frac{\nu}{2}, \frac{\nu\lambda}{2}\right),  \nonumber \\
  \theta &\sim &\operatorname{Beta}(a_{0}, b_{0}),  \nonumber
}
where \eqn{j=1, \ldots, p}, and \eqn{\nu}, \eqn{\lambda}, \eqn{a_0} and \eqn{b_0} are hyper-parameters.

}
\value{
The output is a list of parameters of variational distributions and point estimates, including: 
\item{mu}{Vector of means of the slab component of the posterior of \eqn{\beta_j}'s.}
\item{sigma_j_sq}{Vector of variances of the slab component of the posterior of \eqn{\beta_j}'s.}
\item{phi}{Vector of posterior inclusion probabilities of \eqn{\beta_j}'s.}
\item{sigma2}{Point estimates of \eqn{\sigma^2}.}
\item{theta}{Point estimates of \eqn{\theta}.}
\item{post.save}{A list contains intermediate estimates of mu, sigma_j_sq, phi, sigma2, and theta.}
\item{x.mean}{The column means of the design matrix.}
\item{x.sd}{The column sample standard deviation of the design matrix.}
\item{remove.var}{The removed variables' id, i.e., if a variable is a constant, it will be removed.}
\item{y.mean}{Sample mean of y.}
\item{x}{The design matrix if \code{keep.data=TRUE}. If \code{Normalize=TRUE}, \code{x} is the matrix after normalization.}
\item{y}{The response variable subtracted by its mean if \code{keep.data=TRUE}.}
}

\references{
%% ~put references to the literature/web site here ~
}
\author{
Xichen Huang}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{BVS.linear.1}, \code{BVS.linear.cv}, \code{predict.BVS.linear}
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
# Model Fit
bvs.fit = BVS.linear.2(x, y, v1=1)
# result
# with(bvs.fit, cbind(mu, sigma_j_sq, phi))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line