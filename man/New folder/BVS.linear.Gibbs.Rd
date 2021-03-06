\name{BVS.linear.Gibbs}
\alias{BVS.linear.Gibbs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BVS.linear.Gibbs(x, y, v1 = 1, a0 = 2, b0 = 2, lambda = 1, nu = 1, iter.max = 10^5, thres = 10^-3, normalize = FALSE, burn.in = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{v1}{
%%     ~~Describe \code{v1} here~~
}
  \item{a0}{
%%     ~~Describe \code{a0} here~~
}
  \item{b0}{
%%     ~~Describe \code{b0} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{nu}{
%%     ~~Describe \code{nu} here~~
}
  \item{iter.max}{
%%     ~~Describe \code{iter.max} here~~
}
  \item{thres}{
%%     ~~Describe \code{thres} here~~
}
  \item{normalize}{
%%     ~~Describe \code{normalize} here~~
}
  \item{burn.in}{
%%     ~~Describe \code{burn.in} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, y, v1 = 1, a0 = 2, b0 = 2, lambda = 1, nu = 1, iter.max = 10^5, 
    thres = 10^-3, normalize = FALSE, burn.in = 1000) 
{
    x = as.matrix(x)
    y = as.vector(y)
    if (normalize) {
        x.norm = Norm(x)
        n = x.norm$n
        p = x.norm$p
        x = x.norm$x
        x.mean = x.norm$x.mean
        x.sd = x.sd
        remove.var = x.norm$remove.var
    }
    else {
        n = nrow(x)
        p = ncol(x)
        x.mean = NULL
        x.sd = NULL
        remove.var = NULL
    }
    y.mean = mean(y)
    y = y - y.mean
    sigma2 = 1/rgamma(1, nu/2, (nu * lambda)/2)
    theta = rbeta(1, a0, b0)
    gamma = rbinom(p, 1, theta)
    beta = rnorm(p, 0, sqrt(v1 * sigma2)) * gamma
    res = Inf
    Ent = ent(gamma)
    kappa = drop(y - x \%*\% beta)
    xx = colSums(x^2)
    mu = rep(0, p)
    sigma_j_sq = rep(0, p)
    phi = rep(0, p)
    post.save = list(beta = matrix(NA, p, iter.max), gamma = matrix(NA, 
        p, iter.max), sigma2 = rep(NA, iter.max), theta = rep(NA, 
        iter.max))
    t = 0
    while (t < iter.max && res > thres) {
        t = t + 1
        for (j in 1:p) {
            kappa = kappa + x[, j] * beta[j]
            mu[j] = (x[, j] \%*\% kappa)/(xx[j] + 1/v1)
            sigma_j_sq[j] = sigma2/(xx[j] + 1/v1)
            BF = log(theta/(1 - theta)/sqrt(1 + v1 * xx[j])) + 
                mu[j]^2/sigma_j_sq[j]/2
            phi[j] = 1/(1 + exp(-BF))
            gamma[j] = rbinom(1, 1, phi[j])
            if (gamma[j] == 1) {
                beta[j] = rnorm(1, mu[j], sqrt(sigma_j_sq[j]))
            }
            else {
                beta[j] = 0
            }
            kappa = kappa - x[, j] * beta[j]
        }
        sigma2 = 1/rgamma(1, (n + nu + sum(gamma))/2, (sum(kappa^2) + 
            sum(beta[which(gamma != 0)]^2)/v1 + nu * lambda)/2)
        theta = rbeta(1, sum(gamma) + a0, p - sum(gamma) + b0)
        post.save$beta[, t] = beta
        post.save$gamma[, t] = gamma
        post.save$sigma2[t] = sigma2
        post.save$theta[t] = theta
        if ((t\%\%1000 == 0) && (t > burn.in)) {
            Ent.old = Ent
            Ent = ent(rowMeans(post.save$gamma[, (t - burn.in):t]))
            res = abs(Ent - Ent.old)
        }
    }
    post.save$beta = post.save$beta[, 1:t]
    post.save$gamma = post.save$gamma[, 1:t]
    post.save$sigma2 = post.save$sigma2[1:t]
    post.save$theta = post.save$theta[1:t]
    list(post.save = post.save, x.mean = x.mean, x.sd = x.sd, 
        remove.var = remove.var, y.mean = y.mean)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
