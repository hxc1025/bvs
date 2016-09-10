\name{rMVN}
\alias{rMVN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rMVN(n = 1, mu, Sigma, L = NULL, U = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
}
  \item{Sigma}{
%%     ~~Describe \code{Sigma} here~~
}
  \item{L}{
%%     ~~Describe \code{L} here~~
}
  \item{U}{
%%     ~~Describe \code{U} here~~
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
function (n = 1, mu, Sigma, L = NULL, U = NULL) 
{
    p = length(mu)
    if (is.null(L) && is.null(U)) {
        L = t(chol(Sigma))
    }
    else if (!is.null(L)) {
        L = L
    }
    else if (!is.null(U)) {
        L = t(U)
    }
    x = mu + L \%*\% matrix(rnorm(n * p), p, n)
    t(x)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
