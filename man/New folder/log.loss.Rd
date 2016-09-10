\name{log.loss}
\alias{log.loss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
log.loss(p.hat, y.test, cut = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p.hat}{
%%     ~~Describe \code{p.hat} here~~
}
  \item{y.test}{
%%     ~~Describe \code{y.test} here~~
}
  \item{cut}{
%%     ~~Describe \code{cut} here~~
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
function (p.hat, y.test, cut = 0.01) 
{
    for (i in 1:length(p.hat)) {
        if (p.hat[i] < cut) {
            p.hat[i] = cut
        }
        else if (p.hat[i] > 1 - cut) {
            p.hat[i] = 1 - cut
        }
    }
    -mean(ifelse(y.test, log(p.hat), log(1 - p.hat)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
