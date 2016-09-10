\name{beta.convert.back}
\alias{beta.convert.back}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
beta.convert.back(BVS, phi.cut = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BVS}{
%%     ~~Describe \code{BVS} here~~
}
  \item{phi.cut}{
%%     ~~Describe \code{phi.cut} here~~
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
function (BVS, phi.cut = 0.5) 
{
    beta.bar.temp = with(BVS, ifelse(phi > phi.cut, mu, 0))
    beta.bar.temp = with(BVS, beta.bar.temp/x.sd)
    p2 = length(BVS$remove.var)
    if (p2 != 0) {
        beta.bar = rep(NA, BVS$p + p2)
        beta.bar[-BVS$remove.var] = beta.bar.temp
    }
    else {
        beta.bar = beta.bar.temp
    }
    beta.bar
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
