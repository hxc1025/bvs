\name{sim.data.1}
\alias{sim.data.1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate Simulation Data of Simulation 1
}
\description{
Generate simulation data of simulation 1 in paper. (Fan and Li data).
}
\usage{
sim.data.1(n = 40, sigma = 3, beta = c(3, 1.5, 0, 0, 2, 0, 0, 0), rho = 0.5, Sigma = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Sample size.
}
  \item{sigma2}{
Standard deviation of the error term.
}
  \item{beta}{
True coefficients.
}
  \item{rho}{
Correlation coefficient. Check Fan and Li for more details. 
}
  \item{Sigma}{
Covariance matrix of the features.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A data frame containing:
  \item{x}{Design matrix.}
  \item{y}{response variable.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Xichen
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
train.data = sim.data.1()
x = train.data$x
y = train.data$y
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
