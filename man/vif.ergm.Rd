\name{vif.ergm}
\alias{vif.ergm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to diagnose multicollinearity in ERGM.}
\description{
Computes the variance inflation factor from the covariance matrix of the ERGM estimator.}
\usage{
vif.ergm(my.ergm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{my.ergm}{an ergm object.
%%     ~~Describe \code{x} here~~
}
}
\details{
Computes the variance inflation factor using the covariance matrix of the ERGM estimator. Values above 20 indicate the onset of moderate collinearity problems, values above 100 indicate that strong collinearity likely exists and may be problematizing the convergence of the Markov chain to the stationarity distribution. High VIFs alongside a near-singular Hessian matrix is a clear indicator of multicollinearity. Consider changing model specification in such cases.
}
\value{
Returns a matrix of VIFs for each model parameter. Higher values indicate that the variable in one column of the data matrix is highly correlated with a linear combination of the remaining model covariates. High VIFs can also arise from a large bivariate correlation.
}
\references{
Duxbury, Scott W. 2018. "Diagnosing Multicollinearity in Exponential Random Graph Models" \emph{Sociological Methods and Research.}
}
\author{
Scott W. Duxbury, PhD Candidate, The Ohio State University, Department of Sociology.
}


\examples{
\dontshow{
require(ergm)

set.seed(21093)
a1<-network::as.network(matrix(c(rbinom(10, 1,.3),
            rbinom(10, 1,.3),
             rbinom(10, 1,.3),
              rbinom(10, 1,.3),
               rbinom(10, 1,.3),
            rbinom(10, 1,.3),
            rbinom(10, 1,.3),
            rbinom(10, 1,.3),
            rbinom(10, 1,.3),
            rbinom(10, 1,.3)),
          nrow=10,ncol=10))

network::set.vertex.attribute(a1,"var.1",rbinom(10,1,.3))
a<-ergm(a1~edges+nodeifactor("var.1")+nodeofactor("var.1"))
vif.ergm(a)
}
set.seed(21093)

library(ergm)
data("faux.mesa.high")

my.ergm2<-ergm(faux.mesa.high~edges+
                nodecov("Grade")+
                nodefactor("Race")+
                nodefactor("Sex")+
                nodematch("Grade")+
                nodematch("Sex")+
                nodematch("Race"))

vif.ergm(my.ergm2)

#to compute the tolerance instead of the VIF

tol<-1/(vif.ergm(my.ergm2))

#to compute R squared
r.sq<-1-tol

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~ergm }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~networks }% __ONLY ONE__ keyword per line
\keyword{ ~collinearity}
