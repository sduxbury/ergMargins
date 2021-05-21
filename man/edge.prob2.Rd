\name{edge.prob2}
\alias{edge.prob2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to compute tie probabilities for each dyad.
%%  ~~function to do ... ~~
}
\description{Computes tie probabilities for each dyad. Largely a wrapper for edgeprob() in btergm.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
edge.prob2(model,verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{ model is the ergm object.
%%     ~~Describe \code{x} here~~
}
  \item{verbose}{Boolean operator. If set to T, updates user on progress of edge.prob2. May be useful for large networks.}
}
\details{
Computes tie probabilities for each dyad in the ERGM sample space.This function is largely a wrapper for the edgeprob function in btergm, but can be applied to curved ERGMs.
}
\value{
Returns a dataframe where weach row is a dyad and each column is a term in the model, an identification variable for the vertices/dyads, or the estimated probability of observing a tie. See edgeprob() in btergm for more details.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\author{
Scott Duxbury, Assistant Professor, University of North Carolina--Chapel Hill, Department of Sociology.
original author is Philip Leifeld at University of Essex.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[btergm]{edgeprob}}
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
edge.prob2(a)

}
\donttest{
library(ergm)


data("faux.dixon.high")


my.ergm<-ergm(faux.dixon.high~edges+
                nodeifactor("sex")+
                nodeofactor("sex")+
                nodematch("sex")+
                mutual)
edge.prob2(my.ergm)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~ERGM }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~tie probabilities }% __ONLY ONE__ keyword per line
