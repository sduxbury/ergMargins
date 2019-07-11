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
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a dataframe where weach row is a dyad and each column is a term in the model, an identification variable for the vertices/dyads, or the estimated probability of observing a tie. See edgeprob() in btergm for more details.
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
Scott Duxbury, original author is Philip Leifeld at University of Essex.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[btergm]{edgeprob}}
}
\examples{
library(statnet)


data("faux.dixon.high")


set.seed(21093)

my.ergm<-ergm(faux.dixon.high~edges+
                nodeicov("grade")+
                nodeocov("grade")+
                nodeifactor("race")+
                nodeofactor("race")+
                nodeifactor("sex")+
                nodeofactor("sex")+
                absdiff("grade")+
                nodematch("sex")+
                mutual,
              ergm.control=control.ergm(parallel=4,
                                        parallel.type="PSOCK"))
edge.prob2(my.ergm)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~ERGM }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~tie probabilities }% __ONLY ONE__ keyword per line
