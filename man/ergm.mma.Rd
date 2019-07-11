\name{ergm.mma}
\alias{ergm.mma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to conduct marginal effects mediation analysis in ERGM.
}
\description{
ergm.mma is a function to conduct mediation analysis in ERGM using marginal effects. It computes the total, direct, and indirect average marginal effects for two nested ERGMs fit to the same data. Can be used with discrete or continuous mediators and for multiple mediator analyses.}
\usage{
ergm.mma(restricted.model,full.model,direct.effect,mediator)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{restricted.model}{the ergm without the mediator variable.
%%     ~~Describe \code{x} here~~
}
  \item{full.model}{the ergm including the mediator variable.}

  \item{mediator}{a character string containing the name of the mediator variable. When examining joint mediation with multiple mediators, this is a vector with each entry containing the character strings for each mediator variable.}
  \item{direct.effect}{a character string containing the name of the direct effect.}
}
\details{
ergm.mma conducts marginal effects mediation analysis for nested ERGMs fit to the same network by computing the change in average marginal effects between models. The AMEs are estimated on the scale of tie probabilities, not the linear component (log-odds). Standard errors are computed using the Delta method.

A few notes on using ergm.mma:
\itemize{

\item ergm.mma relies on AMEs, which can be computationally intensive in large networks or high dimensional models.

\item ergm.mma can be applied to discrete mediators or multiple mediator analysis.


}

}
\value{
Returns a table including the point estimates and standard errors for total, direct, and indirect effects. A single indirect effect is reported for joint mediation analyses. All AMEs are computed on the scale of tie probabilities, not the linear component (log-odds).

}
\references{
Duxbury, Scott W. "Mediation and Moderation in Statistical Network Models." Working paper.

Barrett, Tyson. "Marginal Mediation: A New Framework for Interpretable Mediation Effects." Dissertation.

Karlson, Kristian Bernt, Holm, Anders, and Richard Breen. 2012. "Comparing Regression Coefficients Between Same-sample Nested Models Using Logit and Probit: A New Method." \emph{Sociological Methodology} 42: 286 - 313.
}
\author{
Scott Duxbury, PhD Candidate, the Ohio State University, Department of Sociology.
}
\note{
Note that ergm.mma may have restrictive run times in large networks or high dimensional models.

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ergm.AME}}
}
\examples{
\donttest{
library(ergm)


data("faux.mesa.high")

faux.mesa.high

set.seed(21093)

my.ergm1<-ergm(faux.mesa.high~edges+
            nodecov("Grade")+
            nodefactor("Race")+
            nodefactor("Sex")+
            nodematch("Grade")+
            nodematch("Sex")+
            nodematch("Race"),
            estimate = "MLE")


set.seed(21093)

my.ergm2<-ergm(faux.mesa.high~edges+
                nodecov("Grade")+
                nodefactor("Race")+
                nodefactor("Sex")+
                nodematch("Grade")+
                nodematch("Sex")+
                nodematch("Race")+
               gwesp(.5,fixed=TRUE))


###gwesp as single mediator
ergm.mma(my.ergm1,my.ergm2,mediator="gwesp.fixed.0.5",
             direct.effect="nodecov.Grade")


###multiple mediators

set.seed(21093)

my.ergm3<-ergm(faux.mesa.high~edges+
                 nodecov("Grade")+
                 nodefactor("Race")+
                 nodefactor("Sex")+
                 nodematch("Grade")+
                 nodematch("Sex")+
                 nodematch("Race")+
                 gwdsp(1,fixed=TRUE)+
                 gwesp(.5,fixed=TRUE),
                 #increase MC sample to facilitate convergence
               control = control.ergm(MCMC.samplesize=5000),
               constraints = ~bd(maxout=5)
)

##joint mediation for gwdsp and gwesp
ergm.mma(my.ergm2,my.ergm3,mediator=c("gwdsp.fixed.1","gwesp.fixed.0.5"),
         direct.effect=names(my.ergm3$coef[8]))
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~ergm }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~networks }% __ONLY ONE__ keyword per line
\keyword{ ~mediation }% __ONLY ONE__ keyword per line
\keyword{ ~marginal effects }% __ONLY ONE__ keyword per line
