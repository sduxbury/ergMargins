\name{ergm.mod.mma}
\alias{ergm.mod.mma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to assess mediated moderation using marginal effects in ERGM.
%%  ~~function to do ... ~~
}
\description{
ergm.mod.mma is a function to asssess mediated moderation in ERGM using average marginal effects (AME), e.g. when the main effect is a part of an interaction. If joint=TRUE, the function computes the total, partial, and indirect effects for all parts of the interaction. If int.eff=TRUE, the function computes the same values, but only for the moderating effect (e.g., only the moderated portion). Both joint=TRUE and int.eff=TRUE can be set simultaneously. }
\usage{
ergm.mod.mma(restricted.model,full.model,var1, var2, inter,mediator,
                   at.2=NULL,joint=FALSE,int.eff=FALSE)}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{restricted.model}{the ergm without the mediator variable.
%%     ~~Describe \code{x} here~~
}
  \item{full.model}{the ergm including the mediator variable.}

\item{var1}{the character string containing the name of the main effect of the direct effect.}
\item{var2}{the character string containing the name of the moderator variable for the direct effect.}
\item{inter}{the character string containing the name of the interaction for the direct effect.}
\item{at.2}{the levels of var2 at which to compute the average marginal effects. If at.2=NULL, the function computes the AME at each unique value of var2. at.2 can also be set to meaningful values, such as the mean of var2 or other discrete values of interest. If var2 is discrete, leaving at.2 NULL will compute the AME for each discrete category of the moderator variable. It is usually safe to leave at.2=NULL unless examining a large network where var2 has a large number of unique values, in which case computation times can be restrictive.}
  \item{mediator}{a character string containing the name of the mediator variable. When examining joint mediation with multiple mediators, this is a vector with each entry containing the character strings for each mediator variable.}
  \item{joint}{logical parameter indicates whether to compute the mediated composite effect.}
  \item{int.eff}{logical parameter indicating whether to compute the mediated effect \emph{only} for the interaction, net of the main effect. Either int.eff or joint must be set to TRUE for function to work.}

}
\details{
ergm.mod.mma conducts mediation analyses using marginal effects when the main effect is part of an interaction. There are three ways to construe the mediating pathways. First, a researcher may be interested in whether a mediator variable mediates the main effect in an interaction, net of the moderator. To assess this, use ergm.mma for var1. This will compute the AME when the effect of the moderator is zero. Second, a researcher may be interested in the total or composite effect of an interaction, including the moderator and main effect. To assess this, set joint=TRUE. This computes the indirect effect for the AME at every level of at.2. The third scenario is when a researcher is interested in only the interaction effect. The interaction effect is defined here as the difference in AMEs between levels of a moderator (e.g., levels of at.2). Thus, the indirect effect is the third difference, or the difference in second differences between the restricted and full models. An example is whether triad closure mediates the effect of matched nodal characteristics, net of the main effect for nodal attributes (e.g., net of nodefactor/nodecov). To assess this, set int.eff=TRUE. Both joint=TRUE and int.eff=TRUE can be specified at the same time, but at least one must be specified for the function to work.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a list object containing the total, partial, and indirect AMEs at each level of at.2. If joint=TRUE is specified, returns:
\item{marginal.effects}{a list containing the indirect, total, and partial effects for each level of at.2}
\item{summary.stats}{a table reporting the mean second difference and mean absolute Wald statistic for the indirect AMEs.}

If int.eff=TRUE is specified, returns a list with elements:
\item{partial.marginal.effects}{the marginal effects for each level of at.2 in the full model.}
\item{total.marginal.effects}{the marginal effects for each level of at.2 in the restricted model.}
\item{partial.second.diffs}{the second differences between levels of at.2 in the full model.}
\item{total.second.diffs}{the second differences between levels of at.2 in the restricted model.}
\item{third.diffs}{the difference in second differences between the partial and full models.}
\item{summary.output}{the mean third difference and mean absolute Wald statistic for the third differences.}

If both int.eff=TRUE and joint=TRUE, returns a list with elements:
\item{joint.effect}{a list containing all elements when joint=TRUE}
\item{moderator.effect}{a list containing all elements when int.eff=TRUE.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Duxbury, Scott W. "Mediation and Moderation in Statistical Network Models." Working paper.
%% ~put references to the literature/web site here ~
}
\author{
Scott Duxbury, PhD Candidate, the Ohio State University, Department of Sociology.

%%  ~~who you are~~
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ergm.AME}}
\code{\link{ergm.MEM}}
\code{\link{ergm.mma}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
v1<-rbinom(10,1,.3)
v2<-rbinom(10,1,.3)

network::set.vertex.attribute(a1,"var.1",v1)
network::set.vertex.attribute(a1,"var.2",v2)


a<-ergm(a1~edges+nodeifactor("var.1")+nodematch("var.1")+nodeofactor("var.1"))
ab<-ergm(a1~edges+nodeifactor("var.1")+nodematch("var.1")+nodeofactor("var.1")+
           nodeifactor("var.2"))

ergm.mod.mma(a,ab,var1="nodeifactor.var.1.1",var2="nodeofactor.var.1.1",
            inter="nodematch.var.1",mediator="nodeifactor.var.2.1",int.eff=TRUE)
}
\donttest{

library(ergm)


data("faux.dixon.high")

faux.dixon.high


my.ergm<-ergm(faux.dixon.high~edges+
                nodeifactor("sex")+
                nodeofactor("sex")+
                nodematch("sex"))

my.ergmb<-ergm(faux.dixon.high~edges+
                nodeifactor("sex")+
                nodeofactor("sex")+
                nodematch("sex")+
                mutual)

summary(my.ergmb)


sex.mod.med<-ergm.mod.mma(my.ergm,my.ergmb,var1="nodeifactor.sex.2",var2="nodeofactor.sex.2",
             inter="nodematch.sex",mediator="mutual",int.eff = TRUE,joint=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~ergm }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~networks }% __ONLY ONE__ keyword per line
\keyword{ ~mediation }% __ONLY ONE__ keyword per line
\keyword{ ~marginal effects }% __ONLY ONE__ keyword per line
\keyword{ ~moderation }% __ONLY ONE__ keyword per line
\keyword{ ~mediated moderation }% __ONLY ONE__ keyword per line
