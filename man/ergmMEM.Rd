\name{ergm.MEM}
\alias{ergm.MEM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute marginal effects at means in ERGM.
}
\description{
ERGM function to compute marginal effects at means for main effects or at various levels of a moderator. Also provides tests of significance for second differences for interaction terms. Syntax similar to margins by Thomas Leeper and Stata margins command. See ergm.AME to compute average marginal effect.
}
\usage{
ergm.MEM(model,var1,var2=NULL,inter=NULL,at.2=NULL, return.dydx=F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{the ergm object to be analyzed.
%%     ~~Describe \code{x} here~~
}
\item{var1}{the character string containing the name of the focal variable on which to compute the MEM. If var2 is NULL, this is the simple MEM.}
\item{var2}{the character string containing the name of the moderator variable. Default value is NULL.}
\item{inter}{the character string containing the name of the interaction between var1 and var2. Default value is NULL.}
\item{at.2}{the levels of var2 at which to compute the MEMs. If var2=NULL, this parameter is ignored. If var2!=NULL and at.2=NULL, the function computes the AME at each unique value of var2. at.2 can also be set to meaningful values, such as the mean of var2 or other discrete values of interest. If var2 is discrete, leaving at.2 NULL will compute the AME for each discrete category of the moderator variable. It is usually safe to leave at.2=NULL unless examining a large network where var2 has a large number of unique values, in which case computation times can be restrictive.}
\item{return.dydx}{logical parameter indicating whether to return the marginal effects for each dyad observation. Defaults to FALSE.}
}
\details{
Largely a wrapper for ergm.AME. The only difference is that marginal effects are computed at variable means.

Computes the MEM for variables specified. If var2 and inter are left at NULL, this is the simple MEM. All MEMs are computed as tie probabilities (not on the linear scale). Standard errors are computed using the Delta method. If a moderator is specified using var2 and inter, a table of MEMs is returned at each level of at.2. If at.2=NULL, this is computed for every unique level of var2.

When at.2 and var2 are specified, a table of second differences is also returned. The table of second differences assesses whether the difference in marginal effects is significant between levels of a moderator (values of at.2). If the moderator is binary, the second difference tests the overall significance of an interaction. The tests statistic is a Wald test.

If more than 2 values of at.2 are specified, the function also computes the average absolute Wald test statistic and average second difference. The average second difference measures the average change in the MEM of var1 between levels of a moderator. The average absolute Wald test is the mean absolute value of the Wald tests for all second differences. It assesses whether, on average, there is a significant difference between the second differences and zero. A significant result means the interaction is significant.
}
\value{
If var2 and inter=NULL, then a table is returned with the MEM for var1, its standard error, and test statistic. Note again that all MEMs are in the metric of tie probabilities, not log-odds. If a moderator is specified, a list is returned with the following elements:
\item{`Marginal effects at means`}{is the table of MEMs at each level of at.2. If at.2=NULL, the table contains unique entries for every unique value of at.2.}
\item{`Second differences`}{is a table containing point estimates, standard errors, and Wald tests for the difference in MEMs between each level of at.2. It has a number of rows equal to length(at.2)-1}
\item{`Aggregate output`}{is a table of summary statistics when at.2 has length > 2. Reports the mean second difference and mean absolute Wald statistic for the second differences. If at.2 has length<3, it is not returned}
\item{`Marginal effects`}{if return.dydx=T, returns the marginal effects (partial derivatives with respect to var1) for each dyad in the ERGM sample space.}

}
\references{
Duxbury, Scott W. "Mediation and Moderation in Statistical Network Models." Working paper.

Long, J. Scott, and Sarah Mustillo. 2018. "Using Predictions and Marginal Effects to Compare Groups in Regression Models for Binary Outcomes." \emph{Sociological Methods and Research}

Mize, Trenton D. 2019. "Best Practices for Estimating, Interpreting, and Presenting Nonlinear Interaction Effects." \emph{Sociological Science}
}
\author{
Duxbury, Scott W. PhD Candidate, the Ohio State University Department of Sociology.
}
\note{
Syntax is loosely based on the margins package by Thomas Leeper and Stata's margins commands. See ergm.AME to compute average marginal effects. Currently only considers two-ways interactions.}


\seealso{
\code{\link{ergm.AME}}
\code{\link{ergm.mma}}
\code{\link[margins]{margins}}
}
\examples{
\donttest{

library(ergm)


data("faux.dixon.high")


set.seed(21093)

my.ergm<-ergm(faux.dixon.high~edges+
                nodeicov("grade")+
                nodeocov("grade")+
                nodeifactor("sex")+
                nodeofactor("sex")+
                absdiff("grade")+
                nodematch("sex")+
                mutual)

#main effect
main.eff<-ergm.MEM(my.ergm,"mutual")


##with moderator
sex.inter<-ergm.MEM(my.ergm,"nodeifactor.sex.2","nodeofactor.sex.2","nodematch.sex")

##with moderator at specific levels
grade.inter<-ergm.MEM(my.ergm,"nodeicov.grade","nodeocov.grade","absdiff.grade",at.2=c(9,12))

##when moderator changes from mean to mean +sd
tie.probs<-edge.prob2(my.ergm)
m<-mean(tie.probs$nodeocov.grade)
m.sd<-m+sd(tie.probs$nodeocov.grade)

grade.inter.at.means<-ergm.MEM(my.ergm,"nodeicov.grade","nodeocov.grade",
                            "absdiff.grade",at.2=c(m,m.sd))

}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~ergm }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~networks }% __ONLY ONE__ keyword per line
\keyword{ ~marginal effects }% __ONLY ONE__ keyword per line
\keyword{ ~moderation }% __ONLY ONE__ keyword per line
