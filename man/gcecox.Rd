\name{gcecox}
\alias{gcecox}
\title{
Fit Generalized Competing Event Model Based on Proportional Hazards Regression
}
\usage{
gcecox(Time, Ind, Cov, N, M, t)
}
\arguments{
  \item{Time}{
survival time for event(s) of interest.
}
  \item{Ind}{
the status indicators including the primary event(s) of interest, competing event(s) of interest, and all kind of event(s) of interest, normally 0 = alive, 1 = dead from the specific event(s) of interest.
}
  \item{Cov}{
a data frame containing all covariates.
}
\item{N}{the number of bootstrap replicates}

\item{M}{the number of bins for \eqn{\omega} or \eqn{\omega+} plots.}

\item{t}{survival time point for \eqn{\omega} or \eqn{\omega+} plots.}
}

\value{
\item{$coef1}{generalized competing event model coefficients (log (\eqn{\omega} ratio))}
\item{$coef2}{generalized competing event model coefficients (log (\eqn{\omega+} ratio))}
\item{$result1}{result table for generalized competing event model containing exponential of coefficients (\eqn{\omega} ratio) and 95\% confidence intervals}
\item{$result2}{result table for generalized competing event model containing exponential of coefficients (\eqn{\omega+} ratio) and 95\% confidence intervals}
\item{$omegaplot1}{\eqn{\omega} plot for generalized  competing event model}
\item{$omegaplot2}{\eqn{\omega+} plot for generalized  competing event model}
\item{$omegaplot3}{plot of \eqn{\omega} vs time}
\item{$omega}{predicted \eqn{\omega} values}
\item{$omegaplus}{predicted \eqn{\omega+} values}
\item{$riskscore1}{predicted risk scores for \eqn{\omega}}
\item{$riskscore2}{predicted risk scores for \eqn{\omega+}}

}

\description{
Fit a generalized competing event model by using Cox proportational hazards regression model
with \code{coxph} function in \code{survival} package.
}
\details{
The \strong{gcerisk} package is designed to help investigators optimize risk-stratification methods for competing risks data, such as described in
Carmona R, Gulaya S, Murphy JD, Rose BS, Wu J, Noticewala S, McHale MT, Yashar CM, Vaida F, Mell LK. Validated competing event model for the stage I-II endometrial cancer population.
Int J Radiat Oncol Biol Phys. 2014;89:888-98. Standard risk models typically estimate the effects of one or more covariates on either
a single event of interest (such as overall mortality, or disease recurrence), or a composite set of events (e.g., disease-free survival, which combines events of interest with death from any cause).
This method is inefficient in stratifying patients who may be simultaneously at high risk for the event of interest but low risk for competing events, and who thus stand to gain the most from strategies to modulate the event of interest.
Compared to standard risk models, GCE models better stratify patients at higher (lower) risk for an event of interest and lower (higher) risk of competing events. GCE models focus on differentiating subjects based on
the ratio of the cumulative hazard (or cumulative hazard of the subdistribution) for the event of interest to the cumulative hazard (or cumulative hazard of the subdistribution) for all events (\eqn{\omega}),
and the ratio of the cumulative hazard (or cumulative hazard of the subdistribution) for the event of interest to the cumulative hazard (or cumulative hazard of the subdistribution) for competing events (\eqn{\omega+}).

The \code{gcecox} function produces model estimates and confidence intervals from a generalized competing event model based on the Cox PH model for cause-specific hazards. The model assumes proportional hazards for the composite set of events.

The function returns \eqn{\omega} and \eqn{\omega+} ratio estimates for the chosen covariates, with 95\% confidence intervals, and plots \eqn{\omega} and \eqn{\omega+} at time t within M ordered subsets of subjects as a function of increasing risk (based on the linear predictor, i.e. the inner product of a subject's data vector and the coefficient vector).
}

\examples{
# sample data to test
data(Sample)
test <- Sample
rm(list=setdiff(ls(), "test"))
test <- transform(test, LRF_OR_DF_FLAG = as.numeric(test$LRFFLAG | test$DFFLAG))
test <- transform(test, CMFLAG = as.numeric(test$OSFLAG & !test$LRFFLAG & !test$DFFLAG))
test <- transform(test, ACMFLAG = as.numeric(test$LRF_OR_DF_FLAG | test$CMFLAG))

Time <- test$OSMO/12
Ind <- data.frame(test$LRF_OR_DF_FLAG, test$CMFLAG, test$ACMFLAG)
Cov <- test[,c(3,4,6,15)]
N <- 100
M <- 5
t <- 5

fit <- gcecox(Time, Ind, Cov, N, M, t)

}
\references{
\itemize{
\item Carmona R, Gulaya S, Murphy JD, Rose BS, Wu J, Noticewala S, McHale MT, Yashar CM, Vaida F, Mell LK. (2014) Validated competing event model for the stage I-II endometrial cancer population. Int J Radiat Oncol Biol Phys.89:888-98.
\item Carmona R, Green GB, Zakeri K, Gulaya S, Xu B, Verma R, Williamson C, Rose BS, Murphy JD, Vaida F, Mell LK. (2015) Novel method to stratify elderly patients with head and neck cancer. J Clin Oncol 33 (suppl; abstr 9534).
\item Carmona R, Zakeri K, Green GB, Triplett DP, Murphy JD, Mell LK. (2015) Novel method to stratify elderly patients with prostate cancer. J Clin Oncol 33 (suppl; abstr 9532).
}
}
\author{
Hanjie Shen, Ruben Carmona, Loren Mell
}

