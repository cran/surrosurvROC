\name{surrosurvROC}
\alias{surrosurvROC}
\title{Surrogate Survival ROC}
\description{Nonparametric and semiparametric estimations of the time-dependent ROC curve for an incomplete failure time data with surrogate failure time endpoints}
\usage{
  surrosurvROC(DATA, method, pred.time, wt=NULL, span=NULL, b.rep=200)
}
\arguments{
  \item{DATA}{data frame, consisting of Marker: Predictior or marekr value; Survival time; Status: Event indicator (1: event; 0: censoring); STime: Surroagte survival Time; SStatus: Surrogate event indicator (1: event; 0: censoring)}
  \item{method}{"KNN"" for nonparametric model using nearest neighborhood kernel; "COX"" for semiparametric proportional hazard model}
  \item{pred.time}{Prediction time of the ROC curve}
  \item{wt}{Weight, such as inverse probablity weighting}
  \item{span}{Smoothing bandwidth parameter for KNN}
  \item{b.rep}{Number of bootstrap}
}
\details{It provides a more efficient time-dependent ROC curve for an incomplete failure time data, when surrogate failure time endpoints are additionally observed for all subjects.
}

\references{Yunro Chung and Yingye Zheng, Evaluating Prognostic Accuracy of Biomarkers for Incomplete and
Right-Censored Data with Surrogate Outcome (in progress)}
\author{Yunro Chung [cre]}

\examples{
DATA=data.frame(
  Time=   c(1,2,5,3,9,NA,8,9,10,NA,NA,NA,6,4,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
  Status= c(1,1,0,0,1,NA,1,1,0, NA,NA,NA,0,0,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
  STime=  c(3,2,4,2,8,5,8,7,11,1,8,9,3,5,2,5,10,3,5,8,5,2,4,6,7),
  SStatus=c(1,0,1,0,1,1,1,0,0,1,1,1,1,0,1,1,0,0,1,0,1,0,1,0,0),
  Marker= c(1,5,1,2,3,1,2,3,4,5,9,8,5,7,3,4,2,5,3,4,7,5,9,3,8)
)

#COX at year 3
RES1=surrosurvROC(DATA, method="COX", pred.time=3)
print(RES1)

#KNN at year 3
nobs=sum(!is.na(DATA$Time))
span=0.25*nobs^(-0.20)
RES2=surrosurvROC(DATA, method="KNN",pred.time=3,span=span)
print(RES2)
}

\keyword{biomarkers, inverse probability weighing, bootstrap}
