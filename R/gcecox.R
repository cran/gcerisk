#' Fit Generalized Competing Event Model Based on Proportional Hazards Regression
#'
#' @description Fit a generalized competing event model by using Cox proportational hazards regression model
#' with \code{coxph} function in \code{survival} package.
#' @param formula1 a formula object for event(s) of interest, with a survival response returned by \code{Surv}
#' function on the left, and the covariate terms on the right.
#' @param formula2 a formula object for event(s) of interest, with a survival response returned by \code{Surv}
#' function on the left, and the covariate terms on the right.
#' @param formula3 a formula object for competing event(s), with a survival response returned by \code{Surv}
#' function on the left, and the covariate terms on the right.
#' @param formula4 a formula object for the composite set of all events, with a survival response returned by \code{Surv}
#' function on the left, and the covariate terms on the right.
#' @param survtime survival time for event(s) of interest.
#' @param surv1 a formula object for event(s) of interest, with a survival response returned by \code{Surv}
#' function on the left, and 1 on the right.
#' @param surv2 a formula object for competing event(s), with a survival response returned by \code{Surv}
#' function on the left, and 1 on the right.
#' @param data a data frame containing variables named in formula.
#' @param ca the status indicator, normally 0 = alive, 1 = dead from the primary event(s) of interest.
#' @param cm the status indicator, normally 0 = alive, 1 = dead from the competing event(s) of interest.
#' @param all the status indicator, normally 0 = alive, 1 = dead from the all kind of event(s) of interest.
#' @param N the number of bootstrap replicates
#' @param M the number of bins for \eqn{\omega} or \eqn{\omega+} plots.
#' @param t survival time point for \eqn{\omega} or \eqn{\omega+} plots.
#' @details The \strong{gcerisk} package is designed to help investigators optimize risk-stratification methods for competing risks data, such as described in
#' Carmona R, Gulaya S, Murphy JD, Rose BS, Wu J, Noticewala S, McHale MT, Yashar CM, Vaida F, Mell LK. Validated competing event model for the stage I-II endometrial cancer population.
#' Int J Radiat Oncol Biol Phys. 2014;89:888-98. Standard risk models typically estimate the effects of one or more covariates on either
#' a single event of interest (such as overall mortality, or disease recurrence), or a composite set of events (e.g., disease-free survival, which combines events of interest with death from any cause).
#' This method is inefficient in stratifying patients who may be simultaneously at high risk for the event of interest but low risk for competing events, and who thus stand to gain the most from strategies to modulate the event of interest.
#' Compared to standard risk models, GCE models better stratify patients at higher (lower) risk for an event of interest and lower (higher) risk of competing events. GCE models focus on differentiating subjects based on
#' the ratio of the cumulative hazard (or cumulative hazard of the subdistribution) for the event of interest to the cumulative hazard (or cumulative hazard of the subdistribution) for all events (\eqn{\omega}),
#' and the ratio of the cumulative hazard (or cumulative hazard of the subdistribution) for the event of interest to the cumulative hazard (or cumulative hazard of the subdistribution) for competing events (\eqn{\omega+}).
#'
#' The \code{gcecox} function produces model estimates and confidence intervals from a generalized competing event model based on the Cox PH model for cause-specific hazards. The model assumes proportional hazards for the composite set of events.
#'
#' The function returns \eqn{\omega} and \eqn{\omega+} ratio estimates for the chosen covariates, with 95\% confidence intervals, and plots \eqn{\omega} and \eqn{\omega+} at time t within M ordered subsets of subjects as a function of increasing risk (based on the linear predictor, i.e. the inner product of a subject's data vector and the coefficient vector).
#' @import survival cmprsk ggplot2 stats
#' @examples
#' # sample data to test
#' data(Sample)
#' test <- Sample
#' rm(list=setdiff(ls(), "test"))
#' test <- transform(test, LRF_OR_DF_FLAG = as.numeric(test$LRFFLAG | test$DFFLAG))
#' test <- transform(test, CMFLAG = as.numeric(test$OSFLAG & !test$LRFFLAG & !test$DFFLAG))
#' test <- transform(test, ACMFLAG = as.numeric(test$LRF_OR_DF_FLAG | test$CMFLAG))
#'
#' formula1 <- Surv(survtime, status) ~ age + smoke20 + etohheavy + BMI +
#' dage + dsmoke20 + detohheavy + dBMI + strata(fail)
#' formula2 <- Surv(OSMO, LRF_OR_DF_FLAG) ~ age + smoke20 + etohheavy + BMI
#' formula3 <- Surv(OSMO, CMFLAG) ~ age + smoke20 + etohheavy + BMI
#' formula4 <- Surv(OSMO, ACMFLAG) ~ age + smoke20 + etohheavy + BMI
#' surv1 <- Surv(OSMO, LRF_OR_DF_FLAG) ~ 1
#' surv2 <- Surv(OSMO, CMFLAG) ~ 1
#' survtime <- test$OSMO
#' ca <- test$LRF_OR_DF_FLAG
#' cm <- test$CMFLAG
#' all <- test$ACMFLAG
#'
#' N <- 100
#' M <- 5
#' t <- 60
#'
#' fit <- gcecox(formula1, formula2, formula3, formula4, survtime,
#' surv1, surv2, test, ca, cm, all, N, M, t)

#' @author Hanjie Shen, Ruben Carmona, Loren Mell
#' @references
#' \itemize{
#' \item Carmona R, Gulaya S, Murphy JD, Rose BS, Wu J, Noticewala S, McHale MT, Yashar CM, Vaida F, Mell LK. (2014) Validated competing event model for the stage I-II endometrial cancer population. Int J Radiat Oncol Biol Phys.89:888-98.
#' \item Carmona R, Green GB, Zakeri K, Gulaya S, Xu B, Verma R, Williamson C, Rose BS, Murphy JD, Vaida F, Mell LK. (2015) Novel method to stratify elderly patients with head and neck cancer. J Clin Oncol 33 (suppl; abstr 9534).
#' \item Carmona R, Zakeri K, Green GB, Triplett DP, Murphy JD, Mell LK. (2015) Novel method to stratify elderly patients with prostate cancer. J Clin Oncol 33 (suppl; abstr 9532).
#' }
#' @return
#' \item{$coef1}{generalized competing event model coefficients (log (\eqn{\omega} ratio))}
#' \item{$coef2}{generalized competing event model coefficients (log (\eqn{\omega+} ratio))}
#' \item{$result1}{result table for generalized competing event model containing exponential of coefficients (\eqn{\omega} ratio) and 95\% confidence intervals}
#' \item{$result2}{result table for generalized competing event model containing exponential of coefficients (\eqn{\omega+} ratio) and 95\% confidence intervals}
#' \item{$omegaplot1}{\eqn{\omega} plot for generalized  competing evet model}
#' \item{$omegaplot2}{\eqn{\omega+} plot for generalized  competing evet model}
#' \item{$omegaplot3}{plot of \eqn{\omega} vs time}
#' \item{$omega}{predicted \eqn{\omega} values}
#' \item{$omegaplus}{predicted \eqn{\omega+} values}
#' @export


#### gce function by cox
gcecox <- function(formula1, formula2, formula3, formula4, survtime, surv1, surv2, data, ca, cm, all, N, M, t)
{
  # coefficients
  covnames = attr(terms(formula1), "term.labels")[1:((length(attr(terms(formula1), "term.labels"))-1)/2)]

  data$id <- 1:dim(data)[1]
  dat0 <- data[,covnames]
  dat <- data.frame(id=data$id,dat0)
  dat$survtime <- survtime
  dat$survival <- all
  dat$survival2 <- 0
  dat$ca <- ca
  dat$ca2 <- cm
  dat$ca[dat$ca != 1 & dat$ca2 != 1] <- 1

  dat2 <- reshape(dat, varying = list(c((3+length(covnames)):(4+length(covnames))),c((5+length(covnames)):(6+length(covnames)))),
                  idvar = "id",  direction = "long",
                  v.names=c("status","fail"))
  dat2 <- dat2[order(dat2$id),]

  dcovnames <- NULL
  for (k in 1:length(covnames)){
    dcovnames <- paste("d",covnames[k],sep = "")
    dat2[dcovnames] <- dat2$fail*dat2[,covnames[k]]
  }

  fit <- coxph(formula1, data = dat2)
  Beta12 <- coef(fit)[(length(covnames)+1):length(coef(fit))]

  fit.ca <- coxph(formula2, data = data)
  fit.all <- coxph(formula4, data = data)

  Beta.ca <- coef(fit.ca)
  Beta.all <- coef(fit.all)
  Betanew <- Beta.ca - Beta.all

  # variance
  var.ca <- diag(fit.ca$var)
  var.all <- diag(fit.all$var)

  coefvar1.boot <- matrix(nrow = N,ncol = length(covnames), data = NA)
  colnames(coefvar1.boot) <- covnames
  for (j in 1:N){
    for (i in 1:length(covnames)){
      Beta.ca.dist <- sample(rnorm(1000, Beta.ca[i],sqrt(var.ca)[i]), 1000, replace = T)
      Beta.all.dist <- sample(rnorm(1000, Beta.all[i],sqrt(var.all)[i]), 1000, replace = T)
      Betanew.dist <- Beta.ca.dist - Beta.all.dist
      coefvar1.boot[j,i] <- var(Betanew.dist)
    }
  }
  coefvar1 <- matrix(nrow = length(covnames),ncol = 1, data = NA)
  rownames(coefvar1) <- covnames
  colnames(coefvar1) <- " "
  for (i in 1:length(covnames)){
    coefvar1[i,] <- mean(coefvar1.boot[,i])
  }

  coefvar2 = diag(vcov(fit))[(length(covnames)+1):length(coef(fit))]


  # 95% CI
  ci1 <- matrix(NA,length(covnames),2)
  colnames(ci1) <- c("2.5%","97.5%")
  rownames(ci1) <- covnames
  df <- dim(data)[1]
  for (i in 1:length(covnames)) {
    ci1[i,1] <- (Betanew[i]+c(-1,1)*qt(0.975, df = df)*sqrt(coefvar1[i]))[1]
    ci1[i,2] <- (Betanew[i]+c(-1,1)*qt(0.975, df = df)*sqrt(coefvar1[i]))[2]
  }

  ci2 <- matrix(NA,length(covnames),2)
  colnames(ci2) <- c("2.5%","97.5%")
  rownames(ci2) <- covnames
  df <- dim(data)[1]
  for (i in 1:length(covnames)) {
    ci2[i,1] <- (Beta12[i]+c(-1,1)*qt(0.975, df = df)*sqrt(coefvar2[i]))[1]
    ci2[i,2] <- (Beta12[i]+c(-1,1)*qt(0.975, df = df)*sqrt(coefvar2[i]))[2]
  }


  # Omega and Omega plus plots
  riskscorenew <- numeric(dim(data)[1])
  for (i in 1:length(covnames)){
    riskscorenew <- riskscorenew + Betanew[i]*with(data,get(covnames[i]))
  }
  riskscore12 <- numeric(dim(data)[1])
  for (i in 1:length(covnames)){
    riskscore12 <- riskscore12 + Beta12[i]*with(data,get(covnames[i]))
  }

  normCER1 <- (riskscorenew - mean(riskscorenew))/sd(riskscorenew)
  normCER2 <- (riskscore12 - mean(riskscore12))/sd(riskscore12)


  om.fit1 <- summary(survfit(surv1, data))
  om.fit2 <- summary(survfit(surv2, data))
  surv.ca <- om.fit1$surv[om.fit1$time >= 0  & om.fit1$time <= t]
  surv.cm <- om.fit2$surv[om.fit2$time >= 0  & om.fit2$time <= t]
  hca <- -log(surv.ca)
  hcm <- -log(surv.cm)
  w0plus <- mean(hca)/mean(hcm)
  wplus <- w0plus*exp(riskscore12)
  w <- wplus/(1+wplus)

  z1 <- qplot(quantile(normCER1,prob=seq(0,1,1/M))[2:(M+1)],quantile(w,prob=seq(0,1,1/M))[2:(M+1)], xlab = "Risk Score", ylab = expression(omega))
  z2 <- qplot(quantile(normCER2,prob=seq(0,1,1/M))[2:(M+1)],quantile(wplus,prob=seq(0,1,1/M))[2:(M+1)], xlab = "Risk Score", ylab = expression(omega^"+"))


  # Omega vs Time plot
  omegas <- numeric(t)

  for (i in 1:t){
    om.fit1 <- summary(survfit(surv1, data))
    H.hat.cancer.death <- -log(om.fit1$surv)
    yca <- H.hat.cancer.death[is.finite(H.hat.cancer.death)]
    xca <- om.fit1$time[is.finite(H.hat.cancer.death)]
    fitlm <- lm(yca~xca)
    point <- data.frame(xca = i)
    H.hat.catime <- predict(fitlm, point, interval ="prediction")[1]

    om.fit2 <- summary(survfit(surv2,data))
    H.hat.competing.event <- -log(om.fit2$surv)
    ycm <- H.hat.competing.event[is.finite(H.hat.competing.event)]
    xcm <- om.fit2$time[is.finite(H.hat.competing.event)]
    fitlm2 <- lm(ycm~xcm )
    point2 <- data.frame(xcm = i)
    H.hat.cmtime <- predict(fitlm2, point2, interval ="prediction")[1]


    omegas[i] <- H.hat.catime/(H.hat.catime + H.hat.cmtime)
  }

  y3 <- omegas
  x3 <- 1:t
  z3 <- qplot(x3,y3, ylim = c(0,2*max(y3)), xlab = "Time", ylab = expression(omega))

  # Result tables
  table1 <- matrix(0,length(covnames),3)
  rownames(table1) <- covnames
  colnames(table1) <- c("exp(coef)", "lower .95", "upper .95")
  table1[,1] <- round(exp(Betanew),5)
  table1[,2] <- round(exp(ci1),5)[,1]
  table1[,3] <- round(exp(ci1),5)[,2]

  table2 <- matrix(0,length(covnames),3)
  rownames(table2) <- covnames
  colnames(table2) <- c("exp(coef)", "lower .95", "upper .95")
  table2[,1] <- round(exp(Beta12),5)
  table2[,2] <- round(exp(ci2),5)[,1]
  table2[,3] <- round(exp(ci2),5)[,2]

  return(list(coef1 = Betanew, coef2 = Beta12, result1 = table1, result2 = table2, omegaplot1 = z1,
              omegaplot2 = z2, omegaplot3 = z3, omega = w, omegaplus = wplus))

}

