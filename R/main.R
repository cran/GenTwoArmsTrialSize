#' @rdname getSize
#' @aliases getSizeMean
#' @aliases getSizeProp
#' @aliases getSizeTTE
#' @aliases getSizeOrd
#'
#' @title General Formulas for Sample Size Calculation
#' @description This function computes the sample size required for two arms clinical trials with continuous (`getSizeMean`), binary (`getSizeProp`), TTE (`getSizeTTE`), and ordinal outcome measure (`getSizeOrd`). Four hypothesis tests are available under two allocation designs.
#' @param design allocation method (`parallel` or `crossover`).
#' @param test four hypothesis tests: `equality`, `noninferiority`, `superiority`, and `equivalence`.
#' @param alpha level of significance.
#' @param beta type II error.
#' @param sigma pooled standard deviation of two groups.
#' @param k ratio of control to treatment.
#' @param delta delta margin in test hypothesis.
#' @param TTE target treatment effect or effect size.
#' @param rho vector of length 2, positive noncompliance rates of two arms.
#' @param r projected proportion of trial uniform loss of follow-up.
#' @return sample size per arm.
#' @export
getSizeMean <- function(design=c("parallel", "crossover"),
                        test=c("equality", "noninferiority", "superiority", "equivalence"),
                        alpha=0.05, beta=0.20, sigma, k=1, delta=0, TTE,
                        rho=c(0.05, 0.07), r=0.1){

  stopifnot(design %in% c("parallel", "crossover"),
            test %in% c("equality", "noninferiority", "superiority", "equivalence"),
            dplyr::between(alpha,0,1), dplyr::between(beta,0,1),
            dplyr::between(sigma,0,Inf), dplyr::between(k,0,Inf),
            dplyr::between(delta,-Inf,Inf), dplyr::between(TTE,-Inf,Inf),
            length(rho) == 2,
            dplyr::between(rho[1],0,1), dplyr::between(rho[2],0,1),
            dplyr::between(r,0,1))

  m <- (1-rho[1]-rho[2])*TTE

  if(design == "parallel" & test == "equality") n_trt <- ceiling((TrialSize::TwoSampleMean.Equality(alpha, beta, sigma, k, margin=m))/(1-r))

  if(design == "parallel" & test == "noninferiority") n_trt <- ceiling((TrialSize::TwoSampleMean.NIS(alpha, beta, sigma, k, delta, margin=m))/(1-r))

  if(design == "parallel" & test == "superiority") n_trt <- ceiling((TrialSize::TwoSampleMean.NIS(alpha, beta, sigma, k, delta, margin=m))/(1-r))

  if(design == "parallel" & test == "equivalence") n_trt <- ceiling((TrialSize::TwoSampleMean.Equivalence(alpha, beta, sigma, k, delta, margin=m))/(1-r))

  if(design == "crossover" & test == "equality") n_trt <- ceiling((TrialSize::TwoSampleCrossOver.Equality(alpha, beta, sigma, margin=m))/(1-r))

  if(design == "crossover" & test == "noninferiority" ) n_trt <- ceiling((TrialSize::TwoSampleCrossOver.NIS(alpha, beta, sigma,  delta, margin=m))/(1-r))

  if(design == "crossover" & test == "superiority" ) n_trt <- ceiling((TrialSize::TwoSampleCrossOver.NIS(alpha, beta, sigma, delta, margin=m))/(1-r))

  if(design == "crossover" & test == "equivalence" ) n_trt <- ceiling((TrialSize::TwoSampleCrossOver.Equivalence(alpha, beta, sigma, delta, margin=m))/(1-r))

  ans <- stats::setNames(
    data.frame(matrix(data = c(n_trt,k*n_trt), ncol=2, nrow=1), row.names="Size"),
    c("n_trt", "n_ctl"))
  return(ans)

}

#' @examples
#' # Ex 1. (n_trt=91, n_ctl=91)
#' getSizeMean(design="parallel", test="equality", alpha=0.05, beta=0.20,
#'   sigma=0.10, k=1, delta=0, TTE=0.05, rho=c(0.05, 0.07), r=0.1)
#'
#' getSizeMean(design="parallel", test="noninferiority", alpha=0.05,
#'  beta=0.20, sigma=0.10, k=1, delta=-0.05, TTE=0, rho=c(0.05, 0.07), r=0.1)
#'
#' # Ex 3. (n_trt=1022, n_ctl=1022)
#' getSizeMean(design="parallel", test="superiority", alpha=0.05, beta=0.20,
#'   sigma=0.10, k=1, delta=0.05, TTE=0.07, rho=c(0.05, 0.07), r=0.1)
#'
#' # Ex 4. (n_trt=113, n_ctl=113)
#' getSizeMean(design="parallel", test="equivalence", alpha=0.05, beta=0.20,
#'   sigma=0.10, k=1, delta=0.05, TTE=0.01, rho=c(0.05, 0.07), r=0.1)
#'
#' # Ex 5. (n_trt=23, n_ctl=23)
#' getSizeMean(design="crossover", test="equality", alpha=0.05, beta=0.20,
#'  sigma=0.10, k=1, delta=0, TTE=0.05, rho=c(0.05, 0.07), r=0.1)
#'
#' # Ex 6. (n_trt=14, n_ctl=14)
#' getSizeMean(design="crossover", test="noninferiority", alpha=0.05,
#'  beta=0.20, sigma=0.10, k=1, delta=-0.05, TTE=0, rho=c(0.05, 0.07), r=0.1)
#'
#' # Ex 7. (n_trt=21, n_ctl=21)
#' getSizeMean(design="crossover", test="superiority", alpha=0.05, beta=0.20,
#'  sigma=0.10, k=1, delta=0.05, TTE=0.01, rho=c(0.05, 0.07), r=0.1)
#'
#' # Ex 8. (n_trt=29, n_ctl=29)
#' getSizeMean(design="crossover", test="equivalence", alpha=0.05, beta=0.20,
#'  sigma=0.10, k=1, delta=0.05, TTE=0.01, rho=c(0.05, 0.07), r=0.1)
#'

#' @rdname getSize
#' @param varsigma (varsigma1 > 0, varsigma2 > 0) := (p1, p2) probability of mean response in control and treatment arms; (`varsigma1` > 0, `varsigma2` > 0) := (`sigma`, `sigma`) pooled standard deviation of two groups or their difference (sigma>0)
#' @param seqnumber Number of crossover sequences: 0 if parallel; 1+ if crossover (seqnumber>=0)
#' @export
getSizeProp <- function(design=c("parallel", "crossover"),
                        test=c("equality", "noninferiority", "superiority", "equivalence"),
                        alpha=0.05, beta=0.20, varsigma, k=1,
                        seqnumber, delta=0, TTE,
                        rho=c(0.05, 0.07), r=0.1){

  stopifnot(design %in% c("parallel", "crossover"),
            test %in% c("equality", "noninferiority", "superiority", "equivalence"),
            dplyr::between(alpha,0,1), dplyr::between(beta,0,1),
            dplyr::between(varsigma[1],0,Inf), dplyr::between(varsigma[2],0,Inf),
            dplyr::between(k,0,Inf), dplyr::between(seqnumber,0,Inf),
            dplyr::between(delta,-Inf,Inf), dplyr::between(TTE,-Inf,Inf),
            length(rho) == 2,
            dplyr::between(rho[1],0,1), dplyr::between(rho[2],0,1),
            dplyr::between(r,0,1))


  varsigmaa1 <- rho[1]*varsigma[2]+(1-rho[1])*varsigma[1]
  varsigmaa2 <- rho[2]*varsigma[1]+(1-rho[2])*varsigma[2]

  if(design == "parallel" & test == "equality") n_trt <- ceiling( (TrialSize::TwoSampleProportion.Equality(alpha, beta, p1=varsigmaa1, p2=varsigmaa2, k=k))/(1-r))

  if(design == "parallel" & test == "noninferiority") n_trt <- ceiling((TrialSize::TwoSampleProportion.NIS(alpha, beta, p1=varsigmaa1, p2=varsigmaa2, k=k, delta=(1-rho[1]-rho[2])*TTE, margin=delta))/(1-r))

  if(design == "parallel" & test == "superiority") n_trt <- ceiling( (TrialSize::TwoSampleProportion.NIS(alpha, beta, p1=varsigmaa1, p2=varsigmaa2, k=k, delta=(1-rho[1]-rho[2])*TTE, margin=delta))/(1-r))

  if(design == "parallel" & test == "equivalence") n_trt <- ceiling( (TrialSize::TwoSampleProportion.Equivalence(alpha, beta,  p1=varsigmaa1, p2=varsigmaa2, k, delta=(1-rho[1]-rho[2])*TTE, margin=delta))/(1-r))

  if(design == "crossover" & test == "equality") n_trt <- ceiling((TrialSize::TwoSampleSeqCrossOver.Equality(alpha, beta, sigma=(varsigma[1]**2), sequence=seqnumber, delta=(1-rho[1]-rho[2])*TTE))/(1-r))

  if(design == "crossover" & test == "noninferiority") n_trt <- ceiling((TrialSize::TwoSampleSeqCrossOver.NIS(alpha, beta, sigma=(varsigma[1]**2), sequence=seqnumber, delta=delta, margin=(1-rho[1]-rho[2])*TTE))/(1-r))

  if(design == "crossover" & test == "superiority") n_trt <- ceiling((TrialSize::TwoSampleSeqCrossOver.NIS(alpha, beta, sigma=(varsigma[1]**2), sequence=seqnumber, delta=delta, margin=(1-rho[1]-rho[2])*TTE))/(1-r))

  if(design == "crossover" & test == "equivalence") n_trt <- ceiling((TrialSize::TwoSampleSeqCrossOver.Equivalence(alpha, beta, sigma=(varsigma[1]**2), sequence=seqnumber, delta=delta, margin=(1-rho[1]-rho[2])*TTE))/(1-r))

  ans <- stats::setNames(
    data.frame(matrix(data = c(n_trt,k*n_trt),ncol=2, nrow=1), row.names="Size"),
    c("n_trt", "n_ctl"))

  return(ans)
}

#' @examples
#' # Ex 1. (n_trt=102, n_ctl=102)
#' getSizeProp(design="parallel", test="equality", alpha=0.05, beta=0.20,
#'  varsigma=c(0.65, 0.85), k=1, seqnumber=0, delta=0, TTE=0,
#'  rho=c(0.05, 0.07), r=0.1)
#'
#' # Ex 2. (n_trt=33, n_ctl=33)
#' getSizeProp(design="parallel", test="noninferiority", alpha=0.05, beta=0.20,
#'  varsigma=c(0.65,0.85), k=1, seqnumber=0, delta=-0.10, TTE=0.20,
#'  rho=c(0.05, 0.07), r=0.1)
#'
#' # Ex 3. (n_trt=157, n_ctl=157)
#' getSizeProp(design="parallel", test="superiority", alpha=0.05, beta=0.20,
#'  varsigma=c(0.65,0.85), k=1, seqnumber=0, delta=0.05, TTE=0.20,
#'  rho=c(0.05, 0.07), r=0.1)
#'
#' # Ex 4. (n_trt=137, n_ctl=137)
#' getSizeProp(design="parallel", test="equivalence", alpha=0.05, beta=0.20,
#'  varsigma=c(0.75,0.80), k=1, seqnumber=0, delta=0.20, TTE=0.05,
#'  rho=c(0.05, 0.07), r=0.1)
#'
#' # Ex 5. (n_trt=36, n_ctl=36)
#' getSizeProp(design="crossover", test="equality", alpha=0.05, beta=0.20,
#'  varsigma=c(0.5,0.5), k=1, seqnumber=2, delta=0, TTE=0.20,
#'  rho=c(0.05, 0.07), r=0.1)
#'
#' # Ex 6. (n_trt=22, n_ctl=22)
#' getSizeProp(design="crossover", test="noninferiority", alpha=0.05,
#'  beta=0.20, varsigma=c(0.5,0.5), k=1, seqnumber=2, delta=-0.20, TTE=0,
#'  rho=c(0.05, 0.07), r=0.1)
#'
#' # Ex 7. (n_trt=86, n_ctl=86)
#' getSizeProp(design="crossover", test="superiority", alpha=0.05, beta=0.20,
#'  varsigma=c(0.5,0.5), k=1, seqnumber=2, delta=0.10, TTE=0,
#'  rho=c(0.05, 0.07), r=0.1)
#'
#' # Ex 8. (n_trt=30, n_ctl=30)
#' getSizeProp(design="crossover", test="equivalence", alpha=0.05, beta=0.20,
#'  varsigma=c(0.5,0.5), k=1, seqnumber=2, delta=0.20, TTE=0,
#'  rho=c(0.05, 0.07), r=0.1)


#' @rdname getSize
#' @param varlambda (varlambda1>0,varlambda2>0):=(lam1,lam2) hazard rates in control and treatment arms
#' @param ttotal total trial time(ttoal>0)
#' @param taccrual accrual time period (taccrual>0)
#' @param gamma parameter of exponential distribution(gamma>=0)
#' @export
getSizeTTE <- function(design=c("parallel", "crossover"),
                       test=c("equality", "noninferiority", "superiority", "equivalence"),
                       alpha=0.05, beta=0.20, varlambda, k=1,
                       ttotal, taccrual, gamma, delta=0,
                       rho=c(0.05, 0.07), r=0.1){

  stopifnot(design %in% c("parallel", "crossover"),
            test %in% c("equality", "noninferiority", "superiority", "equivalence"),
            dplyr::between(alpha,0,1), dplyr::between(beta,0,1),
            dplyr::between(varlambda[1],0,Inf), dplyr::between(varlambda[2],0,Inf),
            dplyr::between(k,0,Inf), dplyr::between(ttotal,0,Inf),
            dplyr::between(taccrual,0,Inf), dplyr::between(gamma,0,Inf),
            dplyr::between(delta,-Inf,Inf), length(rho) == 2,
            dplyr::between(rho[1],0,1),
            dplyr::between(rho[2],0,1), dplyr::between(r,0,1))

  varlambdaa1 <- rho[1]*varlambda[2]+(1-rho[1])*varlambda[1]
  varlambdaa2 <- rho[2]*varlambda[1]+(1-rho[2])*varlambda[2]

  if ( design=="parallel" & test=="equality" ) n_trt <- ceiling((TrialSize::TwoSampleSurvival.Equality(alpha, beta, lam1=varlambdaa1, lam2=varlambdaa2, k, ttotal, taccrual, gamma))/(1-r))

  if ( design == "parallel" & test == "noninferiority") n_trt <- ceiling((TrialSize::TwoSampleSurvival.NIS(alpha, beta, varlambdaa1, varlambdaa2, k, ttotal, taccrual, gamma, delta))/(1-r))

  if(design == "parallel" & test == "superiority") n_trt <- ceiling((TrialSize::TwoSampleSurvival.NIS(alpha, beta, varlambdaa1, varlambdaa2, k, ttotal, taccrual, gamma, delta))/(1-r))

  if(design == "parallel" & test == "equivalence") n_trt <- ceiling((TrialSize::TwoSampleSurvival.Equivalence(alpha, beta, varlambdaa1, varlambdaa2, k, ttotal, taccrual, gamma, delta))/(1-r))

  if(design == "crossover" & test %in% c("equality", "noninferiority", "superiority", "equivalence")){
    message("Check back next version!")
    n_trt <- NA
  }
  ans <- stats::setNames(data.frame(matrix(data = c(n_trt,k*n_trt), ncol=2, nrow = 1),row.names = "Size"), c("n_trt", "n_ctl"))

  return(ans)

}

#' @examples
#' # Ex 1. (n_trt=56, n_ctl=56)
#' getSizeTTE(design="parallel", test="equality", alpha=0.05, beta=0.20,
#'  varlambda=c(1,2), k=1, ttotal=3, taccrual=1, gamma=0.00001, delta=0,
#'  rho=c(0.05, 0.07), r=0.1)
#'
#' # Ex 2. (n_trt=30, n_ctl=30)
#' getSizeTTE(design="parallel", test="noninferiority", alpha=0.05, beta=0.20,
#'  varlambda=c(1,2), k=1, ttotal=3, taccrual=1, gamma=0.00001, delta= -0.2,
#'  rho=c(0.05, 0.07), r=0.1)
#'
#' # Ex 3. (n_trt=74, n_ctl=74)
#' getSizeTTE(design="parallel", test="superiority", alpha=0.05, beta=0.20,
#' varlambda=c(1,2), k=1, ttotal=3, taccrual=1, gamma=0.00001, delta=0.20,
#' rho=c(0.05, 0.07), r=0.1)
#'
#' # Ex 4. (n_trt=84, n_ctl=84)
#' getSizeTTE(design="parallel", test="equivalence", alpha=0.05, beta=0.20,
#'  varlambda=c(1,1), k=1, ttotal=3, taccrual=1, gamma=0.00001, delta=0.5,
#'  rho=c(0.05, 0.07), r=0.1)
#'
#' # Ex 5. (Check back next version)
#' getSizeTTE(design="crossover", test="equality", alpha=0.05, beta=0.20,
#'  varlambda=c(1,1), k=1, ttotal=3, taccrual=1, gamma=0.00001, delta=0.5,
#'  rho=c(0.05, 0.07), r=0.1)
#'
#' # Ex 6. (Check back next version)
#' getSizeTTE(design="crossover", test="noninferiority", alpha=0.05,
#'  beta=0.20, varlambda=c(1,1), k=1, ttotal=3, taccrual=1, gamma=0.00001,
#'  delta=0.5, rho=c(0.05, 0.07), r=0.1)
#'
#' # Ex 7. (Check back next version)
#' getSizeTTE(design="crossover", test="superiority", alpha=0.05, beta=0.20,
#'  varlambda=c(1,1), k=1, ttotal=3, taccrual=1, gamma=0.00001, delta=0.5,
#'  rho=c(0.05, 0.07), r=0.1)
#'
#' # Ex 8. (Check back next version)
#' getSizeTTE(design="crossover", test="equivalence", alpha=0.05, beta=0.20,
#'  varlambda=c(1,1), k=1, ttotal=3, taccrual=1, gamma=0.00001, delta=0.5,
#'  rho=c(0.05, 0.07), r=0.1)


#' @rdname getSize
#' @param varcatprob list of two probability vectors per treatment arm
#(varcatprob1>0,varcatprob2>0):=(vector of category probs in control arm,vector of category probs in control arm)
#' @param theta log odds ratio of outcome in treatment arm versus control arm
#' @export
getSizeOrd <- function(design=c("parallel", "crossover"),
                       test=c("equality", "noninferiority", "superiority", "equivalence"),
                       alpha=0.05, beta=0.20, varcatprob, k=1,
                       theta, delta=0, rho=c(0.05, 0.07), r=0.1){

  stopifnot(design %in% c("parallel", "crossover"),
            test %in% c("equality", "noninferiority", "superiority", "equivalence"),
            dplyr::between(alpha,0,1), dplyr::between(beta,0,1),
            dplyr::between(min(varcatprob[[1]]),0,Inf),
            dplyr::between(max(varcatprob[[1]]),-Inf,1),
            dplyr::between(min(varcatprob[[2]]),0,Inf),
            dplyr::between(max(varcatprob[[2]]),-Inf,1),
            dplyr::between(k,0,Inf), dplyr::between(theta,-Inf,Inf),
            dplyr::between(delta,-Inf,Inf), dplyr::between(rho[1],0,1),
            dplyr::between(rho[2],0,1), dplyr::between(r,0,1))

  varcatprobb1 <- rho[1]*varcatprob[[2]]+(1-rho[1])*varcatprob[[1]]
  varcatprobb2 <- rho[2]*varcatprob[[1]]+(1-rho[2])*varcatprob[[2]]
  varcatprobave <- (varcatprobb1+varcatprobb2)/2

  if(design == "parallel" & test == "equality") n_total <- ceiling((Hmisc::posamsize(varcatprobave,exp((1-rho[1]-rho[2])*theta),(1/(1+k)),alpha,1-beta))$n/(1-r))

  if(design == "parallel" & test == "noninferiority") n_total <- ceiling((Hmisc::posamsize(varcatprobave,exp((1-rho[1]-rho[2])*theta-delta),(1/(1+k)),2*alpha,1-beta))$n/(1-r))

  if(design == "parallel" & test == "superiority" ) n_total <- ceiling((Hmisc::posamsize(varcatprobave,exp((1-rho[1]-rho[2])*theta-delta),(1/(1+k)),2*alpha,1-beta))$n/(1-r))

  if(design == "parallel" & test == "equivalence") n_total <- ceiling((Hmisc::posamsize(varcatprobave,exp(abs((1-rho[1]-rho[2])*theta)-delta),(1/(1+k)),2*alpha,1-(beta/2)))$n/(1-r))

  if(design == "crossover" & test %in% c("equality", "noninferiority", "superiority", "equivalence")){
    message("Check back next version!")
    n_total <- NA
  }

  ans <- stats::setNames(
    data.frame(matrix(data = c(ceiling((n_total/(k+1))), ceiling((k*n_total/(k+1)))), ncol=2, nrow=1), row.names="Size"),
    c("n_trt", "n_ctl"))
  return(ans)

}

#' @examples
#' # Ex 1. (n_trt=135, n_ctl=135)
#' getSizeOrd(design="parallel", test="equality", alpha=0.05, beta=0.10,
#'  varcatprob = list(c(0.2,0.5,0.2,0.1), c(0.378,0.472,0.106,0.044)),
#'  k=1, theta=0.887, delta=0, rho=c(0.05, 0.07), r=0.1)
#'
#' # Ex 2. (Check back next version)
#' getSizeOrd(design="crossover", test="equality", alpha=0.05, beta=0.10,
#'  varcatprob = list(c(0.2,0.5,0.2,0.1), c(0.378,0.472,0.106,0.044)),
#'  k=1, theta=0.887, delta=0, rho=c(0.05, 0.07), r=0.1)
#'


