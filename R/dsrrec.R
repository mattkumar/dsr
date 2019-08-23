#' Compute Directly Standardized Rates for Recurrent Events
#'
#' Computes directly standardized rates for recurrent events by subgroup with confidence intervals.
#'
#' @param data A data frame with counts and unit-times summarized by the standardization variables.
#' @param event A variable within the input data that corresponds to the event counts.
#' @param fu A variable within the input data that corresponds to the unit-time.
#' @param subgroup A variable within the input data frame for which rates are calculated by.
#' @param ... Variables(s) within the input data that for which rates are to be standardized by. The input data and ref data should both be summarized by these.
#' @param refdata A data frame with population unit-times summarized by the standardization variables. The unit-time variable name must named pop.
#' @param mp A constant to multiply rates by (e.g. mp=1000 for rates per 1000).
#' @param sig The desired level of confidence in computing confidence intervals. The default is 0.95 for 95 percent CIs.
#' @param decimals Round estimates to a desired decimal place.
#' @references Stukel, T. A., Glynn, R. J., Fisher, E. S., Sharp, S. M., Lu-Yao, G and Wennberg, J. E. (1994). Standardized rates of recurrent outcomes. Statistics in Medicine, 13, 1781-1791.
#' @references Fay, M.P., & Feuer, E.J. (1997). Confidence intervals for directly standardized rates: a method based on the gamma distribution. Statistics in Medicine, 16, 791-801.
#' @import stats
#' @import dplyr
#' @import rlang
#' @import utils
#' @examples
#' #An example of directly standardized rates for recurrent events
#'
#'library(frailtypack)
#'library(dplyr)
#'library(dsr)
#'data(readmission)
#'
#'#Make an individual level dataset with total event counts and total observation times
#'treadm <- as.data.frame(readmission %>%
#'                          group_by(id) %>%
#'                          filter(max(enum)==enum ) %>%
#'                          mutate(events=enum-1, time=t.stop) %>%
#'                          select(id, events, time, sex, dukes))
#'
#'#Make the standard pop
#'tref <- as.data.frame(treadm %>%
#'                      group_by(sex) %>%
#'                      mutate(pop=sum(time)) %>%
#'                      select(sex, pop) %>%
#'                      distinct(sex, pop))
#'
#'#Get directly standardized rates (age-adjusted) for readmissions by Dukes' tumor grade.
#'analysis <- dsrrec(data=treadm,
#'                   event=events,
#'                   fu=time,
#'                   refdata=tref,
#'                   subgroup=dukes,
#'                   sex,
#'                   mp=1000,
#'                   decimals=3)
#' @export

dsrrec <- function(data,
                    event,
                    fu,
                    subgroup,
                    ...,
                    refdata,
                    sig=0.95,
                    mp,
                    decimals
                    )
  {

  subgroup  <- enquo(subgroup)
  event <- enquo(event)
  fu <- enquo(fu)
  stdvars  <- quos(...)


  #Compute element: k
  #calculate elements

  #k parameter
  tk <- as.data.frame(  data %>%
                        group_by(!!!stdvars,!!subgroup) %>%
                        mutate(n=sum(!!event), followup=sum(!!fu), k=((mean(!!event)^2) / (var(!!event)-mean(!!event)))) %>%
                        select(!!!stdvars,!!subgroup, k, n, followup) %>%
                        distinct(!!!stdvars,!!subgroup, k, n, followup)
                      )

  #pop-ref
  divisor <- sum(refdata$pop)
  ttmp <- tk %>% left_join(refdata) %>% mutate(total=sum(pop), wts=pop/divisor)


  #get variances
  ttest <- ttmp %>% group_by(!!!stdvars,!!subgroup) %>% mutate(
                                                                std_rate=((n/followup) * wts),
                                                                nb_var = as.numeric((wts)^2 * ((n/(followup)^2) + ((n)^2)/(k * (followup)^3))),
                                                                gam_var= as.numeric((wts^2)*(n/(followup)^2)))

  #calculate ci's, rates
  tmp1 <- ttest %>% group_by(!!subgroup) %>%
    mutate(s_rate=sum(std_rate),
           std_rate=mp*(sum(std_rate)),

           n=sum(n),
           d=sum(followup),
           cr_rate=(n/d)*mp,
           nb_lcl=mp*(s_rate + qnorm((1-sig)/2) * sqrt(sum(nb_var))),
           nb_ucl=mp*(s_rate - qnorm((1-sig)/2) * sqrt(sum(nb_var))),
           gam_lcl=mp*(qgamma((1-sig)/2, shape=s_rate^2/sum(gam_var))/(s_rate/sum(gam_var))),
           gam_ucl=mp*(qgamma(1-((1-sig)/2), shape=1+(s_rate^2/sum(gam_var)))/(s_rate/sum(gam_var))))  %>%
    select(!!subgroup, n, d, cr_rate, std_rate, nb_lcl, nb_ucl, gam_lcl, gam_ucl) %>%
    distinct(!!subgroup, n, d, cr_rate, std_rate, nb_lcl, nb_ucl, gam_lcl, gam_ucl)


  tmp1$cr_rate  <- round(tmp1$cr_rate,  digits=decimals)
  tmp1$std_rate <- round(tmp1$std_rate, digits=decimals)
  tmp1$nb_lcl   <- round(tmp1$nb_lcl,   digits=decimals)
  tmp1$nb_ucl   <- round(tmp1$nb_ucl,   digits=decimals)
  tmp1$gam_lcl  <- round(tmp1$gam_lcl,  digits=decimals)
  tmp1$gam_ucl  <- round(tmp1$gam_ucl,  digits=decimals)

  colnames(tmp1) <-  c('Subgroup', 'Numerator','Denominator',
                       paste('Cr Rate (per ',mp,')',sep=''),
                       paste('Std Rate (per ',mp,')',sep=''),
                       paste(sig*100,'% LCL (NB)',sep=''),
                       paste(sig*100,'% UCL (NB)',sep=''),
                       paste(sig*100,'% LCL (Gamma)',sep=''),
                       paste(sig*100,'% UCL (Gamma)',sep=''))


  tmp1 <- as.data.frame(tmp1)



}
