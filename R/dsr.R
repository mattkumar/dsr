#' Compute Directly Standardized Rates
#'
#' Computes crude and directly standardized rates by subgroup with confidence intervals.
#'
#' @param data A data frame with counts and unit-times summarized by the standardization variables.
#' @param event A variable within the input data that corresponds to the event counts.
#' @param fu A variable within the input data that corresponds to the unit-time.
#' @param subgroup A variable within the input data frame for which rates are calculated by.
#' @param ... Variables(s) within the input data that for which rates are to be standardized by. The input data and ref data should both be summarized by these.
#' @param refdata A data frame with population unit-times summarized by the standardization variables. The unit-time variable name must named pop.
#' @param mp A constant to multiply rates by (e.g. mp=1000 for rates per 1000).
#' @param method Choose between normal, lognormal and gamma confidence intervals for crude and standardized rates. The default method is normal.
#' @param sig The desired level of confidence in computing confidence intervals. The default is 0.95 for 95 percent CIs.
#' @param decimals Round estimates to a desired decimal place.
#' @references Fay, M.P., & Feuer, E.J. (1997). Confidence intervals for directly standardized rates: a method based on the gamma distribution. Statistics in Medicine,16, 791-801.
#' @references Elandt-Johnson, R. C., and Johnson, N. L. (1980). Survival Models and Data Analysis. New York: John Wiley & Sons.
#' @references Chiang C. Standard error of the age-adjusted death rate. US Department of Health, Education and Welfare: Vital Statistics Special Reports 1961;47:271-285.
#' @references Schoenbach, V., and Rosamond W. (2000) Understanding the fundamentals of epidemiology: An evolving text.
#' @import stats
#' @import dplyr
#' @import rlang
#' @import utils
#' @examples
#' #An example of calculating directly standardized rates
#' #Data from Table 1, Page 132 of Schoenbach (2000)
#'
#' #State specific death counts and fu
#' df_study <- data.frame(state=rep(c('Miami',"Alaska"), c(5,5)),
#'                       age=rep(c('00-14','15-24','25-44','45-64','65+'),2),
#'                       deaths=c(136,57,208,1016,3605,59,18,37,90,81),
#'                       fu=c(114350,80259,133440,142670,92168,37164,20036,32693,14947,2077))
#'
#' #US standard population
#' df_ref  <- data.frame(age=c('00-14','15-24','25-44','45-64','65+'),
#'                      pop=c(23961000,15420000,21353000,19601000,10685000))
#'
#' #Directly Standardized Rates (per 1000) - 95% CI's using the gamma method
#' my_results <- dsr(data=df_study,
#'                  event=deaths,
#'                  fu=fu,
#'                  subgroup=state,
#'                  age,
#'                  refdata=df_ref,
#'                  method="gamma",
#'                  sig=0.95,
#'                  mp=1000,
#'                  decimals=4)
#' #View results
#' my_results
#' @export

dsr <- function(data, event, fu, subgroup, ..., refdata, mp, method="normal", sig=0.95, decimals ) {

  subgroup <- enquo(subgroup)
  event <- enquo(event)
  fu <- enquo(fu)
  stdvars <- quos(...)



  #Compute crude and standardized rates and variances
  all_data_st = data %>%
    left_join(refdata) %>%
    group_by(!!subgroup) %>%
    mutate(n=sum(!!event),
           d=sum(!!fu),
           cr_rate=n/d,
           cr_var=n/d^2,
           wts=pop/sum(pop),
           st_rate=sum(wts*(!!event/!!fu)),
           st_var=sum(as.numeric((wts^2)*(!!event/(!!fu )^2)))) %>%
    distinct(!!subgroup, .keep_all=TRUE) %>%
    select(!!subgroup, n, d, cr_rate, cr_var, st_rate, st_var)



  #Compute Confidence Intervals (CI) according to method. The default is 'normal'
  if(method=="gamma") {

    tmp1 =   all_data_st %>%
      mutate(
        c_rate=mp*cr_rate,
        c_lower=mp*qgamma((1-sig)/2, shape=cr_rate^2/(cr_var))/(cr_rate/cr_var),
        c_upper=mp*qgamma(1-((1-sig)/2), shape=1+cr_rate^2/(cr_var))/(cr_rate/cr_var),
        s_rate=mp*st_rate,
        s_lower=mp*qgamma((1-sig)/2, shape=st_rate^2/st_var)/(st_rate/st_var),
        s_upper=mp*qgamma(1-((1-sig)/2), shape=1+(st_rate^2/st_var))/(st_rate/st_var)) %>%
      select(!!subgroup, n, d, c_rate, c_lower, c_upper, s_rate, s_lower, s_upper)


  } else if (method=="normal") {


    tmp1 =   all_data_st %>%
      mutate(
        c_rate=mp*cr_rate,
        c_lower=mp*(cr_rate+qnorm((1-sig)/2)*sqrt(cr_var)),
        c_upper=mp*(cr_rate-qnorm((1-sig)/2)*sqrt(cr_var)),
        s_rate=mp*st_rate,
        s_lower=mp*(st_rate+qnorm((1-sig)/2)*sqrt(st_var)),
        s_upper=mp*(st_rate-qnorm((1-sig)/2)*sqrt(st_var))) %>%
      select(!!subgroup, n, d, c_rate, c_lower, c_upper, s_rate, s_lower, s_upper)

  } else if (method=="lognormal") {


    tmp1 =   all_data_st %>%
      mutate(
        c_rate=mp*cr_rate,
        c_lower=mp*exp((log(cr_rate)+qnorm((1-sig)/2)*sqrt(cr_var)/(cr_rate))),
        c_upper=mp*exp((log(cr_rate)-qnorm((1-sig)/2)*sqrt(cr_var)/(cr_rate))),
        s_rate=mp*st_rate,
        s_lower=mp*exp((log(st_rate)+qnorm((1-sig)/2)*sqrt(st_var)/(st_rate))),
        s_upper=mp*exp((log(st_rate)-qnorm((1-sig)/2)*sqrt(st_var)/(st_rate)))) %>%
      select(!!subgroup, n, d, c_rate, c_lower, c_upper, s_rate, s_lower, s_upper)

  }


  #Clean up and output

  tmp1$c_rate  <- round(tmp1$c_rate,  digits=decimals)
  tmp1$c_lower <- round(tmp1$c_lower, digits=decimals)
  tmp1$c_upper <- round(tmp1$c_upper, digits=decimals)
  tmp1$s_rate  <- round(tmp1$s_rate, digits=decimals)
  tmp1$s_lower <- round(tmp1$s_lower, digits=decimals)
  tmp1$s_upper <- round(tmp1$s_upper, digits=decimals)

  colnames(tmp1) <-  c('Subgroup', 'Numerator','Denominator',
                       paste('Crude Rate (per ',mp,')',sep=''),
                       paste(sig*100,'% LCL (Crude)',sep=''),
                       paste(sig*100,'% UCL (Crude)',sep=''),
                       paste('Std Rate (per ',mp,')',sep=''),
                       paste(sig*100,'% LCL (Std)',sep=''),
                       paste(sig*100,'% UCL (Std)',sep=''))


  tmp1 <- as.data.frame(tmp1)

}
