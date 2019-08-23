#' Compare Directly Standardized Rates by Ratios or Differences.
#'
#' Compare directly standardized rates by ratios or differences.
#' @param data A data frame with counts and unit-times summarized by the standardization variables.
#' @param event A variable within the input data that corresponds to the event counts.
#' @param fu A variable within the input data that corresponds to the unit-time.
#' @param subgroup A variable within the input data frame for which rates are calculated by.
#' @param ... Variables(s) within the input data that for which rates are to be standardized by. The input data and ref data should both be summarized by these.
#' @param refdata A data frame with population unit-times summarized by the standardization variables. The unit-time variable name must named pop.
#' @param estimate Choose between difference or ratio in comparing directly standardized rates.
#' @param refgroup A level of the subgroup variable taken to be the reference in computing rate ratios or differences.
#' @param mp A constant to multiply rates by (e.g. mp=1000 for rates per 1000).
#' @param sig The desired level of confidence in computing confidence intervals. The default is 0.95 for 95 percent CIs.
#' @param decimals Round estimates to a desired decimal place.
#' @references Fay, M.P., & Feuer, E.J. (1997). Confidence intervals for directly standardized rates: a method based on the gamma distribution. Statistics in Medicine,16, 791-801.
#' @references Elandt-Johnson, R. C., and Johnson, N. L. (1980). Survival Models and Data Analysis. New York: John Wiley & Sons.
#' @references Chiang C. Standard error of the age-adjusted death rate. US Department of Health, Education and Welfare: Vital Statistics Special Reports 1961;47:271-285.
#' @references Schoenbach, V., and Rosamond W. (2000) Understanding the fundamentals of epidemiology: An evolving text.
#' @import stats
#' @import stats
#' @import dplyr
#' @import rlang
#' @import utils
#' @examples
#' #An example of comparing directly standardized rates
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
#' #Directly Standardized Rate Ratio (per 1000) - 95% log-normal CI's, Alaska as the refernce
#' my_results2 <- dsrr(data=df_study,
#'                    event=deaths,
#'                    fu=fu,
#'                    subgroup=state,
#'                    age,
#'                    refdata=df_ref,
#'                    refgroup="Alaska",
#'                    estimate="ratio",
#'                    sig=0.95,
#'                    mp=1000,
#'                    decimals=4)
#' #View results
#' my_results2
#' @export

dsrr <- function(data, event, fu, subgroup, ..., refdata, estimate, refgroup, mp, sig=0.95, decimals) {

  refgroup <- enquo(refgroup)
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
    select(!!subgroup, n, d, st_rate, st_var)




  #Select a Referent
  ref = all_data_st %>% filter(!! subgroup == !!  quo_name(refgroup)) %>% select(reference=!! subgroup, st_rate, st_var)

  #Compute comparisons based on desired estimate
  if(estimate=="ratio") {

    out = all_data_st %>%
      mutate(ratio=(st_rate/ref$st_rate),
             se=sqrt(((1/(st_rate^2)) * st_var ) + ((1/(ref$st_rate^2)) * ref$st_var)),
             lower = exp(log(ratio + qnorm((1-sig)/2) * se)),
             upper = exp(log(ratio - qnorm((1-sig)/2) * se)),
             reference = ref$reference,
             s_rate=mp*st_rate
      ) %>%
      select(!!subgroup, reference, s_rate, ratio, lower, upper)

    #Clean up and output


    out$s_rate  <- round(out$s_rate, digits=decimals)
    out$lower <- round(out$lower, digits=decimals)
    out$upper <- round(out$upper, digits=decimals)
    out$ratio <- round(out$ratio, digits=decimals)

    colnames(out) <-  c('Comparator', 'Reference',
                        paste('Std Rate (per ',mp,')',sep=''),
                        paste('Rate Ratio (RR)',sep=''),
                        paste(sig*100,'% LCL (RR)',sep=''),
                        paste(sig*100,'% UCL (RR)',sep='')
    )

    out <- as.data.frame(out)


  } else if (estimate=="difference") {

    out = all_data_st %>%
      mutate(ratio=(st_rate-ref$st_rate),
             se=sqrt((st_var) + (ref$st_var)),
             lower = mp*(ratio + qnorm((1-sig)/2) * se),
             upper = mp*(ratio - qnorm((1-sig)/2) * se),
             reference = ref$reference,
             s_rate=mp*st_rate,
             difference=mp*ratio
      ) %>%
      select(!!subgroup, reference, s_rate, difference, lower, upper)

    #Clean up and output

    out$s_rate  <- round(out$s_rate, digits=decimals)
    out$lower <- round(out$lower, digits=decimals)
    out$upper <- round(out$upper, digits=decimals)
    out$difference <- round(out$difference, digits=decimals)

    colnames(out) <-  c('Comparator', 'Reference',
                        paste('Std Rate (per ',mp,')',sep=''),
                        paste('Rate Difference (RD)',sep=''),
                        paste(sig*100,'% LCL (RD)',sep=''),
                        paste(sig*100,'% UCL (RD)',sep='')
    )

    out <- as.data.frame(out)




  }

}
