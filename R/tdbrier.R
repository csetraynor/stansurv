#'  Calculate tdBrier
#'
#' These functions calculate the survival analysis metric measured of a
#' system compared to a hold-out test set. The measurement and the "truth"
#' have a survival time and a censoring indicator 0/1 indicating if the event
#' result or the event.
#'
#'
#' The Brier score is defined as the squared distance between the
#' expected survival probability and the observed survival.
#' Therefore, it measures the discrepancy between observation
#' and model-based prediction.
#'
#' The integrated Brier Score summarises the Brier Score over the range
#' of observed events.Similar to the original Brier score [40] the iBrier:
#' ranges from 0 to 1; the model with an out-of-training sample value closer
#' to 0 outperforms the rest.
#' @aliases tdbrier tdbrier.model.list tdbrier.int.matrix
#' stdbrier.int.reference
#' @param data For the default functions, a datframe containing survival
#' (time), and status (0:censored/1:event), and the explanatory variables.
#' @param mod Coxph model object fitted with coxph (survival).
#' @return A tdBrier object
#' @seealso [iBrier]
#' @keywords brier
#' @examples
#' require(survival)
#' require(dplyr)
#' data(lung)
#' lung <- lung %>%
#' mutate(status = (status == 2))
#'
#' mod <- coxph(Surv(time, status)~ age, data = lung)
#'
#' tdbrier <- get_tdbrier(lung, mod)
#' integrate_tdbrier(tdroc)
#'
#' @export tdbrier
#' @references
#'
#'  Ulla B. Mogensen, Hemant Ishwaran, Thomas A. Gerds (2012).
#' Evaluating Random Forests for Survival Analysis Using Prediction Error
#' Curves. Journal of Statistical Software, 50(11), 1-23.
#' URL http://www.jstatsoft.org/v50/i11/.
#' @export tdbrier

tdbrier <- function(holdout, fit, time = "time", status = "status", ncores = 1L){
  #select timepoints with observed event
  test.time <- get_obs.time(holdout)
  timepoints <- obs.times <- seq(0, max(test.time), length.out = 100)

  if(any( class(fit)  == "coxph") ){
    #get probabilites
    probs <- pec::predictSurvProb(fit, newdata = holdout, times = timepoints)
    cox.brier <- suppressMessages( pec::pec(probs, Surv(time, status) ~ 1,
                                            data = holdout,
                                            maxtime = max(timepoints),
                                            exact = FALSE,
                                            exactness = length(timepoints) - 1 ) )
    apperror <- cox.brier$AppErr$matrix
    time <- timepoints
    ref <- cox.brier$AppErr$Reference
    ibrier <- integrate_tdbrier(cox.brier)
    object <- cox.brier

    rsquaredbs <- suppressWarnings ( pec::R2(object = cox.brier)$AppErr )
    out <- list(object, rsquaredbs, ibrier, apperror, ref, time)
    names(out) <- c("brier.object", "rsquaredbs","ibrier", "apperror", "reference","time")
    out

  } else {

    if(any( class(fit$stanfit) == "stanfit" ) ){

      surv.prob <- pred_surv(fit = fit,
                             timepoints = timepoints)

      long_test <-  gen_new_frame(dat = test, timepoints = timepoints)

      pred_frame <- pred_surv(long_x = long_test,
                              fit = fit,
                              ncores = ncores)


      mean.surv <- apply(pred_frame, 1, mean)
      long_test$surv <- mean.surv

      # test <- long_test %>%
      #   dplyr::select(time, surv, patient_id) %>%
      #   dplyr::group_by(patient_id) %>%
      #   dplyr::slice( n()) %>%
      #   dplyr::ungroup()

      probs <- get_survProb(newdat = long_test)

      brier <- suppressMessages( pec::pec(probs, Surv(time, status) ~ 1,
                                          data = test,
                                          maxtime = max(timepoints),
                                          exact = FALSE,
                                          exactness = length(ncol(probs) ) -1  ) )


      brier
      ibrier <- integrate_tdbrier(brier)
      brier$AppErr$matrix

      models.brier <- parallel::mclapply(seq_along(colnames(pred_frame)), function(i){
        newdat$surv <- as.vector( pred_frame[ ,i] )

        probs <- get_survProb(newdat = newdat)

        brier <- suppressMessages( pec::pec(probs, Surv(time, status) ~ 1,
                                            data = test,
                                            maxtime = max(timepoints),
                                            exact = FALSE,
                                            exactness = 99L) )
        brier
      }, mc.cores = ncores )

      models.brier.error <- lapply(models.brier, function(b){
        b$AppErr$matrix
      })
      models.ibrier <- sapply(models.brier, function(b){
        integrate_tdbrier(b)
      })
      models.brier.error <- do.call(cbind, models.brier.error )
      apperror <- apply(models.brier.error, 1, mean)
      time <- timepoints
      ibrier <- mean(models.ibrier)

      rsquaredbs <- lapply(models.brier, function(b){
        unlist( suppressWarnings( pec::R2(object = b)$AppErr ) )
      })
      rsquaredbs <- do.call(cbind, rsquaredbs)
      rsquaredbs <- apply(rsquaredbs, 1, mean)
      rsquaredbs <- as.list(rsquaredbs)
      names(rsquaredbs) <- c("time", "model")
      out <- list(rsquaredbs, ibrier, apperror, time)
      names(out) <- c( "rsquaredbs","ibrier", "apperror","time")

      out

    } else {
      "Print object must be either coxph or stanreg"
    }
  }
}

#' @rdname tdbrier
#' @export
integrate_tdbrier <-
  function(x, ...) {
    stop <- max(x$time[!is.na(x$AppErr$matrix)])
    ibrier <- pec::crps(x, models = "matrix", times = stop)[1]
    ibrier <- unlist(ibrier)
    return(ibrier)
  }
#' @export
#' @rdname tdbrier
integrate_tdbrier_reference <-
  function(x, ...) {
    stop <- max(x$time[!is.na(x$AppErr$Reference)])
    ibrier <- pec::crps(x, models = "Reference", times = stop)[1]
    ibrier <- unlist(ibrier)
    return(ibrier)
  }


