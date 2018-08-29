# Function to predict the survival for each individual
#'
#' Predict survival from a stan_surv object
#'
#' @export link_surv
#' @importFrom Hmisc approxExtrap
#'
#' @examples
#' pbc2 <- survival::pbc
#' pbc2 <- pbc2[!is.na(pbc2$trt),]
#' pbc2$status <- as.integer(pbc2$status > 0)
#' m1 <- stan_surv(survival::Surv(time, status) ~ trt, data = pbc2)
#'
#' df <- flexsurv::bc
#' m2 <- stan_surv(survival::Surv(rectime, censrec) ~ group,
#'                 data = df, cores = 1, chains = 1, iter = 2000,
#'                 basehaz = "fpm", iknots = c(6.594869,  7.285963 ),
#'                 degree = 2, prior_aux = normal(0, 2, autoscale = F))
#'
pred_surv2 <- function(fit, testx = NULL, timepoints = NULL, time = "time", status = "status", ncores = 1L, method = "gp"){
  obs.time <- get_obs.time(fit$data)
  if(is.null(timepoints)){
    timepoints <- test.time <- get_obs.time(testx)
  }
  basehaz.samples <- extract_basehaz_draws(fit)
  spline.basis <- extract_splines_bases(fit)

  if(!check_null_model(fit)){
    beta_draws <- extract_beta_draws(fit);
    varis.obs <- get_predictors(fit);

    basehaz.post <- lapply(seq_along(1:nrow(basehaz.samples)), function(s){
      sapply(seq_along(1:nrow(spline.basis)), function(n){
      basehaz.samples[s, ]  %*% spline.basis[n, ]
    })})

    basehaz.post <- do.call(cbind, basehaz.post)
    basehaz.post <- apply(X = basehaz.post, MARGIN = 1, FUN = mean)
    # obs.dat <- data.frame(haz = basehaz,
                           # time = obs.time )
    basehaz.post <- unlist(basehaz.post)
    if(method == "linear" | method == "constant"){
      # surv.approx <- approx(x = obs.time, y = surv.base, xout = timepoints, method = method)$y # only interpolation
      haz.approx <- Hmisc::approxExtrap(x = obs.time, y = basehaz.post, xout = timepoints, method = method)$y # with extrapolation
    } else if(method == "gp" | method == "gaussian process"){
      haz.bgp <- tgp:: bgp(X= obs.time, Z= basehaz.post, verb = 0)
      X <- data.frame(x1 = obs.time)
      XX <- data.frame(x1 = timepoints)
      haz.bgp$Xsplit <- rbind(X, XX)
      haz.approx <- predict(object = haz.bgp, XX = timepoints)$ZZ.mean
    } else if(method == "tgp" | method == "treed gaussian process"){
      haz.btgp <- tgp:: btgp(X= obs.time,
                             Z= basehaz.post, verb = 0)
      X <- data.frame(x1 = obs.time)
      XX <- data.frame(x1 = timepoints)
      haz.btgp$Xsplit <- rbind(X, XX)
      haz.approx <- predict(object = haz.btgp, XX = timepoints)$ZZ.mean
    } else {
      stop2(paste0("method ", method, " not implemented"))
    }
    surv.approx <- link_surv_base.helper(b = haz.approx)
    surv.approx[surv.approx > 1] <- 1
    surv.approx[surv.approx < 0] <- 0

ind.post <- lapply(seq_along(1:nrow(testx)),
                               function(i){
  surv.post <- sapply(seq_along(1:nrow(beta_draws)), function(s){
    lin_post <- varis.obs[i, ] %*% beta_draws[s, ]
    link_surv_lin.helper(p = lin_post, s = surv.approx)
  })
  apply(surv.post, 1, mean)
})

ind.post <- do.call(rbind, ind.post)
rownames(ind.post) <- rownames(testx)
  } else {
    # For null model
    surv.post <- lapply(seq_along(1:nrow(basehaz.samples)), function(s){
      basehaz.post <- sapply(seq_along(1:nrow(spline.basis)), function(n){
        basehaz.samples[s, ]  %*% spline.basis[n, ]
      })
      surv_base <- link_surv_base.helper(b = basehaz.post)
      surv.approx <- Hmisc::approxExtrap(x = obs.time, y = surv_base, xout = timepoints, method = method)$y
      # surv.approx <- approx(x = obs.time, y = surv_base, xout = timepoints, method = method)$y
      surv.approx[surv.approx > 1] <- 1
      surv.approx[surv.approx < 0] <- 0
      surv.approx
    })
    surv.post <- do.call(cbind, surv.post)
    surv.post <- apply(surv.post, 1, median)
    #No difference in survival expected in the null model
    ind.post <- matrix(rep(surv.post, nrow(testx)),
                       nrow = nrow(testx),
                       ncol = n_distinct(timepoints),
                       byrow = TRUE)
  }
  return(ind.post)
}

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
#' @export tdbrier
#' @references
#'
#'  Ulla B. Mogensen, Hemant Ishwaran, Thomas A. Gerds (2012).
#' Evaluating Random Forests for Survival Analysis Using Prediction Error
#' Curves. Journal of Statistical Software, 50(11), 1-23.
#' URL http://www.jstatsoft.org/v50/i11/.
#' @export tdbrier

tdbrier <- function(holdout, fit, time = "time", status = "status", ncores = 1L, method = "gp"){
  #select timepoints with observed event
  holdout <- prepare_surv(holdout, time = time, status = status)
  timepoints =  seq(min(holdout[[time]]), max(holdout[[time]]), length.out = 100)

  if(any( class(fit)  == "coxph") ){
    #get probabilites
    probs <- pec::predictSurvProb(fit, newdata = holdout, times = timepoints)
    cox.brier <- suppressMessages( pec::pec(probs, Surv(time, status) ~ 1,
                                            data = holdout,
                                            maxtime = max(timepoints),
                                            exact = FALSE,
                                            exactness = 99L ) )
    apperror <- cox.brier$AppErr$matrix
    time <- timepoints
    ref <- cox.brier$AppErr$Reference
    ibrier <- integrate_tdbrier(cox.brier)
    ibrier_ref <- integrate_tdbrier_reference(cox.brier)
    object <- cox.brier

    rsquaredbs <- suppressWarnings ( pec::R2(object = cox.brier, times = cox.brier$maxtime, reference = 1 ) )
    out <- list(object, rsquaredbs, ibrier, ibrier_ref, apperror, ref, time)
    names(out) <- c("brier.object", "rsquaredbs","ibrier", "ibrieref", "apperror", "reference","time")
    out

  } else if (any( class(fit$stanfit) == "stanfit" ) ){
    ind.pred <- pred_surv2(fit = fit, testx = holdout, timepoints = timepoints, ncores = ncores, method = method)

    ba.brier <- suppressMessages( pec::pec(ind.pred, Surv(time, status) ~ 1,
                                           data = holdout,
                                           maxtime = max(timepoints),
                                           exact = FALSE,
                                           exactness = 99L ) )
    apperror <- ba.brier$AppErr$matrix
    time <- timepoints
    ref <- ba.brier$AppErr$Reference
    ibrier <- integrate_tdbrier(ba.brier)
    ibrier_ref <- integrate_tdbrier_reference(ba.brier)
    object <- ba.brier

    rsquaredbs <- suppressWarnings ( pec::R2(object = ba.brier, times = ba.brier$maxtime, reference = 1) )
    out <- list(object, rsquaredbs, ibrier, ibrier_ref, apperror, ref, time)
    names(out) <- c("brier.object", "rsquaredbs","ibrier", "ibrier_ref", "apperror", "reference","time")
    out
  } else {
    stop2(paste0("Error: No method for evaluating predicted probabilities from objects in class", class(fit)))
  }
}

#' @rdname tdbrier
#' @export
integrate_tdbrier <-
  function(x, ...) {
    ibrier <- pec::crps(x, models = "matrix", times = x$maxtime)
    return(ibrier)
  }
#' @export
#' @rdname tdbrier
integrate_tdbrier_reference <-
  function(x, ...) {
    ibrier <- pec::crps(x, models = "Reference", times =  x$maxtime)
    return(ibrier)
  }


#' Get c-index
#'
#'
#' The efficacy of the survival model can be measured by
#'  the concordance statistic
#'
#' @param data For the default functions, a datframe containing survival
#' (time), and status (0:censored/1:event), and the explanatory variables.
#' @param mod Coxph model object fitted with coxph (survival).
#' @return A cindex object
#' @seealso [coxph]
#' @keywords cindex
#' @examples
#' require(survival)
#' require(dplyr)
#' data(lung)
#' lung <- lung %>%
#' mutate(status = (status == 2))
#'
#' mod <- coxph(Surv(time, status)~ age, data = lung)
#'
#' get_cindex(lung, mod)
#'
#' @export get_cindex
#' @author Carlos S Traynor
#' @references
#'
#'  Terry M. Therneau and Patricia M. Grambsch (2000).
#'   _Modeling Survival Data: Extending the Cox Model_.
#'   Springer, New York. ISBN 0-387-98784-3.
#'   @export get_cindex
#'
get_cindex <- function(data, mod,...)
  UseMethod("get_cindex")

#' @export
#' @rdname get_cindex
get_cindex <-
  function(holdout, fit, time = "time", status = "status", ncores = 1L, method = "gp"){

    #select timepoints with observed event
    holdout <- prepare_surv(holdout, time = time, status = status)
    timepoints =  seq(min(get_obs.time(holdout)), max(get_obs.time(holdout)), length.out = 100)

    if(any( class(fit)  == "coxph") ){
      #get probabilites
      probs <- pec::predictSurvProb(fit, newdata = holdout, times = timepoints)

      statistic <- suppressMessages( pec::cindex(probs,
                            Surv(time, status) ~ 1,
                                 data = holdout,
                           pred.times = timepoints,
                      eval.times = timepoints ) )

      apperror <- unlist(statistic$AppCindex)
      time <- timepoints
      concordance <-  mean(apperror)
      ref <- 0.5
      out <- list(concordance, apperror, ref, time, statistic)
      names(out) <- c( "concordance", "apperror", "reference","time",
                       "object")

      out

    } else if (any( class(fit$stanfit) == "stanfit" ) ){
      ind.pred <- pred_surv2(fit = fit, testx = holdout, timepoints = timepoints, ncores = ncores, method = method)
      statistic <- suppressMessages( pec::cindex(ind.pred,
                                 Surv(time, status) ~ 1,
                                data = holdout,
                             pred.times = timepoints,
                         eval.times = timepoints ) )

      apperror <- unlist(statistic$AppCindex)
      time <- timepoints
      concordance <-  mean(apperror)
      ref <- 0.5
      out <- list(concordance, apperror, ref, time, statistic)
      names(out) <- c( "concordance", "apperror", "reference","time", "object")

      out

    } else {
      stop2(paste0("Error: No method for evaluating predicted probabilities from objects in class", class(fit)))
    }

  }


#'  Calculate survival ROC
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
#' @export get_survroc
#' @author Carlos S Traynor
#' @references
#'
#'  Ulla B. Mogensen, Hemant Ishwaran, Thomas A. Gerds (2012).
#' Evaluating Random Forests for Survival Analysis Using Prediction Error
#' Curves. Journal of Statistical Software, 50(11), 1-23.
#' URL http://www.jstatsoft.org/v50/i11/.
get_survroc <- function(test, fit, time = "time", status = "status", ncores = ncores){

  #select timepoints with observed event
  holdout <- prepare_surv(holdout, time = time, status = status)
  timepoints =  seq(min(get_obs.time(holdout)), max(get_obs.time(holdout)), length.out = 100)

  if(any( class(fit)  == "coxph") ){
    probs <- predict(fit, newdata = holdout, type = "lp")

    statistic <- tdROC::tdROC(X = probs[!is.na(probs)],
                              Y = holdout[[time]][!is.na(probs)],
                              delta = holdout[[status]][!is.na(probs)],
                              tau = max(holdout[[time]][holdout[[status]] ]),
                              n.grid = 1000)

    iroc <- integrate_tdroc(statistic)

    sens <-  statistic$ROC$sens

    spec <-  statistic$ROC$spec

    grid <- statistic[[1]]$ROC$grid

    out <- list(iroc, sens, spec, grid)
    names(out) <- c( "iroc", "sens", "spec","grid")

    out



  } else if (any( class(fit$stanfit) == "stanfit" ) ){

    if(!check_null_model(fit)){
      beta_draws <- extract_beta_draws(fit);
      varis.obs <- get_predictors(fit);

      x <- holdout[fit$formula$allvars[-(1:2)]]

      X <- data.frame(id = seq_along(1:nrow(x)))

      for(i in seq_along(colnames(x)) ){

        if(is.numeric(x[[i]])){
          matrixna <- x[i]
          X <- cbind(X, matrixna)
        } else{
          x[[i]] <- as.factor(x[[i]])
          coluna <- x[i]
          matrixna <- model.matrix(~ ., coluna)[ ,-1]
          matrixna <- as.data.frame(matrixna)
          if('matrixna' %in% colnames(matrixna) ){
            colnames(matrixna) <- paste0(colnames(coluna), "1")
          }
          X <- cbind(X, matrixna)
        }
      }
      X <- X[-match("id", colnames(X) )]
      mod.samples <- mod.samples[colnames(X)]


      roc.out <- parallel::mclapply(seq_along(1:nrow(mod.samples)), function(samp){
        probs <- as.matrix(X) %*% as.vector(unlist( mod.samples[samp, ] ))
        tdROC::tdROC(X = probs[!is.na(probs)],
                     Y = test[[time]][!is.na(probs)],
                     delta = test[[status]][!is.na(probs)],
                     tau = max(test[[time]][test[[status]] ]),
                     n.grid = 1000)
      }, mc.cores = ncores )

      iroc <- sapply(roc.out, function(statistic){
        integrate_tdroc(statistic)
      })

      iroc <- mean(unlist( iroc))

      sens <- lapply(roc.out, function(sensibility){
        sensibility$ROC$sens
      })

      sens <- do.call(cbind, sens)
      sens <- apply(sens, 1, mean)

      spec <- lapply(roc.out, function(specificity){
        specificity$ROC$spec
      })

      spec <- do.call(cbind, spec)
      spec <- apply(spec, 1, mean)

      grid <- roc.out[[1]]$ROC$grid

      out <- list(iroc, sens, spec, grid)
      names(out) <- c( "iroc", "sens", "spec","grid")

      out

    } else {
      "Print object must be either coxph or stanreg"
    }
  }
}
#' @export
#' @rdname integrate_tdroc
integrate_tdroc  <-
  function( fit, ...) {
    fit$AUC[1] %>% unlist
  }
