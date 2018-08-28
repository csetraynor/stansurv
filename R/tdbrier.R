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
  timepoints =  seq(min(get_obs.time(testx)), max(get_obs.time(testx)), length.out = 100)

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

    basehaz <- do.call(cbind, basehaz.post)
    basehaz <- apply(X = basehaz, MARGIN = 1, FUN = median)

    obs.dat <- data.frame(haz = basehaz,
                           time = obs.time )
    if(method == "linear" | method == "constant"){
      # surv.approx <- approx(x = obs.time, y = surv.base, xout = timepoints, method = method)$y # only interpolation
      haz.approx <- Hmisc::approxExtrap(x = obs.time, y = basehaz, xout = timepoints, method = method)$y # with extrapolation
    } else if(method == "gp" | method == "gaussian process"){
      haz.bgp <- tgp:: bgp(X= obs.time, Z= basehaz)
      X <- data.frame(x1 = obs.time)
      XX <- data.frame(x1 = haz.approx)
      haz.bgp$Xsplit <- rbind(X, XX)
      haz.approx <- predict(object = haz.bgp, XX = timepoints)$ZZ.mean
    } else if(method == "gtp" | method == "treed gaussian process"){
      haz.btgp <- tgp:: btgp(X= obs.time, Z= basehaz)
      X <- data.frame(x1 = obs.time)
      XX <- data.frame(x1 = haz.approx)
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
  apply(surv.post, 1, median)
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
pred_surv <- function(fit, testx = NULL, timepoints = NULL, method = "linear", time = "time", status = "status", ncores = 1L){
  obs.time <- get_obs.time(fit$data)
  if(method != "linear" && method != "constant"){
    stop2(paste0("method ", method, " not implement.") )
  }
  if(is.null(timepoints)){
    timepoints <- test.time <- get_obs.time(testx)
  }
  basehaz.samples <- extract_basehaz_draws(fit)
  spline.basis <- extract_splines_bases(fit)

  if(!check_null_model(fit)){
    beta_draws <- extract_beta_draws(fit);
    varis.obs <- get_predictors(fit)

    ind.post <- parallel::mclapply(seq_along(1:nrow(testx)), function(i){
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
        linear_predictor <- varis.obs[i, ] %*% beta_draws[s, ]
        link_surv_lin.helper(p = linear_predictor, s = surv.approx)
      })
      surv.post <- do.call(cbind, surv.post)
      apply(surv.post, 1, median)
    }, mc.cores = ncores)
    ind.post <- do.call(rbind, ind.post)
    rownames(ind.post) <- rownames(testx)
  } else {
    # For null model
    surv.post <- lapply(seq_along(1:nrow(basehaz.samples)), function(s){
      basehaz.post <- sapply(seq_along(1:nrow(spline.basis)), function(n){
        basehaz.samples[s, ]  %*% spline.basis[n, ]
      })
      surv.base <- link_surv_base.helper(b = basehaz.post)

      if(method == "linear" | method == "constant"){
        # surv.approx <- approx(x = obs.time, y = surv.base, xout = timepoints, method = method)$y # only interpolation
        surv.approx <- Hmisc::approxExtrap(x = obs.time, y = surv.base, xout = timepoints, method = method)$y # with extrapolation
      } else if(method == "gp" | method == "gaussian process"){
        surv.bgp <- tgp:: btgp(X= obs.time, Z= surv.base)
        surv.bgp$Xsplit <- rbind(obs.time, timepoints)
        surv.approx <- predict(object = surv.bgp, XX = timepoints)$ZZ.mean
      } else {
        stop2(paste0("method ", method, " not supported"))
      }

      surv.approx[surv.approx > 1] <- 1
      surv.approx[surv.approx < 0] <- 0
      surv.approx
    })
    surv.post <- do.call(cbind, surv.post)
    surv.post <- apply(surv.post, 1, mean)
    #No difference in survival expected in the null model
    ind.post <- matrix(rep(surv.post, nrow(testx)),
                       nrow = nrow(testx),
                       ncol = n_distinct(timepoints),
                       byrow = TRUE)
  }
  return(ind.post)
}


