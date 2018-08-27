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

tdbrier <- function(holdout, fit, time = "time", status = "status", ncores = 1L){
  #select timepoints with observed event
  holdout <- prepare_surv(holdout, time = time, status = status)
  timepoints <-  seq(min(holdout[[time]]), max(holdout[[time]]), length.out = 100)

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
    object <- cox.brier

    rsquaredbs <- suppressWarnings ( pec::R2(object = cox.brier)$AppErr )
    out <- list(object, rsquaredbs, ibrier, apperror, ref, time)
    names(out) <- c("brier.object", "rsquaredbs","ibrier", "apperror", "reference","time")
    out

  } else if (any( class(fit$stanfit) == "stanfit" ) ){
    ind.pred <- pred_surv(fit = fit, testx = holdout, timepoints = timepoints, ncores = ncores)

    ba.brier <- suppressMessages( pec::pec(ind.pred, Surv(time, status) ~ 1,
                                        data = holdout,
                                        maxtime = max(timepoints),
                                        exact = FALSE,
                                        exactness = 99L ) )
    apperror <- ba.brier$AppErr$matrix
    time <- timepoints
    ref <- ba.brier$AppErr$Reference
    ibrier <- integrate_tdbrier(ba.brier)
    object <- ba.brier

    rsquaredbs <- suppressWarnings ( pec::R2(object = ba.brier)$AppErr )
    out <- list(object, rsquaredbs, ibrier, apperror, ref, time)
    names(out) <- c("brier.object", "rsquaredbs","ibrier", "apperror", "reference","time")
    out
  } else {
    stop2(paste0("Error: No method for evaluating predicted probabilities from objects in class", class(fit)))
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
      surv_base <- link_surv_base.helper(b = basehaz.post)
      surv.approx <- Hmisc::approxExtrap(x = obs.time, y = surv_base, xout = timepoints, method = method)$y
      # surv.approx <- approx(x = obs.time, y = surv_base, xout = timepoints, method = method)$y
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
pred_surv2 <- function(fit, testx = NULL, timepoints = NULL, time = "time", status = "status", ncores = 1L){
  obs.time <- get_obs.time(fit$data)
  if(is.null(timepoints)){
    timepoints <- test.time <- get_obs.time(testx)
  }
  basehaz.samples <- extract_basehaz_draws(fit)
  spline.basis <- extract_splines_bases(fit)

  if(!check_null_model(fit)){
    beta_draws <- extract_beta_draws(fit);
    varis.obs <- get_predictors(fit)

surv.base <- lapply(seq_along(1:nrow(basehaz.samples)), function(s){
  basehaz.post <- sapply(seq_along(1:nrow(spline.basis)), function(n){
    basehaz.samples[s, ]  %*% spline.basis[n, ]
    })
  link_surv_base.helper(b = basehaz.post)
})
surv.base <- do.call(cbind, surv.base)
surv.base.mean <- apply(X = surv.base, MARGIN = 1, FUN = median)

obs.surv <- data.frame(surv = surv.base.mean,
                       time = obs.time )
#st_mod <- rstanarm::stan_gamm4(surv ~ s(time, bs='gp'), data = obs.surv) Gaussian procces??
#st_mod <- rstanarm::stan_glm(surv ~ time, data = obs.surv,  family = gaussian(), prior = cauchy(),  prior_intercept = cauchy() ) fail linear model just plot and see that is not linear
# stan_file <- system.file('stan', 'weibull_survival_null_model.stan', package =  'biostan') try weibull?
options(mc.cores = parallel::detectCores() )
rstan_options(auto_write = TRUE)
fit1 <- brms::brm(surv ~ gp(time), data = obs.surv)

pred.surv.base <- rstanarm::posterior_predict(object = st_mod, newdata = data.frame(time = timepoints) )
pred.surv.base <- exp( apply(pred.surv.base, 2, mean) )
pred.surv.base[pred.surv.base > 1] <- 1
pred.surv.base[pred.surv.base < 0] <- 0

ind.post <- parallel::mclapply(seq_along(1:nrow(testx)),
                               function(i){
  surv.post <- sapply(seq_along(1:nrow(beta_draws)), function(s){
    linear_predictor <- varis.obs[i, ] %*% beta_draws[s, ]
    link_surv_lin.helper(p = linear_predictor, s = pred.surv.base)
  })
  apply(surv.post, 1, mean)
}, mc.cores = ncores)

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


