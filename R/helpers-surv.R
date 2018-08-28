# Function to predict the survival
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
link_surv <- function(fit, testx = NULL, timepoints = NULL, method = "linear", time = "time", status = "status"){
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

    surv.post <- lapply(seq_along(1:nrow(basehaz.samples)),
                        function(s){
      basehaz.post <- sapply(seq_along(1:nrow(spline.basis)),
                             function(n){
        basehaz.samples[s, ]  %*% spline.basis[n, ]
      })
      surv_base <- link_surv_base.helper(b = basehaz.post)

      surv.approx <- Hmisc::approxExtrap(x = obs.time, y = surv_base, xout = timepoints, method = method)$y
      surv.approx[surv.approx > 1] <- 1
      surv.approx[surv.approx < 0] <- 0

      linear_predictor <- varis.obs %*% beta_draws[s, ]
      link_surv_lin.helper(p = linear_predictor, s = surv.approx)
    })
  } else {
    surv.post <- lapply(seq_along(1:nrow(basehaz.samples)),
    function(s){
        basehaz.post <- sapply(seq_along(1:nrow(spline.basis)),
        function(n){
          basehaz.samples[s, ]  %*% spline.basis[n, ]})
        surv_base <- link_surv_base.helper(b = basehaz.post)
        surv.approx <- Hmisc::approxExtrap(x = obs.time, y = surv_base, xout = timepoints, method = method)$y
        surv.approx[surv.approx > 1] <- 1
        surv.approx[surv.approx < 0] <- 0
        surv.approx
  })
  }
  surv.post <- do.call(cbind, surv.post)
  surv.post
}

link_surv_base.helper <- function(b){
    #to translate baseline survival
    exp(-c(b ) )
}

link_surv_lin.helper <- function(p, s){
  s^(exp(c(p ) ) )
}

prepare_surv <- function(d, time = "time", status = "status"){
  if(!is_in("time", colnames(d) ) ){
    d$time <- d[[time]]
  }
  if(!is_in("status", colnames(d) ) ){
    d$status <- d[[status]]
  }
  d <- d[order(d$time), ]
  d
}

get_obs.time <- function(d){
  sort(  d$time[ as.logical(d$status) ] )
}

extract_basehaz_draws <- function(fit){
  as.matrix( rstan::extract(fit$stanfit)$basehaz_coefs )
}

extract_splines_bases <- function(fit){
  unlist( fit$basehaz$spline_basis )
}

check_null_model <- function(fit){
  !( n_distinct(fit$formula$allvars) > 2 )
}

extract_beta_draws <- function(fit){
  if(!check_null_model(fit)){
    as.matrix( rstan::extract(fit$stanfit)$beta )
  } else {
    NULL
  }
}



get_km.frame <- function(obs, strata = NULL, time = "time", status = "status"){

  if(!is.null(strata)){
    form <- as.formula(paste0("Surv(", time, ",", status, ") ~ ", paste(strata, collapse = "+")))
    mle.surv <- survfit(form, data = obs  )
    obs.mortality <- data.frame(time = mle.surv$time,
                                surv = mle.surv$surv,
                                strata = summary(mle.surv, censored=T)$strata)
    zeros <- data.frame(time=0, surv=1, strata=unique(obs.mortality$strata))
    obs.mortality <- rbind(obs.mortality, zeros)

    strata <- str_strata(obs.mortality)
    obs.mortality$strata <- NULL
    obs.mortality <- cbind( obs.mortality, strata)
  }else{
    form <- as.formula(paste0("survival::Surv(", time, ",", status, ") ~ 1"))
    mle.surv <- survival::survfit(form, data = obs  )
    obs.mortality <- data.frame(time = mle.surv$time,
                                surv = mle.surv$surv
                            )
    zeros <- data.frame(time=0, surv=1)
    obs.mortality <- rbind(obs.mortality, zeros)
  }

  obs.mortality
}


