pred_surv <- function(fit, testx = NULL, timepoints = NULL, method = "linear", time = "time", status = "status"){
  obs.time <- get_obs.time(fit$data)
  if(!is.null(timepoints)){
    timepoints <- test.time <- get_obs.time(testx)
  }
  basehaz.samples <- extract_basehaz.samples(fit)
  spline.basis <- extract_splines.bases(fit)

  surv.post <- lapply(seq_along(1:nrow(basehaz.samples)), function(s){
    basehaz.post <- sapply(seq_along(1:nrow(spline.basis)), function(n){
     basehaz.samples[s, ]  %*% spline.basis[n, ]
    })
    link.surv(basehaz.post)
  })

  surv.post <- do.call(cbind, surv.post)

  # not other methods apart from interpolation of surv mean yet
  if(method == "linear" | method == "constant"){
    surv.mean <- apply(surv.post, 1, mean)
    surv.test <- approx(x = obs.time, y = surv.mean, xout = timepoints, method = method)$y
  }
  return(surv.test)
}

link_surv.helper <- function(basehaz, predictionind = NULL ){
  if(!is.null(p)){
    exp( -(b))^exp(p)
  } else {
    exp( -(b))
  }
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

extract_predictors <- function(fit){
  if(!check_null_model(fit)){
    obs.p <- fit$data[fit$formula$allvars[-(1:2)] ]
    obs.vars <- extract_vars(fit)
    factor.vars <- obs.p[grep("factor\\(", obs.vars)]
    factor.vars2 <- obs.p[ ,is.fact] <- sapply(obs.p, is.factor)
    if(ncol(factor.vars) == 0){
      return(obs.p)
    } else {

    }

  } else {
    NULL
  }
}

extract_vars <- function(fit){
  vars <- unlist( strsplit(deparse(fit$formula$rhs), "[+]") )
  vars <- gsub(" ", "", vars) #remove space
  vars
}

get_surv.post <- function(baseline_hazard, spline_basis){
  lapply(seq_along(1:nrow(baseline_hazard)), function(s){
    basehaz_post <- sapply(seq_along(1:nrow(spline_basis)), function(n){
      log( basehaz.samples[s, ] ) %*% spline_basis[n, ]
    })
    link.surv(basehaz_post)
  })
  do.call(cbind, surv_post)
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
