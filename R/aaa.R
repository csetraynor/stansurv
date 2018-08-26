#get predictors from formula

get_predictors <- function(fit){
  if(!check_null_model(fit)){
    obs.p <- extract_predictors(fit)
    factor_vars <- get_factors(fit = fit)
    obs.p[ ,match(colnames(factor_vars), colnames(obs.p))] <- factor_vars
    if(n_distinct(colnames(obs.p)) > 1000){
      obs.p <- preproc_big_mat(obs.p)
    } else if(ncol(obs.p) == 1 && nlevels(obs.p[[1]]) < 3) {
      obs.p <-preproc_one_mat(obs.p)
    } else {
      obs.p <- model.matrix(~ ., obs.p)[ ,-1]
    }
    as.matrix(obs.p)
  } else {
    NULL
  }
}

# Deal with factors in the posterior fit
get_factors <- function(p, v, fit = NULL){
  if(!is.null(fit)){
    p <- extract_predictors(fit)
    v <- extract_vars(fit)
  }
  #from formula
  fact.1 <- p[grep("factor\\(", v)]
  names_fact.1 <- colnames(fact.1)
  if(n_distinct(names_fact.1)){
    fact.1 <- lapply(seq_along(colnames(fact.1) ), function(i) as.factor(fact.1[ ,i]))
    # fact.1 <- lapply(fact.1, addNA)
    if(n_distinct(fact.1) > 1){
      fact.1 <- do.call(cbind.data.frame, fact.1)
    } else {
      fact.1 <- as.data.frame(fact.1)
    }
    colnames(fact.1) <- names_fact.1
    fact.1
  } else {
    NULL
  }
}


extract_predictors <- function(fit){
  if(!check_null_model(fit)){
    fit$data[ fit$formula$allvars[-(1:2)] ]
  } else {
    NULL
  }
}


extract_vars <- function(fit){
  vars <- unlist( strsplit(deparse(fit$formula$rhs), "[+]") )
  vars <- gsub(" ", "", vars) #remove space
  vars
}



preproc_one_mat <- function(o){
  coluna <- as.data.frame(o)
  colnames_coluna <- colnames(coluna)
  matrixna <- model.matrix(~ ., coluna)[ ,-1]
  matrixna <- as.data.frame(matrixna)
  colnames(matrixna) <- colnames_coluna
  as.data.frame(matrixna)
}


preproc_big_mat <- function(big_mat, patient_id = NULL){

  #Convert to factor
  big_mat <- as.data.frame(big_mat)
  big_mat <- lapply(big_mat, as.factor)
  #add NA as factor
  big_mat <- lapply(big_mat, addNA)

  X <- lapply(big_mat, function(i){
    preproc_one_mat(i)
  })
  big_mat <- do.call(cbind.data.frame, X)
  #remove zero var columns
  Filter(function(x)!all(is.na(x)), big_mat)

  return(big_mat)
}

#' String method for strata in surv analysis
#'
#' This function is a string method to facilitate ploting.
#'
#' @param c a character frame \cr
#' @return d a clean dataset
#' @export
#' @importFrom magrittr %>%
str_strata <- function(c){
  test <- strsplit(as.character( c$strata ), "\\,")
  test <-as_data_frame( do.call(rbind, test) )
  colnames(test) <- gsub("=.*| ","",  test[1, ] )
  strata_replace <- function(x){
    x <- gsub(".*=| ", "", x)
    x
  }
  test <- test %>%
    dplyr::mutate_all(funs(strata_replace))
  test
}
