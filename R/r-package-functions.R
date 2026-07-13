#Internal funcs####
#Function creates a vector with parameter names
get_latstat_parNames <- function(data, no_of_resp_options, latentState){
  #Mu pars
  mus <- unique(data[[latentState]])
  mus <- paste("mu.", mus, sep = "")

  #Find the fixed cut
  no_of_cuts <- no_of_resp_options - 1
  fixed_cut <- ifelse(no_of_cuts%% 2 == 0, (ceiling(no_of_cuts/2) + 1), ceiling(no_of_cuts/2))

  #Cutpoint pars
  cuts <- 1:no_of_cuts
  cuts <- cuts[cuts != fixed_cut]
  cuts <- paste("cut.", cuts, sep ="")

  parNames <- c(mus, cuts, "sd")

  return(parNames)
}

#function for running pmwg
run_all_stages <- function(sampler, iterations, particles, ...){
  burned <- pmwg::run_stage(sampler,
                            stage = "burn",
                            iter = iterations,
                            particles = particles,
                            ...)

  adapted <- pmwg::run_stage(burned,
                             stage = "adapt",
                             iter = iterations,
                             particles = particles,
                             ...)


  sampled <- pmwg::run_stage(adapted,
                             stage = "sample",
                             iter = iterations,
                             particles = particles,
                             ...)

  return(sampled)
}

#Function for wrangling mle list output
clean_mle_output <- function(results_list, subject_ids, bounds){
  lower_bound <- bounds[1]
  upper_bound <- bounds[2]

  #Bind all the results together
  results_df <- do.call(rbind, lapply(seq_along(results_list), function(i) {
    c(
      round(results_list[[i]]$par, 3),
      convergence = results_list[[i]]$convergence
    )
  }))


  #MANAGE RESULTS####
  results_df <- as.data.frame(results_df)

  #Add subject IDs
  rownames(results_df) <- subject_ids

  #Record convergence code
  results_df$convergence <- ifelse(results_df$convergence == 0, "Success",
                                   ifelse(results_df$convergence == 1, "Error - Limit Reached", paste0("Error - optim code: ", results_df$convergence)))

  #Override convergence code if bound reached
  check_cols <- setdiff(names(results_df), "convergence")

  has_upper <- apply(results_df[check_cols], 1, function(x) any(x == upper_bound))
  has_lower <- apply(results_df[check_cols], 1, function(x) any(x == lower_bound))

  results_df$convergence <- ifelse(
    has_upper & has_lower, "Both bounds reached",
    ifelse(
      has_upper, "Upper Bound Reached",
      ifelse(
        has_lower, "Lower Bound Reached",
        results_df$convergence
      )
    )
  )

  results_df$subject <- rownames(results_df)
  results_df <- results_df[, c("subject", setdiff(names(results_df), "subject"))]

  return(results_df)
}

#Function to check if all required args for ll_func have been included
check_ll_func_args <- function(ll_func, ...) {
  supplied <- list(...)
  fmls <- formals(ll_func)

  # Identify required args (no default)
  required_args <- names(fmls)[
    sapply(fmls, function(x) identical(x, quote(expr = )))
  ]

  # Remove args handled internally + `...`
  required_args <- setdiff(required_args, c("x", "data", "sample",  "subject_colname","..."))

  # Check missing
  missing_args <- setdiff(required_args, names(supplied))

  if (length(missing_args) > 0) {
    stop(
      paste("Missing required arguments for ll_func:",
            paste(missing_args, collapse = ", ")),
      call. = FALSE
    )
  }
}


#USEABLE FUNCS####
#' Make pmwg priors
#'
#' Creates priors for pmwg sampling
#'
#' @param data A trial-wise data frame
#' @param resp_opts A numeric value representing the number of response options available in the scale
#' @param est_directPref A boolean value reflecting whether a directionPref parameter should be estimated
#' @param latentState_colname A character string representing the latent state identifier variable
#' @param mu_prior A numeric value representing the prior for the mu parameter
#' @param centrePref_prior A numeric value representing the prior for the centrePref_prior parameter
#' @param oddsPref_prior A numeric value representing the prior for the oddsPref parameter
#' @param directionPref_prior A numeric value representing the prior for the directionPref parameter
#' @param threshold_prior A numeric value representing the prior for the thresholds in the latstat model
#' @param diag_prior A numeric value representing the prior for the diagonal
#' @param analysis A character string taking the value of 'lasars' or 'latstat' which indicates which model variant is being used
#' @return A list of priors
#'
#' @examples
#' test_priors <- make_priors(data = example_lasars_data,
#'                            resp_opts = 5,
#'                            latentState_colname = "subscale",
#'                            est_directPref = TRUE,
#'                            mu_prior = 0.5,
#'                            centrePref_prior = 0.3)
#' test_priors
#'
#' @export
make_priors <- function(data,
                        resp_opts,
                        latentState_colname,
                        est_directPref = NULL,
                        mu_prior = 0,
                        centrePref_prior = 0,
                        oddsPref_prior = 0,
                        directionPref_prior = 0,
                        threshold_prior = 0,
                        diag_prior = 1,
                        analysis = 'lasars'){

  if(analysis == 'lasars'){
    #Change the default for oddsPref to be 1/3 if 5 resp opts
    if(resp_opts == 5){
      oddsPref_prior = 1/3
    }

    if(is.null(est_directPref)){
      stop("est_directPref cannot be null in the lasars analysis.")
    }

    if(est_directPref){
      parNames <- c(paste0("mu.", unique(data[[latentState_colname]])), "centrePref", "oddsPref", "directionPref")
      priors <- list(theta_mu_mean = c(rep(mu_prior, length(unique(data[[latentState_colname]]))),
                                       centrePref_prior,
                                       oddsPref_prior,
                                       directionPref_prior),
                     theta_mu_var = diag(rep(diag_prior, length(parNames))))
      names(priors$theta_mu_mean) <- parNames



    } else{
      parNames <- c(paste0("mu.", unique(data[[latentState_colname]])), "centrePref", "oddsPref")
      priors <- list(theta_mu_mean = c(rep(mu_prior, length(unique(data[[latentState_colname]]))),
                                       centrePref_prior,
                                       oddsPref_prior),
                     theta_mu_var = diag(rep(diag_prior, length(parNames))))
      names(priors$theta_mu_mean) <- parNames

    }

    return(priors)

  } else if(analysis == 'latstat'){

    parNames <- get_latstat_parNames(data, resp_opts, latentState_colname)

    priors <- list(theta_mu_mean = c(rep(mu_prior, length(unique(data[[latentState_colname]]))),
                                     rep(threshold_prior, (sum(!grepl("mu", parNames)))-1), 1),
                   theta_mu_var = diag(rep(diag_prior, length(parNames))))

    names(priors$theta_mu_mean) <- parNames

    return(priors)

  } else{
    stop(paste0(analysis, " is not a valid analysis type.
                Please enter 'lasars' (default) or 'latstat"))
  }
}




#' Likelihood function for lasars Model
#'
#' The powerhouse of the lasars model - calculates the likelihood of the observed data,
#' given the set of input parameters.
#' If sample = TRUE, input parameters will be used to generate data.
#'
#' @param x A named vector of parameter estimates
#' @param data A trial-wise data frame
#' @param sample A boolean variable reflecting whether the ll func should evaluate likelihood or generate data
#' @param resp_opts A numeric value representing the number of response options available in the scale
#' @param subject_colname A character string representing the subject identifier variable
#' @param response_colname A character string representing the chosen Likert scale response
#' @param latentState_colname A character string representing the latent state identifier variable
#' @param direction_colname A character string representing the reverse-coding variable
#' @param rev_score_id A character string representing the level of the direction variable which indicates reverse-scored items
#' @return A likelihood value
#'
#' @examples
#' \dontrun{
#' ll_func <- lasars_ll_func(x = c("mu.extra" = 0.2,
#'                                "mu.open" = 0.1,
#'                                "mu.consc" = 0.1,
#'                                "mu.agree" = 0.2,
#'                                "mu.neuro" = 0.1,
#'                                "centrePref" = -.3,
#'                                "oddsPref" = 0.2,
#'                                "directPref" = -0.5),
#'                            data = example_lasars_data,
#'                            resp_opts = 5,
#'                            subject_colname = "subject",
#'                            response_colname = "response",
#'                            latentState_colname = "subscale",
#'                            direction_colname = "reverse",
#'                            rev_score_id = "TRUE")
#'}
#'
#' @export
lasars_ll_func <- function(x,
                           data,
                           sample = FALSE,
                           resp_opts,
                           subject_colname,
                           response_colname,
                           latentState_colname,
                           direction_colname = NULL,
                           rev_score_id = NULL) {
  #DEFINE PARAMETERS####
  #Fix sd
  sd = 1

  #DEFINE RESP PARS
  centrePref <- x["centrePref"]

  #get the oddsPref estimate (if it's being estimated)
  if("oddsPref" %in% names(x)){
    oddsPref <- x["oddsPref"]
  } else{
    oddsPref <- 0
  }

  #Define directionPref based on whether directionPref is in the pars vector
  if("directionPref" %in% names(x)){
    if(is.null(direction_colname) | is.null(rev_score_id)){
      stop("Must include direction_colname and rev_score_id inputs is estimating Direction Preference parameter")
    }
    directionPref <- x["directionPref"]
  } else{
    directionPref <- 0
  }

  #Func to make thresholds
  make_thresholds <- function(K, centrePref_alpha, oddsPref_gamma, directionPref_omega) {
    lambda <- numeric(K-1)
    G <- stats::pnorm(oddsPref_gamma)

    lowers <- 1:floor(K/2)
    lambda[lowers] <- exp(centrePref_alpha)*(-(K-2)/2)*G^(lowers-1) - directionPref_omega

    uppers <- ceiling((K+0.1)/2):(K-1)
    lambda[uppers] <- exp(centrePref_alpha)*((K-2)/2)*G^(K-1-uppers) - directionPref_omega

    if ((K%%2)==0) lambda[K/2] <- directionPref_omega

    return(lambda)
  }

  cutpoints <- make_thresholds(resp_opts, centrePref, oddsPref, directionPref)

  #DEFINE MUS
  scales <- sort(unique(data[[latentState_colname]]))

  mu <- c()
  for (i in 1:length(scales)){
    mu[i] <- x[paste0("mu.", scales[i])]
  }

  names(mu)<- scales


  #SAMPLE####
  if (sample){
    #replace response column with NAs
    data[[response_colname]] <- rep(NA, nrow(data))
    samples <- rep(NA, nrow(data))

    for(i in scales) {
      #If estimating directionPref:
      if("directionPref" %in% names(x)){
        for (j in unique(data[[direction_colname]])) {

          #do positively-scored items
          if (j != rev_score_id){
            use_trials <- data[[latentState_colname]] == i & data[[direction_colname]] == j
            samples[use_trials] <- stats::rnorm(sum(use_trials), mu[i], sd)

            #Do reverse-scored items
          } else{
            use_trials <- data[[latentState_colname]] == i & data[[direction_colname]] == j
            samples[use_trials] <- stats::rnorm(sum(use_trials), -mu[i], sd)
          }
        }
        #If not estimating directionPref:
      } else{
        use_trials <- data[[latentState_colname]] == i
        samples[use_trials] <- stats::rnorm(sum(use_trials), mu[i], sd)
      }
    }

    samples <- as.numeric(base::cut(samples, c(-Inf, cutpoints, Inf)))
    data[[response_colname]] <- samples
    return(data)

    #CALC LL
  } else {
    like <- rep(NA, nrow(data))

    for(i in scales) {

      #If estimating directionPref:
      if("directionPref" %in% names(x)){
        for (j in unique(data[[direction_colname]])){

          #Do positively-scored lls
          if(j != rev_score_id){
            use_trials <- data[[latentState_colname]] == i & data[[direction_colname]] == j
            p <- diff(c(0, stats::pnorm(c(cutpoints, Inf), mu[i], sd)))
            like[use_trials] <- p[data[[response_colname]][use_trials]]

            #Do reverse-scored lls
          } else {
            use_trials <- data[[latentState_colname]] == i & data[[direction_colname]] == j
            p <- diff(c(0, stats::pnorm(c(cutpoints, Inf), -mu[i], sd)))
            like[use_trials] <- p[data[[response_colname]][use_trials]]
          }
        }
        #If not estimating directionPref
      } else{
        use_trials <- data[[latentState_colname]] == i
        p <- diff(c(0, stats::pnorm(c(cutpoints, Inf), mu[i], sd)))
        like[use_trials] <- p[data[[response_colname]][use_trials]]
      }
    }

    # double check at this point there are no NAs in like!
    if(any(is.na(like))){
      stop(paste0("NA in subject: ", data[[subject_colname]][1], " likelihood vector."))
    }

    out <- sum(log(pmax(like, 1e-10))) # for protection against log 0 problems

    return(out)
  }
}


#' Likelihood function for Latent State Only Model
#'
#' The powerhouse of the latstat model - calculates the likelihood of the observed data,
#' given the set of input parameters.
#' If sample = TRUE, input parameters will be used to generate data.
#'
#' @param x A named vector of parameter estimates
#' @param data A trial-wise data frame
#' @param sample A boolean variable reflecting whether the ll func should evaluate likelihood or generate data
#' @param resp_opts A numeric value representing the number of response options available in the scale
#' @param subject_colname A character string representing the subject identifier variable
#' @param response_colname A character string representing the chosen Likert scale response
#' @param latentState_colname A character string representing the latent state identifier variable
#' @param direction_colname A character string representing the reverse-coding variable
#' @param rev_score_id A character string representing the level of the direction variable which indicates reverse-scored items
#' @return A likelihood value
#'
#' @examples
#' \dontrun{
#' ll_func <- latstat_ll_func(x = c("mu.extra" = 0.2,
#'                                  "mu.open" = 0.1,
#'                                  "mu.consc" = 0.1,
#'                                  "mu.agree" = 0.2,
#'                                  "mu.neuro" = 0.1,
#'                                  "cut.1" = -0.3,
#'                                  "cut.2" = -0.2,
#'                                  "cut.4" = 0.5,
#'                                  sd = 0.3),
#'                            data = example_lasars_data,
#'                            resp_opts = 5,
#'                            subject_colname = "subject",
#'                            response_colname = "response",
#'                            latentState_colname = "subscale",
#'                            direction_colname = "reverse",
#'                            rev_score_id = "TRUE")
#'}
#' @export
latstat_ll_func <- function(x,
                            data,
                            sample = FALSE,
                            resp_opts,
                            subject_colname,
                            response_colname,
                            latentState_colname,
                            direction_colname,
                            rev_score_id) {

  no_of_cuts <- resp_opts - 1
  fixed_cut <- ifelse(no_of_cuts%% 2 == 0, (ceiling(no_of_cuts/2) + 1), ceiling(no_of_cuts/2))


  #Define sd
  sd = exp(x["sd"])


  ###DEFINE CUTPOINTS
  #assign cut.1
  cut.1 <- x["cut.1"]
  cutpoints <- cut.1

  #define cuts up to fixed cut
  if(fixed_cut != 2){
    #loop through until you reach the fixed cut
    for (point in 2:(fixed_cut-1)) {
      #and define cut points
      assign(paste0("cut.", point), get(paste0("cut.", (point-1))) + exp(x[paste0("cut.", point)]))
      cutpoints <- c(cutpoints, get(paste0("cut.", point)))
      names(cutpoints)[point] <- paste0("cut.", point)
    }
  }

  #define fixed cut
  assign(paste0("cut.", fixed_cut), get(paste0("cut.", (fixed_cut -1))) + 1)
  cutpoints <- c(cutpoints, get(paste0("cut.", fixed_cut)))
  names(cutpoints)[fixed_cut] <- paste0("cut.", fixed_cut)


  #loop through remaining cuts and define them
  for (point in (fixed_cut+1):no_of_cuts) {
    assign(paste0("cut.", point), get(paste0("cut.", (point-1))) + exp(x[paste0("cut.", point)]))
    cutpoints <- c(cutpoints, get(paste0("cut.", point)))
    names(cutpoints)[point] <- paste0("cut.", point)
  }

  #Mus
  scales <- sort(unique(data[[latentState_colname]]))

  mu <- c()
  for (i in 1:length(scales)){
    mu[i] <- x[paste0("mu.", scales[i])]
  }

  names(mu)<- scales

  if(!sample){
    #GET LL
    like <- rep(NA, length(data[[response_colname]]))

    if(!is.null(direction_colname)){
      for(i in scales) {
        for (j in unique(data[[direction_colname]])){
          if(j != rev_score_id){
            use_trials <- data[[latentState_colname]] == i & data[[direction_colname]] == j
            p <- diff(c(0, stats::pnorm(c(cutpoints, Inf), mu[i], sd)))
            like[use_trials] <- p[data[[response_colname]][use_trials]]
          } else if (j == rev_score_id) {
            use_trials <- data[[latentState_colname]] == i & data[[direction_colname]] == j
            p <- diff(c(0, stats::pnorm(c(cutpoints, Inf), -mu[i], sd)))
            like[use_trials] <- p[data[[response_colname]][use_trials]]
          }
        }
      }
    } else{
      for(i in scales){
        use_trials <- data[[latentState_colname]] == i
        p <- diff(c(0, stats::pnorm(c(cutpoints, Inf), mu[i], sd)))
        like[use_trials] <- p[data[[response_colname]][use_trials]]
      }
    }

    out <- sum(log(pmax(like, 1e-10)))

    return(out)


  } else{

    #replace response column with NAs
    data[[response_colname]] <- rep(NA, nrow(data))
    samples <- rep(NA, nrow(data))

    for(i in scales) {
      for (j in unique(data[[direction_colname]])) {
        if (j != rev_score_id){
          #use_trials is a vector which indicates which trials to use this time around
          use_trials <- data[[latentState_colname]] == i & data[[direction_colname]] == j
          #samples gets a random sample from a norm dist for all trials which satisfy this loop's conditions
          samples[use_trials] <- stats::rnorm(sum(use_trials), mu[i], sd)
        } else if (j == rev_score_id) {
          use_trials <- data[[latentState_colname]] == i & data[[direction_colname]] == j
          samples[use_trials] <- stats::rnorm(sum(use_trials), -mu[i], sd)
        }
      }
    }

    #the cut function identifies which bin the sample falls into (eg. between cut.1 and cut.2)
    #as.numeric turns these bins into numbers (eg. 2) which correspond with response options
    samples <- as.numeric(cut(samples, c(-Inf, cutpoints, Inf)))
    data[[response_colname]] <- samples
    return(data)
  }
}


#' Fit the lasars model
#'
#' Estimates latent-state and response-style parameters from
#' Likert-scale survey data.
#'
#' @param data A trial-wise data frame
#' @param resp_opts A numeric value representing the number of response options available in the scale
#' @param est_directPref A boolean value reflecting whether an directionPref parameter should be estimated
#' @param subject_colname A character string representing the subject identifier variable
#' @param response_colname A character string representing the chosen Likert scale response
#' @param latentState_colname A character string representing the latent state identifier variable
#' @param direction_colname A character string representing the reverse-coding variable
#' @param rev_score_id A character string representing the level of the direction variable which indicates reverse-scored items
#' @param priors A list of priors
#' @param sampling_method A character string which reflects the sampling approach to be taken.
#' @param iterations A numeric value reflecting the number of sampling iterations to run
#' @param particles A numeric value reflecting the number of proposed particles on each sampling iteration
#' @param ... Additional parameters to pass into the pmwg run_stage calls
#' @return A pmwg sampled object
#'
#' @examples
#' \donttest{
#' result <- run_lasars(
#'   data = example_lasars_data,
#'   resp_opts = 5,
#'   est_directPref = TRUE,
#'   subject_colname = "subject",
#'   response_colname = "response",
#'   latentState_colname = "subscale",
#'   direction_colname = "reverse",
#'   rev_score_id = "TRUE",
#'   iterations = 100,
#'   particles = 10
#' )
#'}
#'
#' @export
run_lasars <- function(data,
                       resp_opts,
                       est_directPref,
                       subject_colname,
                       response_colname,
                       latentState_colname,
                       direction_colname = NULL,
                       rev_score_id = NULL,
                       priors = NULL,
                       sampling_method = 'pmwg',
                       iterations = 3000,
                       particles = 50,
                       ...){
  #STOPS FOR INCORRECT OR MISSING DATA####
  column_vars <- c(
    subject_colname   = subject_colname,
    response_colname  = response_colname,
    latentState_colname  = latentState_colname
  )

  # find which are not in the data columns
  missing_cols <- names(column_vars)[!column_vars %in% colnames(data)]

  # stop if any are missing
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "The following inputs do not match any column names in `data`: ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }

  #Test for NULLS and NAS
  input_vars <- c(column_vars, "resp_opts" = resp_opts)
  invalid_cols <- names(input_vars)[
    vapply(column_vars, function(x) {
      is.null(x) || (length(x) == 1 && is.na(x))
    }, logical(1))
  ]

  if (length(invalid_cols) > 0) {
    stop(
      paste0(
        "The following inputs are NULL or NA: ",
        paste(invalid_cols, collapse = ", ")
      )
    )
  }

  #No direction
  if(est_directPref & is.null(direction_colname)){
    stop("Must supply 'direction_colname' if estimating directionPref.
         This is the variable which indicates whether items are positively- or reverse-scored.")
  }

  #Check if direction colname is real
  if (!is.null(direction_colname) && !(direction_colname %in% colnames(data))) {
    stop("The following inputs do not match any column names in `data`: direction_colname")
  }

  #No rev_scored_id
  if(est_directPref & is.null(rev_score_id)){
    stop(paste0("Must supply reverse scoring identifier.
         This is the string in your '", direction_colname , "' variable which indicates an item is to be reverse-scored."))
  }

  #Check if rev_score_id is real
  if (!is.null(rev_score_id) && !(rev_score_id %in% unique(data[[direction_colname]]))) {
    stop(paste0("No level: '", rev_score_id, "' found in ", direction_colname, " column."))
  }

  #Stop if more than 2 levels in direction
  if(!is.null(direction_colname) && length(unique(data[[direction_colname]])) > 2){
    stop(paste0("'", direction_colname, "' should only have two levels: one identifying postively-worded items
                and one for reverse-scored items.
                Currently has ", length(unique(data[[direction_colname]])), " levels: ", unique(data[[direction_colname]])))
  }

  #Stop if any responses larger (or smaller) than resp_opts
  if (
    max(data[[response_colname]], na.rm = TRUE) > resp_opts |
    min(data[[response_colname]], na.rm = TRUE) < 1
  ) {
    stop("Responses must be integers in 1:resp_opts")
  }

  #DATA CLEANING####
  data[[response_colname]] <- as.numeric(data[[response_colname]])

  #Remove NA response vals
  data <- data[!is.na(data[[response_colname]]), ]
  data <- data[!is.null(data[[response_colname]]), ]
  rev_score_id <- as.character(rev_score_id)

  if(sampling_method == 'pmwg'){
    #Make default priors if none passed in
    if(is.null(priors)){
      priors <- make_priors(data = data,
                            resp_opts = resp_opts,
                            est_directPref = est_directPref,
                            latentState_colname = latentState_colname)
    }

    #Rename subject column
    names(data)[names(data) == subject_colname] <- "subject"

    #MAKE LL_FUNC####
    ll_func <- purrr::partial(lasars_ll_func,
                              resp_opts = resp_opts,
                              subject_colname = subject_colname,
                              response_colname = response_colname,
                              latentState_colname = latentState_colname,
                              direction_colname = direction_colname,
                              rev_score_id = rev_score_id)

    #RUN PMWG####
    sampler <- pmwg::pmwgs(data = data,
                           pars = names(priors$theta_mu_mean),
                           ll_func = ll_func,
                           prior = priors)

    sampler <- pmwg::init(sampler)

    #Run
    sampled <- run_all_stages(sampler, iterations, particles, ...)
    return(sampled)

  } else if(sampling_method == 'mle'){
    #DEFINE STARTPOINTS (dynamically)####
    #mu pars
    mu_starts <- c(rep(0, length(unique(data[[latentState_colname]]))))
    names(mu_starts) <- c(paste0('mu.', unique(data[[latentState_colname]])))

    #Response pars
    if(est_directPref == TRUE){
      resp_starts <- c("centrePref" = 0,
                       "oddsPref" = 0,
                       "directionPref" = 0)
    } else{
      resp_starts <- c("centrePref" = 0,
                       "oddsPref" = 0)
    }

    start_vals <- c(mu_starts, resp_starts)

    #Set upper and lower bounds
    upper_bound <- 20
    lower_bound <- -20

    #RUN OPTIM###
    #Run optim over each subject's data and add to a list
    subject_ids <- unique(data[[subject_colname]])
    results_list <- lapply(subject_ids, function(subj) {

      this_subj_data <- data[data[[subject_colname]] == subj, ]

      stats::optim(
        par = start_vals,
        fn = lasars_ll_func,
        data = this_subj_data,
        control = list(fnscale = -1),
        latentState_colname = latentState_colname,
        response_colname = response_colname,
        direction_colname = direction_colname,
        rev_score_id = rev_score_id,
        resp_opts = resp_opts,
        method = "L-BFGS-B",
        lower = rep(lower_bound, length(start_vals)),
        upper = rep(upper_bound, length(start_vals))
      )
    })

    results_df <- clean_mle_output(results_list,
                                     subject_ids,
                                     c(lower_bound, upper_bound))

    return(results_df)

  } else{
    stop(paste0("Error: Sampling method ", sampling_method, " is not an option for this analysis.
          Please enter either 'pmwg' (default) or 'mle'."))
  }
}

#' Runs the latent state only model
#'
#' This function runs the model over the data
#'
#' @param data A trial-wise data frame
#' @param resp_opts A numeric value representing the number of response options available in the scale
#' @param subject_colname A character string representing the subject identifier variable
#' @param response_colname A character string representing the chosen Likert scale response
#' @param latentState_colname A character string representing the latent state identifier variable
#' @param direction_colname A character string representing the reverse-coding variable
#' @param rev_score_id A character string representing the level of the direction variable which indicates reverse-scored items
#' @param priors A list of priors
#' @param sampling_method A character string which reflects the sampling approach to be taken. Can take the value 'mle' or 'pmwg'.
#' @param iterations A numeric value reflecting the number of sampling iterations to run
#' @param particles A numeric value reflecting the number of proposed particles on each sampling iteration
#' @param ... Additional parameters to pass into the pmwg run_stage calls
#' @return A pmwg sampled object
#'
#' @examples
#' \dontrun{
#' result <- run_latstat(
#'   data = example_lasars_data,
#'   resp_opts = 5,
#'   subject_colname = "subject",
#'   response_colname = "response",
#'   latentState_colname = "subscale",
#'   direction_colname = "reverse",
#'   rev_score_id = "TRUE"
#' )
#'}
#' @export
#Run original version
run_latstat <- function(data,
                        resp_opts,
                        subject_colname,
                        response_colname,
                        latentState_colname,
                        direction_colname,
                        rev_score_id,
                        priors = NULL,
                        sampling_method = 'pmwg',
                        iterations = 3000,
                        particles = 50,
                        ...){
  #STOPS FOR INCORRECT DATA
  # your named vector
  column_vars <- c(
    subject_colname   = subject_colname,
    response_colname  = response_colname,
    latentState_colname  = latentState_colname,
    direction_colname = direction_colname
  )

  # find which are not in the data columns
  missing_cols <- names(column_vars)[!column_vars %in% colnames(data)]

  # stop if any are missing
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "The following inputs do not match any column names in `data`: ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }

  #Test for NULLS and NAS
  input_vars <- c(column_vars, "resp_opts" = resp_opts)
  invalid_cols <- names(input_vars)[
    vapply(column_vars, function(x) {
      is.null(x) || (length(x) == 1 && is.na(x))
    }, logical(1))
  ]

  if (length(invalid_cols) > 0) {
    stop(
      paste0(
        "The following inputs are NULL or NA: ",
        paste(invalid_cols, collapse = ", ")
      )
    )
  }

  #Check if direction colname is real
  if (!is.null(direction_colname) && !(direction_colname %in% colnames(data))) {
    stop("The following inputs do not match any column names in `data`: direction_colname")
  }

  #No rev_scored_id
  if(!is.null(direction_colname) && is.null(rev_score_id)){
    stop(paste0("Must supply reverse scoring identifier for this analysis.
         This is the string in your '", direction_colname , "' variable which indicates an item is to be reverse-scored.
                If none exists, consider using the lasars model."))
  }

  #Check if rev_score_id is real
  if (!is.null(rev_score_id) && !(rev_score_id %in% unique(data[[direction_colname]]))) {
    stop(paste0("No level: '", rev_score_id, "' found in ", direction_colname, " column."))
  }

  #Stop if more than 2 levels in direction
  if(!is.null(direction_colname) && length(unique(data[[direction_colname]])) > 2){
    stop(paste0("'", direction_colname, "' should only have two levels: one identifying postively-worded items
                and one for reverse-scored items.
                Currently has ", length(unique(data[[direction_colname]])), " levels: ", unique(data[[direction_colname]])))
  }

  #Stop if any responses larger (or smaller) than resp_opts
  if (
    max(data[[response_colname]], na.rm = TRUE) > resp_opts |
    min(data[[response_colname]], na.rm = TRUE) < 1
  ) {
    stop("Responses must be integers in 1:resp_opts")
  }

  data[[response_colname]] <- as.numeric(data[[response_colname]])
  data <- data[!is.na(data[[response_colname]]), ]
  data <- data[!is.null(data[[response_colname]]), ]
  rev_score_id <- as.character(rev_score_id)

  if(sampling_method == 'pmwg'){
    #MAKE PRIORS####
    #Make default priors if none passed in
    if(is.null(priors)){
      priors <- make_priors(data = data,
                            resp_opts = resp_opts,
                            latentState_colname = latentState_colname,
                            analysis = 'latstat')
    }

    #Rename subject column
    names(data)[names(data) == subject_colname] <- "subject"

    #MAKE LL_FUNC####
    ll_func <- purrr::partial(latstat_ll_func,
                              resp_opts = resp_opts,
                              subject_colname = subject_colname,
                              response_colname = response_colname,
                              latentState_colname = latentState_colname,
                              direction_colname = direction_colname,
                              rev_score_id = rev_score_id)


    #RUN PMWG####
    sampler <- pmwg::pmwgs(data = data,
                           pars = names(priors$theta_mu_mean),
                           ll_func = ll_func,
                           prior = priors)

    sampler <- pmwg::init(sampler)

    #Run
    sampled <- run_all_stages(sampler, iterations, particles, ...)
    return(sampled)

  } else if(sampling_method == 'mle'){
    #DEFINE STARTPOINTS (dynamically)####
    parNames <- get_latstat_parNames(data, resp_opts, latentState_colname)
    start_vals <- c(rep(0, length(parNames) - 1), 1)

    names(start_vals) <- parNames

    #Set upper and lower bounds
    upper_bound <- 20
    lower_bound <- -20

    #RUN OPTIM####
    subject_ids <- unique(data[[subject_colname]])

    #Run optim over each subject's data and add to a list
    results_list <- lapply(subject_ids, function(subj) {

      this_subj_data <- data[data[[subject_colname]] == subj, ]

      stats::optim(
        par = start_vals,
        fn = latstat_ll_func,
        data = this_subj_data,
        control = list(fnscale = -1),
        latentState_colname = latentState_colname,
        response_colname = response_colname,
        direction_colname = direction_colname,
        rev_score_id = rev_score_id,
        resp_opts = resp_opts,
        method = "L-BFGS-B",
        lower = rep(lower_bound, length(start_vals)),
        upper = rep(upper_bound, length(start_vals))
      )
    })

    results_df <- clean_mle_output(results_list, subject_ids, c(lower_bound, upper_bound))

    return(results_df)

  } else{
    stop(paste0("Error: Sampling method ", sampling_method, " is not an option for this analysis.
                Please enter either 'pmwg' (default) or 'mle'."))
  }
}

#' Generates posterior samples
#'
#' This function generates posterior-predictive data
#'
#' @param sampled A post-sampling object with parameter estimates. Should be the output of either run_lasars() or run_latstat()
#' @param n the number of posterior samples being generated
#' @param ll_func The likelihood function being used
#' @param sampling_method A character string which reflects the sampling approach to be taken. Can take the value 'mle' or 'pmwg'
#' @param original_df Necessary when sampling_method = 'mle'. Should be the original trial-wise data frame sampling was performed on
#' @param subject_colname A character string representing the subject identifier variable
#' @param ... Additional parameters to pass into the pmwg run_stage calls
#' @return A data frame
#'
#' @examples
#' \dontrun{
#' pp_data <- gen_pp_data(sampled = result)
#' }
#'
#' @export
gen_pp_data <- function(sampled,
                        n = 20,
                        ll_func = sampled$ll_func,
                        sampling_method = 'pmwg',
                        original_df = NULL,
                        subject_colname = "subject",
                        ...){

  #Setup data and iterations
  if(sampling_method == 'pmwg'){
    data <- sampled$data
    subjects <- unique(sampled$data$subject)
    sampled_stage <- length(sampled$samples$stage[sampled$samples$stage == "sample"])
    iterations <- round(seq(from = (sampled$samples$idx - sampled_stage),
                            to = sampled$samples$idx,
                            length.out = n))

  } else if(sampling_method == 'mle'){
    #Check if ll_func has been entered
    if(missing(ll_func)){
      stop("ll_func input must be included for mle analyses")
    }

    #Make sure they've entered a df
    if(is.null(original_df)){
      stop("If using mle sampling method, original_df cannot be NULL.
           This should reflect the data frame structure samples should follow.")
    }

    data <- original_df
    subjects <- unique(sampled$subject)

  } else{
    stop(paste0("Sampling method '", sampling_method, "' is not a valid option. Please choose either 'pmwg' or 'mle'."))
  }

  #Check for missing ll_func args
  check_ll_func_args(ll_func, ...)

  if(!subject_colname %in% colnames(data)){
    stop("Please enter a valid subject_colname value")
  }

  #Create empty pp_data frame
  pp_data <- c()

  #Loop through each person
  for(s in 1:length(subjects)){
    subj <- subjects[s]

    #subset to just their data
    print(paste0("subject ", s))
    this_subj_data <- data[data[[subject_colname]] == subj, ]

    #Get
    if(sampling_method == 'mle'){
      this_subj_pars <- sampled[sampled$subject == subj, setdiff(colnames(sampled), c("subject", "convergence"))]
      this_subj_pars <- unlist(this_subj_pars[1, ], use.names = TRUE)
    }

    this_subj_pp <- c()

    #loop through iterations
    for(i in 1:n){
      if(sampling_method == 'pmwg'){
        this_subj_pars <- sampled$samples$alpha[,s, iterations[i]]

        #Generate data
        this_iter_pp <- ll_func(x = this_subj_pars,
                                data = this_subj_data,
                                sample = TRUE,
                                ...)
      } else{
        #Generate data
        this_iter_pp <- ll_func(x = this_subj_pars,
                                data = this_subj_data,
                                sample = TRUE,
                                subject_colname = subject_colname,
                                ...)
      }

      #Join this data to existing subj data
      this_iter_pp$iter <- i
      this_subj_pp <- rbind(this_subj_pp, this_iter_pp)
    }

    this_subj_pp$subject <- subj

    #Join to overall data
    pp_data <- rbind(pp_data, this_subj_pp)
  }
  return(pp_data)
}





#' Example dataset for lasars
#'
#' A small Big-5 dataset used in examples and vignettes.
#'
#' @format A data frame with 50 participants and 5 latent states:
#' \describe{
#'   \item{subject}{Subject identifier}
#'   \item{item}{Survey item identifier}
#'   \item{response}{Chosen Likert scale response}
#'   \item{subscale}{Latent state identifier}
#'   \item{reverse}{Identifies positively- and reverse-scored items}
#' }
#'
#' @source Open Source Psychometrics
"example_lasars_data"
