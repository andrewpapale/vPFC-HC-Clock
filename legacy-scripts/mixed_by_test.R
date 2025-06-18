mixed_by_test <- function (data, outcomes = NULL, rhs_model_formulae = NULL, model_formulae = NULL, 
          split_on = NULL, external_df = NULL, external_merge_by = NULL, 
          padjust_by = "term", padjust_method = "BY", outcome_transform = NULL, 
          scale_predictors = NULL, ncores = 1L, cl = NULL, refit_on_nonconvergence = 3, 
          tidy_args = list(effects = "fixed", conf.int = TRUE), lmer_control = lmerControl(optimizer = "nloptwrap"), 
          calculate = c("parameter_estimates_reml", "parameter_estimates_ml", 
                        "fit_statistics"), return_models = FALSE, emmeans_spec = NULL, 
          emtrends_spec = NULL) 
{
  requireNamespace("data.table")
  requireNamespace("dplyr")
  requireNamespace("lme4")
  requireNamespace("lmerTest")
  requireNamespace("foreach")
  requireNamespace("doParallel")
  requireNamespace("iterators")
  requireNamespace("broom.mixed")
  requireNamespace("formula.tools")
  if (checkmate::test_data_frame(data)) {
    single_df <- TRUE
  }
  else {
    if (!checkmate::test_character(data)) {
      stop("If data is not a data.frame, it should be a set of file names.")
    }
    checkmate::assert_file_exists(data)
    single_df <- FALSE
  }
  checkmate::assert_character(split_on, null.ok = TRUE, unique = TRUE)
  if (!is.null(padjust_by)) {
    if (!is.list(padjust_by)) {
      padjust_by <- list(padjust_by)
    }
    sapply(padjust_by, checkmate::assert_character, null.ok = TRUE)
  }
  checkmate::assert_string(padjust_method)
  checkmate::assert_integerish(ncores, lower = 1L)
  checkmate::assert_class(cl, "cluster", null.ok = TRUE)
  checkmate::assert_subset(calculate, c("parameter_estimates_reml", 
                                        "parameter_estimates_ml", "fit_statistics", "residuals", 
                                        "fitted"))
  checkmate::assert_logical(return_models, len = 1L)
  checkmate::assert_character(scale_predictors, null.ok = TRUE)
  checkmate::assert_list(emmeans_spec, null.ok = TRUE)
  checkmate::assert_list(emtrends_spec, null.ok = TRUE)
  if (is.logical(refit_on_nonconvergence) && isFALSE(refit_on_nonconvergence)) {
    refit_on_nonconvergence <- 0L
  }
  checkmate::assert_integerish(refit_on_nonconvergence, null.ok = FALSE)
  validate_form_list <- function(ll) {
    if (inherits(ll, "formula")) {
      ll <- list(ll)
    }
    else {
      checkmate::assert_list(ll)
    }
    lapply(ll, checkmate::assert_formula)
    if (is.null(names(ll))) {
      nm <- paste0("model", seq_along(ll))
      message("Using default model names of ", paste(nm, 
                                                     collapse = ", "))
      names(ll) <- nm
    }
    return(ll)
  }
  if (!is.null(model_formulae)) {
    message("Using model_formulae to setup model_set")
    model_formulae <- validate_form_list(model_formulae)
    model_set <- tibble(outcome = sapply(model_formulae, 
                                         formula.tools::lhs.vars), model_name = names(model_formulae)) %>% 
      mutate(rhs = lapply(model_formulae, function(x) {
        update.formula(x, "NULL ~ .")
      }), form = model_formulae)
  }
  else if (!is.null(rhs_model_formulae)) {
    checkmate::assert_character(outcomes, null.ok = FALSE)
    rhs_model_formulae <- validate_form_list(rhs_model_formulae)
    model_set <- expand.grid(outcome = outcomes, rhs = rhs_model_formulae, 
                             stringsAsFactors = FALSE) %>% rowwise() %>% mutate(form = list(update.formula(rhs, 
                                                                                                           paste(outcome, "~ .")))) %>% ungroup() %>% as_tibble() %>% 
      mutate(model_name = names(.$rhs))
  }
  if (!is.null(external_df)) {
    if (checkmate::test_string(external_df)) {
      checkmate::assert_file_exists(external_df)
      external_df <- readRDS(external_df)
    }
    else {
      checkmate::assert_data_frame(external_df)
      if (!is.data.table(external_df)) {
        data.table::setDT(external_df)
      }
      checkmate::assert_subset(external_merge_by, names(external_df))
      data.table::setkeyv(external_df, external_merge_by)
    }
  }
  model_worker <- function(data, model_formula, lmer_control, 
                           outcome_transform = NULL, scale_predictors = NULL, REML = TRUE) {
    if (!is.null(outcome_transform)) {
      lhs <- all.vars(model_formula)[1]
      data[[lhs]] <- outcome_transform(data[[lhs]])
    }
    if (!is.null(scale_predictors)) {
      model_terms <- all.vars(model_formula)
      if (any(model_terms %in% scale_predictors)) {
        for (aa in intersect(scale_predictors, model_terms)) {
          data[[aa]] <- as.vector(scale(data[[aa]]))
        }
      }
    }
    md <- lmerTest::lmer(model_formula, data, control = lmer_control, 
                         REML = REML, na.action = na.exclude)
    if (refit_on_nonconvergence > 0L) {
      rfc <- 0
      while (any(grepl("failed to converge", md@optinfo$conv$lme4$messages)) && 
             rfc < refit_on_nonconvergence) {
        ss <- getME(md, c("theta", "fixef"))
        lmod <- lmer_control
        lmod$optimizer <- "bobyqa"
        md <- update(md, start = ss, control = lmod)
        rfc <- rfc + 1
      }
    }
    return(md)
  }
  if (isTRUE(single_df)) {
    if (!is.data.table(data)) {
      setDT(data)
    }
    df_set <- c("internal")
  }
  else {
    df_set <- data
    if (!".filename" %in% split_on) {
      split_on <- c(".filename", split_on)
    }
  }
  if (!is.null(emmeans_spec)) {
    if (is.null(names(emmeans_spec))) {
      names(emmeans_spec) <- paste("emm", seq_along(emmeans_spec), 
                                   sep = "_")
    }
    else {
      empty_names <- which(names(emmeans_spec) == "")
      names(emmeans_spec)[empty_names] <- paste("emm", 
                                                empty_names, sep = "_")
    }
    emm_metadata <- rbindlist(lapply(emmeans_spec, function(ee) {
      data.frame(ee[c("outcome", "model_name")])
    })) %>% mutate(emm_label = names(emmeans_spec))
    bad_outcomes <- !emm_metadata$outcome %in% model_set$outcome
    if (any(bad_outcomes)) {
      cat("Cannot have an outcome in emmeans that is not present in the models.\nIgnoring this specification!\n")
      print(emm_metadata[bad_outcomes, ])
    }
    bad_models <- !emm_metadata$model_name %in% model_set$model_name
    if (any(bad_models)) {
      cat("Cannot have model names in emmeans that are not present in the model formula names.\nIgnoring this specification!\n")
      print(emm_metadata[bad_models, ])
    }
    keep_models <- !bad_outcomes & !bad_models
    if (sum(keep_models) == 0L) {
      emmeans_spec <- NULL
    }
    else {
      emm_metadata <- emm_metadata[keep_models, , drop = FALSE]
      emm_metadata$emm_number <- 1:nrow(emm_metadata)
      emmeans_spec <- emmeans_spec[keep_models]
    }
  }
  if (!is.null(emtrends_spec)) {
    if (is.null(names(emtrends_spec))) {
      names(emtrends_spec) <- paste("emt", seq_along(emtrends_spec), 
                                    sep = "_")
    }
    else {
      empty_names <- which(names(emtrends_spec) == "")
      names(emtrends_spec)[empty_names] <- paste("emt", 
                                                 empty_names, sep = "_")
    }
    emt_metadata <- rbindlist(lapply(emtrends_spec, function(ee) {
      data.frame(ee[c("outcome", "model_name")])
    })) %>% mutate(emt_label = names(emtrends_spec))
    bad_outcomes <- !emt_metadata$outcome %in% model_set$outcome
    if (any(bad_outcomes)) {
      cat("Cannot have an outcome in emtrends that is not present in the models.\nIgnoring this specification!\n")
      print(emt_metadata[bad_outcomes, ])
    }
    bad_models <- !emt_metadata$model_name %in% model_set$model_name
    if (any(bad_models)) {
      cat("Cannot have model names in emtrends that are not present in the model formula names.\nIgnoring this specification!\n")
      print(emt_metadata[bad_models, ])
    }
    keep_models <- !bad_outcomes & !bad_models
    if (sum(keep_models) == 0L) {
      emtrends_spec <- NULL
    }
    else {
      emt_metadata <- emt_metadata[keep_models, , drop = FALSE]
      emt_metadata$emt_number <- 1:nrow(emt_metadata)
      emtrends_spec <- emtrends_spec[keep_models]
    }
  }
  if (ncores > 1L && is.null(cl)) {
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    on.exit(try(stopCluster(cl)))
  }
  else if (!is.null(cl)) {
    message("Using external cl object for parallelism. mixed_by will not stop the cluster upon completion.")
    message("Use stopCluster() yourself when done with all computation.")
    registerDoParallel(cl, cores = ncores)
  }
  else {
    registerDoSEQ()
  }
  mresults <- vector(mode = "list", length(df_set))
  for (i in seq_along(df_set)) {
    df_i <- df_set[i]
    if (isFALSE(single_df)) {
      message("Reading file: ", df_i)
      if (grepl(".rds$", df_i, ignore.case = TRUE)) {
        data <- readRDS(df_i)
        if (!is.data.table(data)) {
          setDT(data)
        }
      }
      else if (grepl("(.csv|.csv.gz|.csv.bz2|.dat|.txt|.txt.gz|.txt.bz2)", 
                     df_i, ignore.case = TRUE, perl = TRUE)) {
        data <- data.table::fread(df_i, data.table = TRUE)
      }
      else {
        stop("Unable to sort out this data input: ", 
             df_i)
      }
      data[, `:=`(.filename, basename(df_i))]
    }
    if (is.null(split_on)) {
      split_on <- "split"
      data[, `:=`(split, factor(1))]
      has_split <- FALSE
    }
    else {
      has_split <- TRUE
    }
    if (!is.null(external_df)) {
      checkmate::assert_subset(external_merge_by, names(data))
      data <- merge(data, external_df, by = external_merge_by)
    }
    checkmate::assert_subset(model_set$outcome, names(data))
    checkmate::assert_subset(split_on, names(data))
    setkeyv(data, split_on)
    data <- data[, .(dt = list(.SD)), by = split_on]
    message("Starting processing of data splits")
    mresults[[i]] <- foreach(dt_split = iter(data, by = "row"), 
                             .packages = c("lme4", "lmerTest", "data.table", "dplyr", 
                                           "broom.mixed", "emmeans"), .noexport = "data", 
                             .inorder = FALSE) %dopar% {
                               if (nrow(dt_split$dt[[1L]]) == 0L) {
                                 warning("No rows found in split. Skipping")
                                 return(NULL)
                               }
                               split_results <- lapply(seq_len(nrow(model_set)), 
                                                       function(mm) {
                                                         ff <- model_set$form[[mm]]
                                                         ret <- data.table::copy(dt_split)
                                                         ret[, `:=`(outcome, model_set$outcome[[mm]])]
                                                         ret[, `:=`(model_name, names(model_set$rhs)[mm])]
                                                         ret[, `:=`(rhs, as.character(model_set$rhs[mm]))]
                                                         ret[, `:=`(coef_df_reml, list())]
                                                         ret[, `:=`(coef_df_ml, list())]
                                                         ret[, `:=`(residuals, list())]
                                                         ret[, `:=`(fitted, list())]
                                                         ret[, `:=`(emm, list())]
                                                         vv <- all.vars(ff)
                                                         miss_data <- ret$dt[[1]] %>% dplyr::select(!!vv) %>% 
                                                           mutate(any_miss = rowSums(is.na(dplyr::select(., 
                                                                                                         any_of(!!vv)))) > 0)
                                                         n_present <- miss_data %>% dplyr::filter(any_miss == 
                                                                                                    FALSE) %>% nrow()
                                                         if (n_present == 0) {
                                                           ret[, `:=`(dt, NULL)]
                                                           return(ret)
                                                         }
                                                         if ("parameter_estimates_reml" %in% calculate) {
                                                           thism <- model_worker(ret$dt[[1]], ff, lmer_control, 
                                                                                 outcome_transform, scale_predictors, REML = TRUE)
                                                           ret[, `:=`(coef_df_reml, list(do.call(tidy, 
                                                                                                 append(tidy_args, x = thism))))]
                                                         }
                                                         if (any(c("parameter_estimates_ml", "fit_statistics") %in% 
                                                                 calculate)) {
                                                           thism_ml <- model_worker(ret$dt[[1]], ff, 
                                                                                    lmer_control, outcome_transform, scale_predictors, 
                                                                                    REML = FALSE)
                                                           ret[, `:=`(fit_df, list(glance(thism_ml)))]
                                                           ret[, `:=`(coef_df_ml, list(do.call(tidy, 
                                                                                               append(tidy_args, x = thism_ml))))]
                                                         }
                                                         if ("residuals" %in% calculate) {
                                                           ret[, `:=`(residuals, list(resid = residuals(thism)))]
                                                         }
                                                         if ("fitted" %in% calculate) {
                                                           ret[, `:=`(fitted, list(fitted = fitted(thism)))]
                                                         }
                                                         if (!is.null(emmeans_spec)) {
                                                           emm_torun <- emm_metadata %>% filter(outcome == 
                                                                                                  !!model_set$outcome[[mm]] & model_name == 
                                                                                                  !!names(model_set$rhs)[mm])
                                                           if (nrow(emm_torun) > 0L) {
                                                             this_emmspec <- emmeans_spec[emm_torun %>% 
                                                                                            pull(emm_number)]
                                                             emms <- lapply(seq_along(this_emmspec), 
                                                                            function(emm_i) {
                                                                              tidy(do.call(emmeans, c(this_emmspec[[emm_i]], 
                                                                                                      object = thism))) %>% dplyr::bind_cols(emm_torun[emm_i, 
                                                                                                      ] %>% dplyr::select(emm_number, emm_label))
                                                                            })
                                                             names(emms) <- names(this_emmspec)
                                                             ret[, `:=`(emm, list(list(emms)))]
                                                           }
                                                         }
                                                         if (!is.null(emtrends_spec)) {
                                                           emt_torun <- emt_metadata %>% filter(outcome == 
                                                                                                  !!model_set$outcome[[mm]] & model_name == 
                                                                                                  !!names(model_set$rhs)[mm])
                                                           if (nrow(emt_torun) > 0L) {
                                                             this_emtspec <- emtrends_spec[emt_torun %>% 
                                                                                             pull(emt_number)]
                                                             emts <- lapply(seq_along(this_emtspec), 
                                                                            function(emt_i) {
                                                                              tidy(do.call(emtrends, c(this_emtspec[[emt_i]], 
                                                                                                       object = thism))) %>% dplyr::bind_cols(emt_torun[emt_i, 
                                                                                                       ] %>% dplyr::select(emt_number, emt_label))
                                                                            })
                                                             names(emts) <- names(this_emtspec)
                                                             ret[, `:=`(emt, list(list(emts)))]
                                                           }
                                                         }
                                                         ret[, `:=`(dt, NULL)]
                                                         return(ret)
                                                       })
                               split_results <- rbindlist(split_results, fill = TRUE)
                               coef_df_reml <- NULL
                               if ("parameter_estimates_reml" %in% calculate) {
                                 coef_df_reml <- split_results[, coef_df_reml[[1]], 
                                                               by = .(outcome, model_name, rhs)]
                                 coef_df_reml <- cbind(dt_split[, ..split_on], 
                                                       coef_df_reml)
                               }
                               coef_df_ml <- NULL
                               if ("parameter_estimates_ml" %in% calculate) {
                                 coef_df_ml <- split_results[, coef_df_ml[[1]], 
                                                             by = .(outcome, model_name, rhs)]
                                 coef_df_ml <- cbind(dt_split[, ..split_on], coef_df_ml)
                               }
                               fit_df <- NULL
                               if ("fit_statistics" %in% calculate) {
                                 fit_df <- split_results[, fit_df[[1]], by = .(outcome, 
                                                                               model_name, rhs)]
                                 fit_df <- cbind(dt_split[, ..split_on], fit_df)
                               }
                               emm_data <- NULL
                               if (!is.null(emmeans_spec) && "emm" %in% names(split_results)) {
                                 emm_data <- subset(split_results, sapply(emm, 
                                                                          function(x) !is.null(x)), select = c(split_on, 
                                                                                                               "outcome", "model_name", "rhs", "emm"))
                               }
                               emt_data <- NULL
                               if (!is.null(emtrends_spec) && "emt" %in% names(split_results)) {
                                 emt_data <- subset(split_results, sapply(emt, 
                                                                          function(x) !is.null(x)), select = c(split_on, 
                                                                                                               "outcome", "model_name", "rhs", "emt"))
                               }
                               residual_data <- NULL
                               if ("residuals" %in% calculate) {
                                 residual_data <- split_results[, residuals[[1]], 
                                                                by = .(outcome, model_name, rhs)]
                                 residual_data <- cbind(dt_split[, ..split_on], 
                                                        residual_data)
                               }
                               fitted_data <- NULL
                               if ("fitted" %in% calculate) {
                                 fitted_data <- split_results[, fitted[[1]], by = .(outcome, 
                                                                                    model_name, rhs)]
                                 fitted_data <- cbind(dt_split[, ..split_on], 
                                                      fitted_data)
                               }
                               list(coef_df_reml = coef_df_reml, coef_df_ml = coef_df_ml, 
                                    fit_df = fit_df, emm_data = emm_data, emt_data = emt_data, 
                                    residuals = residual_data, fitted = fitted_data)
                             }
    rm(data)
    gc()
  }
  extract_df <- function(nested_list, element, split_on) {
    browser()
    result_df <- rbindlist(lapply(nested_list, function(df_set) {
      rbindlist(lapply(df_set, "[[", element))
    }))
    if (nrow(result_df) > 0L) {
      data.table::setattr(result_df, "split_on", split_on)
      data.table::setorderv(result_df, split_on)
    }
    return(result_df)
  }
  emm_data <- NULL
  emt_data <- NULL
  coef_results_reml <- NULL
  coef_results_ml <- NULL
  fit_results <- NULL
  fitted_data <- NULL
  residual_data <- NULL
  if ("parameter_estimates_reml" %in% calculate) {
    coef_results_reml <- extract_df(mresults, "coef_df_reml", 
                                    split_on)
  }
  if ("parameter_estimates_ml" %in% calculate) {
    coef_results_ml <- extract_df(mresults, "coef_df_ml", 
                                  split_on)
  }
  if ("fit_statistics" %in% calculate) {
    fit_results <- extract_df(mresults, "fit_df", split_on)
  }
  if ("fitted" %in% calculate) {
    fitted_data <- extract_df(mresults, "fitted", split_on)
  }
  if ("residuals" %in% calculate) {
    residual_data <- extract_df(mresults, "residuals", split_on)
  }
  emmeans_list <- NULL
  if (!is.null(emmeans_spec)) {
    emm_data <- extract_df(mresults, "emm_data", split_on)
    if (nrow(emm_data) > 0L) {
      emmeans_list <- lapply(seq_along(emmeans_spec), function(aa) {
        em_sub <- subset(emm_data, outcome == emmeans_spec[[aa]]$outcome & 
                           model_name == emmeans_spec[[aa]]$model_name)
        nn <- names(em_sub)
        other_keys <- nn[nn != "emm"]
        em_name <- names(emmeans_spec)[aa]
        em_sub[, emm[[1]][[em_name]], by = other_keys]
      })
      names(emmeans_list) <- names(emmeans_spec)
    }
  }
  emtrends_list <- NULL
  if (!is.null(emtrends_spec)) {
    emt_data <- extract_df(mresults, "emt_data", split_on)
    if (nrow(emt_data) > 0L) {
      emtrends_list <- lapply(seq_along(emtrends_spec), 
                              function(aa) {
                                em_sub <- subset(emt_data, outcome == emtrends_spec[[aa]]$outcome & 
                                                   model_name == emtrends_spec[[aa]]$model_name)
                                nn <- names(em_sub)
                                other_keys <- nn[nn != "emt"]
                                em_name <- names(emtrends_spec)[aa]
                                em_sub[, emt[[1]][[em_name]], by = other_keys]
                              })
      names(emtrends_list) <- names(emtrends_spec)
    }
  }
  adjust_dt <- function(dt, padjust_by) {
    for (ff in padjust_by) {
      checkmate::assert_subset(ff, names(dt))
      cname <- paste0("padj_", padjust_method, "_", paste(ff, 
                                                          collapse = "_"))
      dt <- dt[, `:=`((cname), p.adjust(p.value, method = padjust_method)), 
               by = ff]
      if (isTRUE(all.equal(dt[["p.value"]], dt[[cname]]))) {
        warning("p-value adjustment: ", paste(ff, collapse = ", "), 
                " yields the same result as the uncorrected p-value. Setting ", 
                cname, " to NA.")
        dt[[cname]] <- NA_real_
      }
    }
    return(dt)
  }
  if (!is.null(padjust_by)) {
    checkmate::assert_subset(padjust_method, c("holm", "hochberg", 
                                               "hommel", "bonferroni", "BH", "BY", "fdr", "none"))
    if (!is.null(coef_results_reml)) 
      coef_results_reml <- adjust_dt(coef_results_reml, 
                                     padjust_by)
    if (!is.null(coef_results_ml)) 
      coef_results_ml <- adjust_dt(coef_results_ml, padjust_by)
  }
  if (isFALSE(has_split)) {
    if (!is.null(coef_results_reml)) {
      coef_results_reml[, `:=`(split, NULL)]
    }
    if (!is.null(coef_results_ml)) {
      coef_results_ml[, `:=`(split, NULL)]
    }
    if (!is.null(fit_results)) {
      fit_results[, `:=`(split, NULL)]
    }
  }
  return(list(coef_df_reml = coef_results_reml, coef_df_ml = coef_results_ml, 
              fit_df = fit_results, emmeans_list = emmeans_list, emtrends_list = emtrends_list, 
              residuals = residual_data, fitted = fitted_data))
}