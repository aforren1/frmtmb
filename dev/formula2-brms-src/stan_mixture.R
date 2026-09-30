function (bterms, prior, threads, normalize, ...) 
{
    out <- list()
    if (!is.mixfamily(bterms$family)) {
        return(out)
    }
    lpdf <- stan_lpdf_name(normalize)
    px <- check_prefix(bterms)
    p <- usc(combine_prefix(px))
    nmix <- length(bterms$family$mix)
    theta_pred <- grepl("^theta", names(bterms$dpars))
    theta_pred <- bterms$dpars[theta_pred]
    theta_fix <- grepl("^theta", names(bterms$fdpars))
    theta_fix <- bterms$fdpars[theta_fix]
    def_thetas <- cglue("  real<lower=0,upper=1> theta{1:nmix}{p};  // mixing proportion\n")
    if (length(theta_pred)) {
        if (length(theta_pred) != nmix - 1) {
            stop2("Can only predict all but one mixing proportion.")
        }
        missing_id <- setdiff(1:nmix, dpar_id(names(theta_pred)))
        str_add(out$model_def) <- glue("  vector[N{p}] theta{missing_id}{p} = rep_vector(0.0, N{p});\n", 
            "  real log_sum_exp_theta{p};\n")
        sum_exp_theta <- glue("exp(theta{1:nmix}{p}[n])", collapse = " + ")
        str_add(out$model_comp_mix) <- glue("  for (n in 1:N{p}) {{\n", 
            "    // scale theta to become a probability vector\n", 
            "    log_sum_exp_theta{p} = log({sum_exp_theta});\n")
        str_add(out$model_comp_mix) <- cglue("    theta{1:nmix}{p}[n] = theta{1:nmix}{p}[n] - log_sum_exp_theta{p};\n")
        str_add(out$model_comp_mix) <- "  }\n"
    }
    else if (length(theta_fix)) {
        if (length(theta_fix) != nmix) {
            stop2("Can only fix no or all mixing proportions.")
        }
        str_add(out$data) <- "  // mixing proportions\n"
        str_add(out$data) <- cglue("  real<lower=0,upper=1> theta{1:nmix}{p};\n")
        str_add(out$pll_args) <- cglue(", real theta{1:nmix}{p}")
    }
    else {
        str_add(out$data) <- glue("  vector[{nmix}] con_theta{p};  // prior concentration\n")
        str_add(out$par) <- glue("  simplex[{nmix}] theta{p};  // mixing proportions\n")
        str_add(out$tpar_prior) <- glue("  lprior += dirichlet_{lpdf}(theta{p} | con_theta{p});\n")
        str_add(out$tpar_def) <- "  // mixing proportions\n"
        str_add(out$tpar_def) <- cglue("  real<lower=0,upper=1> theta{1:nmix}{p};\n")
        str_add(out$tpar_comp) <- cglue("  theta{1:nmix}{p} = theta{p}[{1:nmix}];\n")
        str_add(out$pll_args) <- cglue(", real theta{1:nmix}{p}")
    }
    if (order_intercepts(bterms)) {
        str_add(out$par) <- glue("  ordered[{nmix}] ordered_Intercept{p};  // to identify mixtures\n")
    }
    if (fix_intercepts(bterms)) {
        stopifnot(is_ordinal(bterms))
        gr <- grb <- ""
        groups <- get_thres_groups(bterms)
        if (has_thres_groups(bterms)) {
            gr <- usc(seq_along(groups))
            grb <- paste0("[", seq_along(groups), "]")
        }
        type <- str_if(has_ordered_thres(bterms), "ordered", 
            "vector")
        coef_type <- str_if(has_ordered_thres(bterms), "", "real")
        for (i in seq_along(groups)) {
            str_add_list(out) <- stan_prior(prior, class = "Intercept", 
                coef = get_thres(bterms, group = groups[i]), 
                type = glue("{type}[nthres{p}{grb[i]}]"), coef_type = coef_type, 
                px = px, prefix = "fixed_", suffix = glue("{p}{gr[i]}"), 
                comment = "thresholds fixed over mixture components", 
                normalize = normalize)
        }
    }
    out
}
