function (bframe, prior, stanvars, threads, primitive, normalize, 
    ...) 
{
    stopifnot(is.bframel(bframe))
    out <- list()
    family <- bframe$family
    fixef <- bframe$frame$fe$vars_stan
    sparse <- bframe$frame$fe$sparse
    decomp <- bframe$frame$fe$decomp
    center <- bframe$frame$fe$center
    ct <- str_if(center, "c")
    px <- check_prefix(bframe)
    p <- usc(combine_prefix(px))
    resp <- usc(px$resp)
    lpdf <- stan_lpdf_name(normalize)
    if (length(fixef)) {
        str_add(out$data) <- glue("  int<lower=1> K{p};", "  // number of population-level effects\n", 
            "  matrix[N{resp}, K{p}] X{p};", "  // population-level design matrix\n")
        if (decomp == "none") {
            str_add(out$pll_args) <- glue(", data matrix X{ct}{p}")
        }
        if (sparse) {
            if (decomp != "none") {
                stop2("Cannot use ", decomp, " decomposition for sparse matrices.")
            }
            if (use_threading(threads)) {
                stop2("Cannot use threading and sparse matrices at the same time.")
            }
            str_add(out$tdata_def) <- glue("  // sparse matrix representation of X{p}\n", 
                "  vector[rows(csr_extract_w(X{p}))] wX{p}", 
                " = csr_extract_w(X{p});\n", "  int vX{p}[size(csr_extract_v(X{p}))]", 
                " = csr_extract_v(X{p});\n", "  int uX{p}[size(csr_extract_u(X{p}))]", 
                " = csr_extract_u(X{p});\n")
        }
        b_type <- glue("vector[K{ct}{p}]")
        has_special_prior <- has_special_prior(prior, bframe, 
            class = "b")
        if (decomp == "none") {
            if (has_special_prior) {
                str_add_list(out) <- stan_prior_non_centered(suffix = p, 
                  suffix_K = ct, normalize = normalize)
            }
            else {
                str_add_list(out) <- stan_prior(prior, class = "b", 
                  coef = fixef, type = b_type, px = px, suffix = p, 
                  header_type = "vector", comment = "regression coefficients", 
                  normalize = normalize)
            }
        }
        else {
            stopifnot(decomp == "QR")
            stopif_prior_bound(prior, class = "b", ls = px)
            if (has_special_prior) {
                str_add_list(out) <- stan_prior_non_centered(suffix = p, 
                  suffix_class = "Q", suffix_K = ct, normalize = normalize)
            }
            else {
                str_add_list(out) <- stan_prior(prior, class = "b", 
                  coef = fixef, type = b_type, px = px, suffix = glue("Q{p}"), 
                  header_type = "vector", comment = "regression coefficients on QR scale", 
                  normalize = normalize)
            }
            str_add(out$gen_def) <- glue("  // obtain the actual coefficients\n", 
                "  vector[K{ct}{p}] b{p} = XR{p}_inv * bQ{p};\n")
        }
    }
    order_intercepts <- order_intercepts(bframe)
    if (order_intercepts && !center) {
        stop2("Identifying mixture components via ordering requires ", 
            "population-level intercepts to be present.\n", "Try setting order = 'none' in function 'mixture'.")
    }
    if (center) {
        sub_X_means <- ""
        if (length(fixef)) {
            str_add(out$data) <- glue("  int<lower=1> Kc{p};", 
                "  // number of population-level effects after centering\n")
            sub_X_means <- glue(" - dot_product(means_X{p}, b{p})")
            if (is_ordinal(family)) {
                str_add(out$tdata_def) <- glue("  matrix[N{resp}, Kc{p}] Xc{p};", 
                  "  // centered version of X{p}\n", "  vector[Kc{p}] means_X{p};", 
                  "  // column means of X{p} before centering\n")
                str_add(out$tdata_comp) <- glue("  for (i in 1:K{p}) {{\n", 
                  "    means_X{p}[i] = mean(X{p}[, i]);\n", "    Xc{p}[, i] = X{p}[, i] - means_X{p}[i];\n", 
                  "  }}\n")
            }
            else {
                str_add(out$tdata_def) <- glue("  matrix[N{resp}, Kc{p}] Xc{p};", 
                  "  // centered version of X{p} without an intercept\n", 
                  "  vector[Kc{p}] means_X{p};", "  // column means of X{p} before centering\n")
                str_add(out$tdata_comp) <- glue("  for (i in 2:K{p}) {{\n", 
                  "    means_X{p}[i - 1] = mean(X{p}[, i]);\n", 
                  "    Xc{p}[, i - 1] = X{p}[, i] - means_X{p}[i - 1];\n", 
                  "  }}\n")
            }
        }
        if (!is_ordinal(family)) {
            intercept_type <- "real"
            if (order_intercepts) {
                dp_id <- dpar_id(px$dpar)
                str_add(out$tpar_def) <- glue("  // identify mixtures via ordering of the intercepts\n", 
                  "  real Intercept{p} = ordered_Intercept{resp}[{dp_id}];\n")
                str_add(out$pll_args) <- glue(", real Intercept{p}")
                intercept_type <- ""
            }
            str_add(out$eta) <- glue(" + Intercept{p}")
            str_add(out$gen_def) <- glue("  // actual population-level intercept\n", 
                "  real b{p}_Intercept = Intercept{p}{sub_X_means};\n")
            str_add_list(out) <- stan_prior(prior, class = "Intercept", 
                type = intercept_type, suffix = p, px = px, header_type = "real", 
                comment = "temporary intercept for centered predictors", 
                normalize = normalize)
        }
    }
    if (decomp == "QR") {
        if (!length(fixef)) {
            stop2("QR decomposition requires non-intercept predictors.")
        }
        str_add(out$tdata_def) <- glue("  // matrices for QR decomposition\n", 
            "  matrix[N{resp}, K{ct}{p}] XQ{p};\n", "  matrix[K{ct}{p}, K{ct}{p}] XR{p};\n", 
            "  matrix[K{ct}{p}, K{ct}{p}] XR{p}_inv;\n")
        str_add(out$tdata_comp) <- glue("  // compute and scale QR decomposition\n", 
            "  XQ{p} = qr_thin_Q(X{ct}{p}) * sqrt(N{resp} - 1);\n", 
            "  XR{p} = qr_thin_R(X{ct}{p}) / sqrt(N{resp} - 1);\n", 
            "  XR{p}_inv = inverse(XR{p});\n")
        str_add(out$pll_args) <- glue(", data matrix XQ{p}")
    }
    if (length(fixef) && !primitive) {
        if (sparse) {
            stopifnot(!center && decomp == "none")
            csr_args <- sargs(paste0(c("rows", "cols"), "(X", 
                p, ")"), paste0(c("wX", "vX", "uX", "b"), p))
            eta_fe <- glue(" + csr_matrix_times_vector({csr_args})")
        }
        else {
            sfx_X <- sfx_b <- ""
            if (decomp == "QR") {
                sfx_X <- sfx_b <- "Q"
            }
            else if (center) {
                sfx_X <- "c"
            }
            slice <- stan_slice(threads)
            eta_fe <- glue(" + X{sfx_X}{p}{slice} * b{sfx_b}{p}")
        }
        str_add(out$eta) <- eta_fe
    }
    out
}
