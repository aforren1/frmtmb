prepare_predictions_ac <- 
function (bframe, draws, sdata, oos = NULL, nat_cov = FALSE, 
    new = FALSE, ...) 
{
    out <- list()
    nat_cov <- as_one_logical(nat_cov)
    acframe <- subset2(bframe$frame$ac, nat_cov = nat_cov)
    if (!has_rows(acframe)) {
        return(out)
    }
    stopifnot(is.acframe(acframe))
    out$acframe <- acframe
    p <- usc(combine_prefix(bframe))
    out$N_tg <- sdata[[paste0("N_tg", p)]]
    if (has_ac_class(acframe, "arma")) {
        acframe_arma <- subset2(acframe, class = "arma")
        out$Y <- sdata[[paste0("Y", p)]]
        if (!is.null(oos)) {
            if (any(oos > length(out$Y))) {
                stop2("'oos' should not contain integers larger than N.")
            }
            out$Y[oos] <- NA
        }
        out$J_lag <- sdata[[paste0("J_lag", p)]]
        if (acframe_arma$p > 0) {
            ar_regex <- paste0("^ar", p, "\\[")
            out$ar <- prepare_draws(draws, ar_regex, regex = TRUE)
        }
        if (acframe_arma$q > 0) {
            ma_regex <- paste0("^ma", p, "\\[")
            out$ma <- prepare_draws(draws, ma_regex, regex = TRUE)
        }
    }
    if (has_ac_class(acframe, "cosy")) {
        cosy_regex <- paste0("^cosy", p, "$")
        out$cosy <- prepare_draws(draws, cosy_regex, regex = TRUE)
    }
    if (has_ac_class(acframe, "unstr")) {
        cortime_regex <- paste0("^cortime", p, "__")
        out$cortime <- prepare_draws(draws, cortime_regex, regex = TRUE)
        out$Jtime_tg <- sdata[[paste0("Jtime_tg", p)]]
    }
    if (use_ac_cov_time(acframe)) {
        out$begin_tg <- sdata[[paste0("begin_tg", p)]]
        out$end_tg <- sdata[[paste0("end_tg", p)]]
    }
    if (has_ac_latent_residuals(bframe)) {
        err_regex <- paste0("^err", p, "\\[")
        has_err <- any(grepl(err_regex, colnames(draws)))
        if (has_err && !new) {
            out$err <- prepare_draws(draws, err_regex, regex = TRUE)
        }
        else {
            if (!use_ac_cov_time(acframe)) {
                stop2("Cannot predict new latent residuals ", 
                  "when using cov = FALSE in autocor terms.")
            }
            out$err <- matrix(nrow = nrow(draws), ncol = length(out$Y))
            sderr_regex <- paste0("^sderr", p, "$")
            out$sderr <- prepare_draws(draws, sderr_regex, regex = TRUE)
            for (i in seq_len(out$N_tg)) {
                obs <- with(out, begin_tg[i]:end_tg[i])
                Jtime <- out$Jtime_tg[i, ]
                cov <- get_cov_matrix_ac(list(ac = out), obs, 
                  Jtime = Jtime, latent = TRUE)
                zeros <- rep(0, length(obs))
                .err <- function(s) rmulti_normal(1, zeros, Sigma = cov[s, 
                  , ])
                out$err[, obs] <- rblapply(seq_rows(draws), .err)
            }
        }
    }
    if (has_ac_class(acframe, "sar")) {
        lagsar_regex <- paste0("^lagsar", p, "$")
        errorsar_regex <- paste0("^errorsar", p, "$")
        out$lagsar <- prepare_draws(draws, lagsar_regex, regex = TRUE)
        out$errorsar <- prepare_draws(draws, errorsar_regex, 
            regex = TRUE)
        out$Msar <- sdata[[paste0("Msar", p)]]
    }
    if (has_ac_class(acframe, "car")) {
        acframe_car <- subset2(acframe, class = "car")
        if (new && acframe_car$gr == "NA") {
            stop2("Without a grouping factor, CAR models cannot handle newdata.")
        }
        gcar <- sdata[[paste0("Jloc", p)]]
        Zcar <- matrix(rep(1, length(gcar)))
        out$Zcar <- prepare_Z(Zcar, list(gcar))
        rcar_regex <- paste0("^rcar", p, "\\[")
        rcar <- prepare_draws(draws, rcar_regex, regex = TRUE)
        rcar <- rcar[, unique(gcar), drop = FALSE]
        out$rcar <- rcar
    }
    if (has_ac_class(acframe, "fcor")) {
        out$Mfcor <- sdata[[paste0("Mfcor", p)]]
    }
    out
}
