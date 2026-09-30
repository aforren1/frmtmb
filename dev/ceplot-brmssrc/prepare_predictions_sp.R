prepare_predictions_sp <- 
function (bframe, draws, sdata, new = FALSE, ...) 
{
    stopifnot(is.bframel(bframe))
    out <- list()
    spframe <- bframe$frame$sp
    meframe <- bframe$frame$me
    if (!has_rows(spframe)) {
        return(out)
    }
    p <- usc(combine_prefix(bframe))
    resp <- usc(bframe$resp)
    out$calls <- vector("list", nrow(spframe))
    for (i in seq_along(out$calls)) {
        call <- spframe$joint_call[[i]]
        if (!is.null(spframe$calls_mo[[i]])) {
            new_mo <- paste0(".mo(simo_", spframe$Imo[[i]], ", Xmo_", 
                spframe$Imo[[i]], ")")
            call <- rename(call, spframe$calls_mo[[i]], new_mo)
        }
        if (!is.null(spframe$calls_me[[i]])) {
            new_me <- paste0("Xme_", seq_along(meframe$term))
            call <- rename(call, meframe$term, new_me)
        }
        if (!is.null(spframe$calls_mi[[i]])) {
            is_na_idx <- is.na(spframe$idx2_mi[[i]])
            idx_mi <- paste0("idxl", p, "_", spframe$vars_mi[[i]], 
                "_", spframe$idx2_mi[[i]])
            idx_mi <- ifelse(is_na_idx, "", paste0("[, ", idx_mi, 
                "]"))
            new_mi <- paste0("Yl_", spframe$vars_mi[[i]], idx_mi)
            call <- rename(call, spframe$calls_mi[[i]], new_mi)
        }
        if (spframe$Ic[i] > 0) {
            str_add(call) <- paste0(" * Csp_", spframe$Ic[i])
        }
        out$calls[[i]] <- parse(text = paste0(call))
    }
    bsp_pars <- paste0("bsp", p, "_", spframe$coef)
    out$bsp <- prepare_draws(draws, bsp_pars)
    colnames(out$bsp) <- spframe$coef
    simo_coef <- get_simo_labels(spframe)
    Jmo <- sdata[[paste0("Jmo", p)]]
    out$simo <- out$Xmo <- named_list(simo_coef)
    for (i in seq_along(simo_coef)) {
        J <- seq_len(Jmo[i])
        simo_par <- paste0("simo", p, "_", simo_coef[i], "[", 
            J, "]")
        out$simo[[i]] <- prepare_draws(draws, simo_par)
        out$Xmo[[i]] <- sdata[[paste0("Xmo", p, "_", i)]]
    }
    warn_me <- FALSE
    if (has_rows(meframe)) {
        save_mevars <- any(grepl("^Xme_", colnames(draws)))
        warn_me <- warn_me || !new && !save_mevars
        out$Xme <- named_list(meframe$coef)
        Xme_regex <- paste0("^Xme_", escape_all(meframe$coef), 
            "\\[")
        Xn <- sdata[paste0("Xn_", seq_rows(meframe))]
        noise <- sdata[paste0("noise_", seq_rows(meframe))]
        groups <- unique(meframe$grname)
        for (i in seq_along(groups)) {
            g <- groups[i]
            K <- which(meframe$grname %in% g)
            if (nzchar(g)) {
                Jme <- sdata[[paste0("Jme_", i)]]
            }
            if (!new && save_mevars) {
                for (k in K) {
                  out$Xme[[k]] <- prepare_draws(draws, Xme_regex[k], 
                    regex = TRUE)
                }
            }
            else {
                if (nzchar(g)) {
                  Jme <- as.numeric(factor(Jme))
                  me_dim <- c(nrow(out$bsp), max(Jme))
                }
                else {
                  me_dim <- c(nrow(out$bsp), sdata$N)
                }
                for (k in K) {
                  dXn <- data2draws(Xn[[k]], me_dim)
                  dnoise <- data2draws(noise[[k]], me_dim)
                  out$Xme[[k]] <- array(rnorm(prod(me_dim), dXn, 
                    dnoise), me_dim)
                  remove(dXn, dnoise)
                }
            }
            if (nzchar(g)) {
                for (k in K) {
                  out$Xme[[k]] <- out$Xme[[k]][, Jme, drop = FALSE]
                }
            }
        }
    }
    dim <- c(nrow(out$bsp), sdata[[paste0("N", resp)]])
    vars_mi <- unique(unlist(spframe$vars_mi))
    if (length(vars_mi)) {
        Yl_names <- paste0("Yl_", vars_mi)
        out$Yl <- named_list(Yl_names)
        for (i in seq_along(out$Yl)) {
            vmi <- vars_mi[i]
            dim_y <- c(nrow(out$bsp), sdata[[paste0("N_", vmi)]])
            Y <- data2draws(sdata[[paste0("Y_", vmi)]], dim_y)
            sdy <- sdata[[paste0("noise_", vmi)]]
            if (is.null(sdy)) {
                out$Yl[[i]] <- Y
                if (!new) {
                  Ymi_regex <- paste0("^Ymi_", escape_all(vmi), 
                    "\\[")
                  Ymi <- prepare_draws(draws, Ymi_regex, regex = TRUE)
                  Jmi <- sdata[[paste0("Jmi_", vmi)]]
                  out$Yl[[i]][, Jmi] <- Ymi
                }
            }
            else {
                save_mevars <- any(grepl("^Yl_", colnames(draws)))
                if (save_mevars && !new) {
                  Ymi_regex <- paste0("^Yl_", escape_all(vmi), 
                    "\\[")
                  out$Yl[[i]] <- prepare_draws(draws, Ymi_regex, 
                    regex = TRUE)
                }
                else {
                  warn_me <- warn_me || !new
                  sdy <- data2draws(sdy, dim)
                  out$Yl[[i]] <- rcontinuous(n = prod(dim), dist = "norm", 
                    mean = Y, sd = sdy, lb = sdata[[paste0("lbmi_", 
                      vmi)]], ub = sdata[[paste0("ubmi_", vmi)]])
                  out$Yl[[i]] <- array(out$Yl[[i]], dim_y)
                }
            }
        }
        uni_mi <- na.omit(attr(spframe, "uni_mi"))
        idxl_vars <- paste0("idxl", p, "_", uni_mi$var, "_", 
            uni_mi$idx2)
        out$idxl <- sdata[idxl_vars]
    }
    if (warn_me) {
        warning2("Noise-free latent variables were not saved. ", 
            "You can control saving those variables via 'save_pars()'. ", 
            "Treating original data as if it was new data as a workaround.")
    }
    ncovars <- max(spframe$Ic)
    out$Csp <- vector("list", ncovars)
    for (i in seq_len(ncovars)) {
        out$Csp[[i]] <- sdata[[paste0("Csp", p, "_", i)]]
        out$Csp[[i]] <- data2draws(out$Csp[[i]], dim = dim)
    }
    out
}
