data_ac <- 
function (bframe, data, data2, ...) 
{
    if (!is.null(bframe$sdata$ac)) {
        return(bframe$sdata$ac)
    }
    out <- list()
    N <- nrow(data)
    basis <- bframe$basis$ac
    acframe <- bframe$frame$ac
    stopifnot(is.acframe(acframe))
    if (has_ac_subset(bframe, dim = "time")) {
        gr <- get_ac_vars(acframe, "gr", dim = "time")
        if (isTRUE(nzchar(gr))) {
            tgroup <- as.numeric(factor(data[[gr]]))
        }
        else {
            tgroup <- rep(1, N)
        }
    }
    if (has_ac_class(acframe, "arma")) {
        acframe_arma <- subset2(acframe, class = "arma")
        out$Kar <- acframe_arma$p
        out$Kma <- acframe_arma$q
        if (!use_ac_cov_time(acframe_arma)) {
            max_lag <- max(out$Kar, out$Kma)
            out$J_lag <- as.array(rep(0, N))
            for (n in seq_len(N)[-N]) {
                ind <- n:max(1, n + 1 - max_lag)
                out$J_lag[n] <- sum(tgroup[ind] %in% tgroup[n + 
                  1])
            }
        }
    }
    if (use_ac_cov_time(acframe)) {
        out$N_tg <- length(unique(tgroup))
        out$begin_tg <- as.array(ulapply(unique(tgroup), match, 
            tgroup))
        out$nobs_tg <- as.array(with(out, c(if (N_tg > 1) begin_tg[2:N_tg], 
            N + 1) - begin_tg))
        out$end_tg <- with(out, begin_tg + nobs_tg - 1)
        if (has_ac_class(acframe, "unstr")) {
            time <- get_ac_vars(bframe, "time", dim = "time")
            time_data <- get(time, data)
            new_times <- extract_levels(time_data)
            if (length(basis)) {
                times <- basis$times
                invalid_times <- setdiff(new_times, times)
                if (length(invalid_times)) {
                  stop2("Cannot handle new time points in UNSTR models.")
                }
            }
            else {
                times <- new_times
            }
            out$n_unique_t <- length(times)
            out$n_unique_cortime <- out$n_unique_t * (out$n_unique_t - 
                1)/2
            Jtime <- match(time_data, times)
            out$Jtime_tg <- matrix(0, out$N_tg, max(out$nobs_tg))
            for (i in seq_len(out$N_tg)) {
                out$Jtime_tg[i, seq_len(out$nobs_tg[i])] <- Jtime[out$begin_tg[i]:out$end_tg[i]]
            }
        }
    }
    if (has_ac_class(acframe, "sar")) {
        acframe_sar <- subset2(acframe, class = "sar")
        M <- data2[[acframe_sar$M]]
        rmd_rows <- attr(data, "na.action")
        if (!is.null(rmd_rows)) {
            class(rmd_rows) <- NULL
            M <- M[-rmd_rows, -rmd_rows, drop = FALSE]
        }
        if (!is_equal(dim(M), rep(N, 2))) {
            stop2("Dimensions of 'M' for SAR terms must be equal to ", 
                "the number of observations.")
        }
        out$Msar <- as.matrix(M)
        out$eigenMsar <- eigen(M)$values
        out$N_tg <- 1
    }
    if (has_ac_class(acframe, "car")) {
        acframe_car <- subset2(acframe, class = "car")
        locations <- NULL
        if (length(basis)) {
            locations <- basis$locations
        }
        M <- data2[[acframe_car$M]]
        if (acframe_car$gr != "NA") {
            loc_data <- get(acframe_car$gr, data)
            new_locations <- extract_levels(loc_data)
            if (is.null(locations)) {
                locations <- new_locations
            }
            else {
                invalid_locations <- setdiff(new_locations, locations)
                if (length(invalid_locations)) {
                  stop2("Cannot handle new locations in CAR models.")
                }
            }
            Nloc <- length(locations)
            Jloc <- as.array(match(loc_data, locations))
            if (is.null(rownames(M))) {
                stop2("Row names are required for 'M' in CAR terms.")
            }
            found <- locations %in% rownames(M)
            if (any(!found)) {
                stop2("Row names of 'M' for CAR terms do not match ", 
                  "the names of the grouping levels.")
            }
            M <- M[locations, locations, drop = FALSE]
        }
        else {
            warning2("Using CAR terms without a grouping factor is deprecated. ", 
                "Please use argument 'gr' even if each observation ", 
                "represents its own location.")
            Nloc <- N
            Jloc <- as.array(seq_len(Nloc))
            if (!is_equal(dim(M), rep(Nloc, 2))) {
                if (length(basis)) {
                  stop2("Cannot handle new data in CAR terms ", 
                    "without a grouping factor.")
                }
                else {
                  stop2("Dimensions of 'M' for CAR terms must be equal ", 
                    "to the number of observations.")
                }
            }
        }
        edges_rows <- (Matrix::tril(M)@i + 1)
        edges_cols <- sort(Matrix::triu(M)@i + 1)
        edges <- cbind(rows = edges_rows, cols = edges_cols)
        c(out) <- nlist(Nloc, Jloc, Nedges = length(edges_rows), 
            edges1 = as.array(edges_rows), edges2 = as.array(edges_cols))
        if (acframe_car$type %in% c("escar", "esicar")) {
            Nneigh <- Matrix::colSums(M)
            if (any(Nneigh == 0) && !length(basis)) {
                stop2("For exact sparse CAR, all locations should have at ", 
                  "least one neighbor within the provided data set. ", 
                  "Consider using type = 'icar' instead.")
            }
            inv_sqrt_D <- diag(1/sqrt(Nneigh))
            eigenMcar <- t(inv_sqrt_D) %*% M %*% inv_sqrt_D
            eigenMcar <- eigen(eigenMcar, TRUE, only.values = TRUE)$values
            c(out) <- nlist(Nneigh, eigenMcar)
        }
        else if (acframe_car$type %in% "bym2") {
            c(out) <- list(car_scale = .car_scale(edges, Nloc))
        }
    }
    if (has_ac_class(acframe, "fcor")) {
        acframe_fcor <- subset2(acframe, class = "fcor")
        M <- data2[[acframe_fcor$M]]
        rmd_rows <- attr(data, "na.action")
        if (!is.null(rmd_rows)) {
            class(rmd_rows) <- NULL
            M <- M[-rmd_rows, -rmd_rows, drop = FALSE]
        }
        if (nrow(M) != N) {
            stop2("Dimensions of 'M' for FCOR terms must be equal ", 
                "to the number of observations.")
        }
        out$Mfcor <- M
        out$N_tg <- 1
    }
    if (length(out)) {
        resp <- usc(combine_prefix(bframe))
        out <- setNames(out, paste0(names(out), resp))
    }
    out
}
