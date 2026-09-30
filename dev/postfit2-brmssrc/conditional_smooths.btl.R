conditional_smooths.btl <- 
function (x, fit, smooths, conditions, int_conditions, probs, 
    resolution, too_far, spaghetti, surface, ...) 
{
    stopifnot(is.brmsfit(fit))
    out <- list()
    mf <- model.frame(fit)
    smframe <- frame_sm(x, mf)
    smframe$term <- rm_wsp(smframe$term)
    smterms <- unique(smframe$term)
    if (!length(smooths)) {
        I <- seq_along(smterms)
    }
    else {
        I <- which(smterms %in% smooths)
    }
    for (i in I) {
        smooth <- smterms[i]
        sub_smframe <- subset2(smframe, term = smooth)
        covars <- all_vars(sub_smframe$covars[[1]])
        byvars <- all_vars(sub_smframe$byvars[[1]])
        if (length(covars) > 2) {
            byvars <- c(covars[-(1:2)], byvars)
            covars <- covars[1:2]
        }
        else if (length(covars) == 1 && length(byvars)) {
            covars <- c(covars, byvars[1])
            byvars <- byvars[-1]
        }
        vars <- c(covars, byvars)
        values <- named_list(vars)
        is_numeric <- setNames(rep(FALSE, length(covars)), covars)
        for (cv in covars) {
            is_numeric[cv] <- is.numeric(mf[[cv]])
            is_second_covar <- isTRUE(cv == covars[2])
            if (cv %in% names(int_conditions)) {
                int_cond <- int_conditions[[cv]]
                if (is.function(int_cond)) {
                  int_cond <- int_cond(mf[[cv]])
                }
                values[[cv]] <- int_cond
            }
            else if (is_numeric[cv]) {
                if (!surface && is_second_covar) {
                  mean2 <- mean(mf[[cv]], na.rm = TRUE)
                  sd2 <- sd(mf[[cv]], na.rm = TRUE)
                  values[[cv]] <- (-1:1) * sd2 + mean2
                }
                else {
                  values[[cv]] <- seq(min(mf[[cv]]), max(mf[[cv]]), 
                    length.out = resolution)
                }
            }
            else {
                values[[cv]] <- levels(factor(mf[[cv]]))
            }
        }
        for (cv in byvars) {
            if (cv %in% names(int_conditions)) {
                int_cond <- int_conditions[[cv]]
                if (is.function(int_cond)) {
                  int_cond <- int_cond(mf[[cv]])
                }
                values[[cv]] <- int_cond
            }
            else if (is.numeric(mf[[cv]])) {
                mean2 <- mean(mf[[cv]], na.rm = TRUE)
                sd2 <- sd(mf[[cv]], na.rm = TRUE)
                values[[cv]] <- (-1:1) * sd2 + mean2
            }
            else {
                values[[cv]] <- levels(factor(mf[[cv]]))
            }
        }
        newdata <- expand.grid(values)
        show_surface <- surface && length(covars) == 2 && all(is_numeric)
        if (show_surface && too_far > 0) {
            ex_too_far <- mgcv::exclude.too.far(g1 = newdata[[covars[1]]], 
                g2 = newdata[[covars[2]]], d1 = mf[, covars[1]], 
                d2 = mf[, covars[2]], dist = too_far)
            newdata <- newdata[!ex_too_far, ]
        }
        other_vars <- setdiff(names(conditions), vars)
        newdata <- fill_newdata(newdata, other_vars, conditions)
        eta <- posterior_smooths(x, fit, smooth, newdata, ...)
        effects <- na.omit(covars[1:2])
        cond_data <- add_effects__(newdata[, vars, drop = FALSE], 
            effects)
        second_numeric <- isTRUE(is_numeric[2])
        if (second_numeric && !surface) {
            mde2 <- round(cond_data[[effects[2]]], 2)
            levels2 <- sort(unique(mde2), TRUE)
            cond_data$effect2__ <- factor(mde2, levels = levels2)
            labels2 <- names(int_conditions[[effects[2]]])
            if (length(labels2) == length(levels2)) {
                levels(cond_data$effect2__) <- labels2
            }
        }
        if (length(byvars)) {
            cond_data$cond__ <- rows2labels(cond_data[, byvars, 
                drop = FALSE])
        }
        else {
            cond_data$cond__ <- factor(1)
        }
        spa_data <- NULL
        if (spaghetti && !show_surface && is_numeric[1]) {
            sample <- rep(seq_rows(eta), each = ncol(eta))
            if (length(covars) == 2) {
                sample <- paste0(sample, "_", cond_data[[effects[2]]])
            }
            spa_data <- data.frame(as.numeric(t(eta)), factor(sample))
            colnames(spa_data) <- c("estimate__", "sample__")
            spa_data <- cbind(cond_data, spa_data)
        }
        eta <- posterior_summary(eta, robust = TRUE, probs = probs)
        colnames(eta) <- c("estimate__", "se__", "lower__", "upper__")
        eta <- cbind(cond_data, eta)
        response <- combine_prefix(x, keep_mu = TRUE)
        response <- paste0(response, ": ", smooth)
        points <- mf[, vars, drop = FALSE]
        points <- add_effects__(points, covars)
        attr(eta, "response") <- response
        attr(eta, "effects") <- effects
        attr(eta, "surface") <- show_surface
        attr(eta, "spaghetti") <- spa_data
        attr(eta, "points") <- points
        out[[response]] <- eta
    }
    out
}
