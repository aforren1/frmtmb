prepare_family <- 
function (x) 
{
    stopifnot(is.brmsformula(x) || is.brmsterms(x))
    family <- x$family
    acframe <- frame_ac(x)
    family$fun <- family[["fun"]] %||% family$family
    if (use_ac_cov_time(acframe) && has_natural_residuals(x)) {
        family$fun <- paste0(family$fun, "_time")
    }
    else if (has_ac_class(acframe, "sar")) {
        acframe_sar <- subset2(acframe, class = "sar")
        if (has_ac_subset(acframe_sar, type = "lag")) {
            family$fun <- paste0(family$fun, "_lagsar")
        }
        else if (has_ac_subset(acframe_sar, type = "error")) {
            family$fun <- paste0(family$fun, "_errorsar")
        }
    }
    else if (has_ac_class(acframe, "fcor")) {
        family$fun <- paste0(family$fun, "_fcor")
    }
    family
}
