prepare_predictions_data <- 
function (bframe, sdata, stanvars = NULL, ...) 
{
    resp <- usc(combine_prefix(bframe))
    vars <- c("Y", "trials", "ncat", "nthres", "se", "weights", 
        "denom", "dec", "cens", "rcens", "lb", "ub")
    vars <- paste0(vars, resp)
    vars <- intersect(vars, names(sdata))
    escaped_resp <- escape_all(resp)
    vl_vars <- c("vreal", "vint")
    vl_vars <- regex_or(vl_vars)
    vl_vars <- paste0("^", vl_vars, "[[:digit:]]+", escaped_resp, 
        "$")
    vl_vars <- str_subset(names(sdata), vl_vars)
    vars <- union(vars, vl_vars)
    out <- sdata[vars]
    names(out) <- sub(paste0(escaped_resp, "$"), "", names(out))
    if (length(stanvars)) {
        stopifnot(is.stanvars(stanvars))
        out[names(stanvars)] <- sdata[names(stanvars)]
    }
    out
}
