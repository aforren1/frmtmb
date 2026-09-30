posterior_epred_trunc <- 
function (prep) 
{
    stopifnot(is_trunc(prep))
    lb <- data2draws(prep$data[["lb"]], dim_mu(prep))
    ub <- data2draws(prep$data[["ub"]], dim_mu(prep))
    posterior_epred_trunc_fun <- paste0("posterior_epred_trunc_", 
        prep$family$family)
    posterior_epred_trunc_fun <- try(get(posterior_epred_trunc_fun, 
        asNamespace("brms")), silent = TRUE)
    if (is_try_error(posterior_epred_trunc_fun)) {
        stop2("posterior_epred values on the respone scale not yet implemented ", 
            "for truncated '", prep$family$family, "' models.")
    }
    trunc_args <- nlist(prep, lb, ub)
    do_call(posterior_epred_trunc_fun, trunc_args)
}
