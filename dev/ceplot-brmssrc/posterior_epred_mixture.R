posterior_epred_mixture <- 
function (prep) 
{
    families <- family_names(prep$family)
    prep$dpars$theta <- get_theta(prep)
    out <- 0
    for (j in seq_along(families)) {
        posterior_epred_fun <- paste0("posterior_epred_", families[j])
        posterior_epred_fun <- get(posterior_epred_fun, asNamespace("brms"))
        tmp_prep <- pseudo_prep_for_mixture(prep, j)
        if (length(dim(prep$dpars$theta)) == 3) {
            theta <- prep$dpars$theta[, , j]
        }
        else {
            theta <- prep$dpars$theta[, j]
        }
        out <- out + theta * posterior_epred_fun(tmp_prep)
    }
    out
}
