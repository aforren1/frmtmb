validate_data2 <- 
function (data2, bterms, ...) 
{
    if (is.null(data2)) {
        data2 <- list()
    }
    if (!is.list(data2)) {
        stop2("'data2' must be a list.")
    }
    if (length(data2) && !is_named(data2)) {
        stop2("All elements of 'data2' must be named.")
    }
    dots <- list(...)
    for (i in seq_along(dots)) {
        if (length(dots[[i]])) {
            stopifnot(is.list(dots[[i]]), is_named(dots[[i]]))
            data2[names(dots[[i]])] <- dots[[i]]
        }
    }
    acframe <- frame_ac(bterms)
    sar_M_names <- get_ac_vars(acframe, "M", class = "sar")
    for (M in sar_M_names) {
        data2[[M]] <- validate_sar_matrix(get_from_data2(M, data2))
        attr(data2[[M]], "obs_based_matrix") <- TRUE
    }
    car_M_names <- get_ac_vars(acframe, "M", class = "car")
    for (M in car_M_names) {
        data2[[M]] <- validate_car_matrix(get_from_data2(M, data2))
    }
    fcor_M_names <- get_ac_vars(acframe, "M", class = "fcor")
    for (M in fcor_M_names) {
        data2[[M]] <- validate_fcor_matrix(get_from_data2(M, 
            data2))
        attr(data2[[M]], "obs_based_matrix") <- TRUE
    }
    cov_names <- ufrom_list(get_re(bterms)$gcall, "cov")
    cov_names <- cov_names[nzchar(cov_names)]
    for (cov in cov_names) {
        data2[[cov]] <- validate_recov_matrix(get_from_data2(cov, 
            data2))
    }
    data2
}
