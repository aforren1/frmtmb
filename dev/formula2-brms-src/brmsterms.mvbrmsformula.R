function (formula, ...) 
{
    x <- validate_formula(formula)
    x$rescor <- isTRUE(x$rescor)
    x$mecor <- isTRUE(x$mecor)
    out <- structure(list(), class = "mvbrmsterms")
    out$terms <- named_list(names(x$forms))
    for (i in seq_along(out$terms)) {
        x$forms[[i]]$rescor <- x$rescor
        x$forms[[i]]$mecor <- x$mecor
        x$forms[[i]]$mv <- TRUE
        out$terms[[i]] <- brmsterms(x$forms[[i]], ...)
    }
    list_allvars <- lapply(out$terms, get_allvars)
    out$allvars <- allvars_formula(list_allvars, .env = environment(list_allvars[[1]]))
    lhs_resp <- function(x) deparse0(lhs(x$respform)[[2]])
    out$respform <- paste0(ulapply(out$terms, lhs_resp), collapse = ",")
    out$respform <- formula(paste0("mvbind(", out$respform, ") ~ 1"))
    out$responses <- ufrom_list(out$terms, "resp")
    out$rescor <- x$rescor
    out$mecor <- x$mecor
    out$cov_ranef <- x$cov_ranef
    out
}
