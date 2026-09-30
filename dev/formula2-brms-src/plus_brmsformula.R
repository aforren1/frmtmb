function (e1, e2) 
{
    if (is.function(e2)) {
        e2 <- try(e2(), silent = TRUE)
        if (!is.family(e2)) {
            stop2("Don't know how to handle non-family functions.")
        }
    }
    if (is.family(e2)) {
        e1 <- bf(e1, family = e2)
    }
    else if (is.cor_brms(e2) || inherits(e2, "acformula")) {
        e1 <- bf(e1, autocor = e2)
    }
    else if (inherits(e2, "setnl")) {
        dpar <- attr(e2, "dpar")
        if (is.null(dpar)) {
            e1 <- bf(e1, nl = e2)
        }
        else {
            if (is.null(e1$pforms[[dpar]])) {
                stop2("Parameter '", dpar, "' has no formula.")
            }
            attr(e1$pforms[[dpar]], "nl") <- e2
            e1 <- bf(e1)
        }
    }
    else if (inherits(e2, "setmecor")) {
        e1$mecor <- e2[1]
    }
    else if (is.brmsformula(e2)) {
        e1 <- mvbf(e1, e2)
    }
    else if (inherits(e2, "setrescor")) {
        stop2("Setting 'rescor' is only possible in multivariate models.")
    }
    else if (is.ac_term(e2)) {
        stop2("Autocorrelation terms can only be specified on the right-hand ", 
            "side of a formula, not added to a 'brmsformula' object.")
    }
    else if (!is.null(e2)) {
        e1 <- bf(e1, e2)
    }
    e1
}
