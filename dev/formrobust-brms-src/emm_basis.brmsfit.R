emm_basis.brmsfit <- 
function (object, trms, xlev, grid, vcov., resp = NULL, dpar = NULL, 
    nlpar = NULL, re_formula = NA, epred = FALSE, ...) 
{
    if (is_equal(dpar, "mean")) {
        warning2("dpar = 'mean' is deprecated. Please use epred = TRUE instead.")
        epred <- TRUE
        dpar <- NULL
    }
    epred <- as_one_logical(epred)
    bterms <- .extract_par_terms(object, resp = resp, dpar = dpar, 
        nlpar = nlpar, re_formula = re_formula, epred = epred)
    if (epred) {
        post.beta <- posterior_epred(object, newdata = grid, 
            re_formula = re_formula, resp = resp, incl_autocor = FALSE, 
            ...)
    }
    else {
        req_vars <- all_vars(bterms$allvars)
        post.beta <- posterior_linpred(object, newdata = grid, 
            re_formula = re_formula, resp = resp, dpar = dpar, 
            nlpar = nlpar, incl_autocor = FALSE, req_vars = req_vars, 
            transform = FALSE, offset = FALSE, ...)
    }
    if (anyNA(post.beta)) {
        stop2("emm_basis.brmsfit created NAs. Please check your reference grid.")
    }
    misc <- bterms$.misc
    if (length(dim(post.beta)) == 3) {
        ynames <- dimnames(post.beta)[[3]]
        if (is.null(ynames)) {
            ynames <- as.character(seq_len(dim(post.beta)[3]))
        }
        dims <- dim(post.beta)
        post.beta <- matrix(post.beta, ncol = prod(dims[2:3]))
        misc$ylevs = list(rep.meas = ynames)
    }
    attr(post.beta, "n.chains") <- object$fit@sim$chains
    X <- diag(ncol(post.beta))
    bhat <- apply(post.beta, 2, mean)
    V <- cov(post.beta)
    nbasis <- matrix(NA)
    dfargs <- list()
    dffun <- function(k, dfargs) Inf
    environment(dffun) <- baseenv()
    nlist(X, bhat, nbasis, V, dffun, dfargs, misc, post.beta)
}
