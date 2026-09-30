function (object, formula., mode = c("update", "replace", "keep"), 
    ...) 
{
    mode <- match.arg(mode)
    object <- bf(object)
    up_nl <- get_nl(formula., aol = FALSE)
    if (is.null(up_nl)) {
        up_nl <- get_nl(object)
    }
    formula. <- bf(formula., nl = up_nl)
    up_family <- formula.[["family"]]
    if (is.null(up_family)) {
        up_family <- object[["family"]]
    }
    up_autocor <- attr(formula.$formula, "autocor")
    if (is.null(up_autocor)) {
        up_autocor <- attr(object$formula, "autocor")
    }
    old_form <- object$formula
    up_form <- formula.$formula
    if (mode == "update") {
        new_form <- update(old_form, up_form, ...)
    }
    else if (mode == "replace") {
        new_form <- up_form
    }
    else if (mode == "keep") {
        new_form <- old_form
    }
    flist <- c(object$pforms, object$pfix, formula.$pforms, formula.$pfix)
    bf(new_form, flist = flist, family = up_family, autocor = up_autocor, 
        nl = up_nl)
}
