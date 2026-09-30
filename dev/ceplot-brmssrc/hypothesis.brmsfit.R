hypothesis.brmsfit <- 
function (x, hypothesis, class = "b", group = "", scope = c("standard", 
    "ranef", "coef"), alpha = 0.05, robust = FALSE, seed = NULL, 
    ...) 
{
    if (!is.null(seed)) {
        set.seed(seed)
    }
    contains_draws(x)
    x <- restructure(x)
    group <- as_one_character(group)
    scope <- match.arg(scope)
    if (scope == "standard") {
        if (!length(class)) {
            class <- ""
        }
        class <- as_one_character(class)
        if (nzchar(group)) {
            class <- paste0(class, "_", group, "__")
        }
        else if (nzchar(class)) {
            class <- paste0(class, "_")
        }
        out <- .hypothesis(x, hypothesis, class = class, alpha = alpha, 
            robust = robust, ...)
    }
    else {
        co <- do_call(scope, list(x, summary = FALSE))
        if (!group %in% names(co)) {
            stop2("'group' should be one of ", collapse_comma(names(co)))
        }
        out <- hypothesis_coef(co[[group]], hypothesis, alpha = alpha, 
            robust = robust, ...)
    }
    out
}
