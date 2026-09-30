function (x, empty_ok = TRUE) 
{
    out <- lhs(as_formula(x))
    if (is.null(out)) {
        if (empty_ok) {
            out <- ~1
        }
        else {
            str_x <- formula2str(x, space = "trim")
            stop2("Response variable is missing in formula ", 
                str_x)
        }
    }
    out <- gsub("\\|+[^~]*~", "~", formula2str(out))
    out <- try(formula(out), silent = TRUE)
    if (is_try_error(out)) {
        str_x <- formula2str(x, space = "trim")
        stop2("Incorrect use of '|' on the left-hand side of ", 
            str_x)
    }
    environment(out) <- environment(x)
    out
}
