terms_ac <- 
function (formula) 
{
    autocor <- attr(formula, "autocor")
    out <- c(find_terms(formula, "ac"), find_terms(autocor, "ac"))
    if (!length(out)) {
        return(NULL)
    }
    eterms <- lapply(out, eval2, envir = environment())
    allvars <- unlist(c(from_list(eterms, "time"), from_list(eterms, 
        "gr")))
    allvars <- str2formula(all_vars(allvars))
    out <- str2formula(out)
    attr(out, "allvars") <- allvars
    out
}
