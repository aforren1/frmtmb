function (x, ...) 
{
    out <- list()
    for (dp in names(x$dpars)) {
        c(out) <- conditional_smooths(x$dpars[[dp]], ...)
    }
    for (nlp in names(x$nlpars)) {
        c(out) <- conditional_smooths(x$nlpars[[nlp]], ...)
    }
    out
}
