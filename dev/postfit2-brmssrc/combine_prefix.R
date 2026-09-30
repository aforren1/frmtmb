combine_prefix <- 
function (prefix, keep_mu = FALSE, nlp = FALSE) 
{
    prefix <- check_prefix(prefix, keep_mu = keep_mu)
    if (is_nlpar(prefix) && nlp) {
        prefix$dpar <- "nlp"
    }
    prefix <- lapply(prefix, usc)
    sub("^_", "", do_call(paste0, prefix))
}
