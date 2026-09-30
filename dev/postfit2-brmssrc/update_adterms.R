update_adterms <- 
function (formula, adform, action = c("update", "replace")) 
{
    formula <- as_formula(formula)
    adform <- as_formula(adform)
    action <- match.arg(action)
    if (is.null(lhs(formula))) {
        stop2("Can't update a ond-sided formula.")
    }
    str_formula <- formula2str(formula)
    old_ad <- get_matches("(?<=\\|)[^~]*(?=~)", str_formula, 
        perl = TRUE)
    new_ad_terms <- attr(terms(adform), "term.labels")
    if (action == "update" && length(old_ad)) {
        old_ad <- formula(paste("~", old_ad))
        old_ad_terms <- attr(terms(old_ad), "term.labels")
        old_adnames <- get_matches("^[^\\(]+", old_ad_terms)
        old_adnames <- sub("^resp_", "", old_adnames)
        new_adnames <- get_matches("^[^\\(]+", new_ad_terms)
        new_adnames <- sub("^resp_", "", new_adnames)
        keep <- !old_adnames %in% new_adnames
        new_ad_terms <- c(old_ad_terms[keep], new_ad_terms)
    }
    if (length(new_ad_terms)) {
        new_ad_terms <- paste(new_ad_terms, collapse = "+")
        new_ad_terms <- paste("|", new_ad_terms)
    }
    resp <- gsub("\\|.+", "", deparse0(formula[[2]]))
    out <- formula(paste(resp, new_ad_terms, "~1"))
    out[[3]] <- formula[[3]]
    attributes(out) <- attributes(formula)
    out
}
