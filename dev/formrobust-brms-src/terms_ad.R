terms_ad <- 
function (formula, family = NULL, check_response = TRUE) 
{
    x <- list()
    ad_funs <- lsp("brms", what = "exports", pattern = "^resp_")
    ad_funs <- sub("^resp_", "", ad_funs)
    families <- family_names(family)
    if (is.family(family) && any(nzchar(families))) {
        str_formula <- formula2str(formula)
        ad <- get_matches("(?<=\\|)[^~]*(?=~)", str_formula, 
            perl = TRUE)
        valid_ads <- family_info(family, "ad")
        if (length(ad)) {
            ad_terms <- terms(str2formula(ad))
            if (length(attr(ad_terms, "offset"))) {
                stop2("Offsets are not allowed in addition terms.")
            }
            ad_terms <- attr(ad_terms, "term.labels")
            for (a in ad_funs) {
                matches <- grep(paste0("^(resp_)?", a, "\\(.*\\)$"), 
                  ad_terms)
                if (length(matches) == 1) {
                  x[[a]] <- ad_terms[matches]
                  if (!grepl("^resp_", x[[a]])) {
                    x[[a]] <- paste0("resp_", x[[a]])
                  }
                  ad_terms <- ad_terms[-matches]
                  if (!is.na(x[[a]]) && a %in% valid_ads) {
                    x[[a]] <- str2formula(x[[a]])
                  }
                  else {
                    stop2("Argument '", a, "' is not supported for ", 
                      "family '", summary(family), "'.")
                  }
                }
                else if (length(matches) > 1) {
                  stop2("Each addition argument may only be defined once.")
                }
            }
            if (length(ad_terms)) {
                stop2("The following addition terms are invalid:\n", 
                  collapse_comma(ad_terms))
            }
        }
        if (check_response && "wiener" %in% families && !is.formula(x$dec)) {
            stop2("Addition argument 'dec' is required for family 'wiener'.")
        }
        if (is.formula(x$cat)) {
            x$thres <- x$cat
        }
    }
    x
}
