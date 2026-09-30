get_all_effects_type <- 
function (x, type) 
{
    stopifnot(is.btl(x))
    type <- as_one_character(type)
    regex_type <- regex_sp(type)
    terms <- all_terms(x[[type]])
    out <- named_list(terms)
    for (i in seq_along(terms)) {
        term_parts <- unlist(strsplit(terms[i], split = ":"))
        vars <- vector("list", length(term_parts))
        for (j in seq_along(term_parts)) {
            matches <- get_matches_expr(regex_type, term_parts[j])
            for (k in seq_along(matches)) {
                tmp <- eval2(matches[[k]])
                c(vars[[j]]) <- setdiff(unique(c(tmp$term, tmp$by)), 
                  "NA")
            }
            c(vars[[j]]) <- setdiff(all_vars(term_parts[j]), 
                all_vars(matches))
        }
        vars <- unique(unlist(vars))
        out[[i]] <- str2formula(vars, collapse = "*")
    }
    get_var_combs(alist = out)
}
