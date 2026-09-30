function (x, formula = FALSE, brackets = TRUE) 
{
    if (is.formula(x)) {
        x <- all_terms(x)
    }
    re_pos <- grepl("\\|", x)
    out <- x[re_pos]
    if (brackets && length(out)) {
        out <- paste0("(", out, ")")
    }
    if (formula) {
        out <- str2formula(out)
    }
    out
}
