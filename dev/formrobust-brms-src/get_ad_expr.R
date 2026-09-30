get_ad_expr <- 
function (x, ad, name, type = "vars") 
{
    ad <- as_one_character(ad)
    name <- as_one_character(name)
    type <- as_one_character(type)
    if (is.null(x$adforms[[ad]])) {
        return(NULL)
    }
    out <- eval_rhs(x$adforms[[ad]])[[type]][[name]]
    if (type == "vars" && is_equal(out, "NA")) {
        out <- NULL
    }
    out
}
