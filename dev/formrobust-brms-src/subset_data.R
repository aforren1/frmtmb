subset_data <- 
function (data, bterms) 
{
    if (has_subset(bterms)) {
        subset <- as.logical(get_ad_values(bterms, "subset", 
            "subset", data))
        if (length(subset) != nrow(data)) {
            stop2("Length of 'subset' does not match the rows of 'data'.")
        }
        if (anyNA(subset)) {
            stop2("Subset variables may not contain NAs.")
        }
        check_cross_formula_indexing(bterms)
        data <- data[subset, , drop = FALSE]
        attr(data, "subset") <- subset
    }
    if (!NROW(data)) {
        stop2("All rows of 'data' were removed via 'subset'. ", 
            "Please make sure that variables do not contain NAs ", 
            "for observations in which they are supposed to be used. ", 
            "Please also make sure that each subset variable is ", 
            "TRUE for at least one observation.")
    }
    data
}
