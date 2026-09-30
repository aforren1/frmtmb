get_matches <- 
function (pattern, text, simplify = TRUE, first = FALSE, ...) 
{
    x <- regmatches(text, gregexpr(pattern, text, ...))
    if (first) {
        x <- lapply(x, function(t) if (length(t)) 
            t[1]
        else t)
    }
    if (simplify) {
        if (first) {
            x <- lapply(x, function(t) if (length(t)) 
                t
            else "")
        }
        x <- unlist(x)
    }
    x
}
