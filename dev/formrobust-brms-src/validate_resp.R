validate_resp <- 
function (resp, x, multiple = TRUE) 
{
    if (is.brmsfit(x)) {
        x <- brmsterms(x$formula)$responses
    }
    x <- as.character(x)
    if (!length(x)) {
        return(NULL)
    }
    if (length(resp)) {
        resp <- as.character(resp)
        if (!all(resp %in% x)) {
            stop2("Invalid argument 'resp'. Valid response ", 
                "variables are: ", collapse_comma(x))
        }
        if (!multiple) {
            resp <- as_one_character(resp)
        }
    }
    else {
        resp <- x
    }
    resp
}
