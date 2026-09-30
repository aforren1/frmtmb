autocor.brmsfit <- 
function (object, resp = NULL, ...) 
{
    warning2("Method 'autocor' is deprecated and will be removed in the future.")
    object <- restructure(object)
    resp <- validate_resp(resp, object)
    if (!is.null(resp)) {
        autocor <- object$autocor[resp]
        if (length(resp) == 1) {
            autocor <- autocor[[1]]
        }
    }
    else {
        autocor <- object$autocor
    }
    autocor
}
