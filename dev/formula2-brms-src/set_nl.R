function (nl = TRUE, dpar = NULL, resp = NULL) 
{
    nl <- as_one_logical(nl)
    if (!is.null(dpar)) {
        dpar <- as_one_character(dpar)
    }
    if (!is.null(resp)) {
        resp <- as_one_character(resp)
    }
    structure(nl, dpar = dpar, resp = resp, class = "setnl")
}
