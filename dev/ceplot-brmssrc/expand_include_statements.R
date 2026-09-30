expand_include_statements <- 
function (model) 
{
    path <- system.file("chunks", package = "brms")
    includes <- unique(get_matches("#include '[^']+'", model))
    files <- gsub("(#include )|(')", "", includes)
    for (i in seq_along(includes)) {
        code <- readLines(paste0(path, "/", files[i]))
        code <- paste0(code, collapse = "\n")
        pattern <- paste0(" *", escape_all(includes[i]))
        model <- sub(pattern, code, model)
        model <- gsub(pattern, "", model)
    }
    includes <- unique(get_matches("#includeR `[^`]+`", model))
    calls <- gsub("(#includeR )|(`)", "", includes)
    for (i in seq_along(includes)) {
        code <- eval2(calls[i])
        pattern <- paste0(" *", escape_all(includes[i]))
        model <- sub(pattern, code, model)
        model <- gsub(pattern, "", model)
    }
    model
}
