formula2str <- 
function (formula, rm = c(0, 0), space = c("trim", "rm")) 
{
    if (is.null(formula)) {
        return(NULL)
    }
    formula <- as.formula(formula)
    space <- match.arg(space)
    if (anyNA(rm[2])) 
        rm[2] <- 0
    x <- Reduce(paste, deparse(formula))
    x <- gsub("[\t\r\n]+", " ", x, perl = TRUE)
    if (space == "trim") {
        x <- trim_wsp(x)
    }
    else {
        x <- rm_wsp(x)
    }
    substr(x, 1 + rm[1], nchar(x) - rm[2])
}
