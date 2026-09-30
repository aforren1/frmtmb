hypothesis.default <- 
function (x, hypothesis, alpha = 0.05, robust = FALSE, ...) 
{
    x <- as.data.frame(x)
    .hypothesis(x, hypothesis, class = "", alpha = alpha, robust = robust, 
        ...)
}
