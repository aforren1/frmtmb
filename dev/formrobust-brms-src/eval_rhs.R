eval_rhs <- 
function (formula, data = NULL) 
{
    formula <- as.formula(formula)
    eval(rhs(formula)[[2]], data, environment(formula))
}
