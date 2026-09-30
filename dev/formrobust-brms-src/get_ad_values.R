get_ad_values <- 
function (x, ad, name, data) 
{
    expr <- get_ad_expr(x, ad, name, type = "vars")
    eval2(expr, data)
}
