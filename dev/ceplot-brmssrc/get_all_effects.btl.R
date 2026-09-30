get_all_effects.btl <- 
function (x, ...) 
{
    c(get_var_combs(x[["fe"]], x[["cs"]]), get_all_effects_type(x, 
        "sp"), get_all_effects_type(x, "sm"), get_all_effects_type(x, 
        "gp"))
}
