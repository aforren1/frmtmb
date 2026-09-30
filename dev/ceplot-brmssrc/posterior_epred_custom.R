posterior_epred_custom <- 
function (prep) 
{
    custom_family_method(prep$family, "posterior_epred")(prep)
}
