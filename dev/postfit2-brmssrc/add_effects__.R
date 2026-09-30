add_effects__ <- 
function (data, effects) 
{
    for (i in seq_along(effects)) {
        data[[paste0("effect", i, "__")]] <- eval2(effects[i], 
            data)
    }
    data
}
