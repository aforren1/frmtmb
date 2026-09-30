prepare_predictions_offset <- 
function (bframe, sdata, ...) 
{
    p <- usc(combine_prefix(bframe))
    sdata[[paste0("offsets", p)]]
}
