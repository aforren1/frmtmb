bf <- 
function (formula, ..., flist = NULL, family = NULL, autocor = NULL, 
    nl = NULL, loop = NULL, center = NULL, cmc = NULL, sparse = NULL, 
    decomp = NULL) 
{
    brmsformula(formula, ..., flist = flist, family = family, 
        autocor = autocor, nl = nl, loop = loop, center = center, 
        cmc = cmc, sparse = sparse, decomp = decomp)
}
