.prepare_predictions_gp <- 
function (gpframe, draws, sdata, nug, new, p, i, byj = NULL) 
{
    sfx1 <- escape_all(gpframe$sfx1[[i]])
    sfx2 <- escape_all(gpframe$sfx2[[i]])
    if (is.null(byj)) {
        lvl <- ""
    }
    else {
        lvl <- gpframe$bylevels[[i]][byj]
        sfx1 <- sfx1[byj]
        sfx2 <- sfx2[byj, ]
    }
    j <- usc(byj)
    pi <- paste0(p, "_", i)
    gp <- list()
    gp$cov <- gpframe$cov[i]
    sdgp <- paste0("^sdgp", p, "_", sfx1, "$")
    gp$sdgp <- as.vector(prepare_draws(draws, sdgp, regex = TRUE))
    lscale <- paste0("^lscale", p, "_", sfx2, "$")
    gp$lscale <- prepare_draws(draws, lscale, regex = TRUE)
    zgp_regex <- paste0("^zgp", p, "_", sfx1, "\\[")
    gp$zgp <- prepare_draws(draws, zgp_regex, regex = TRUE)
    Xgp_name <- paste0("Xgp", pi, j)
    Igp_name <- paste0("Igp", pi, j)
    Jgp_name <- paste0("Jgp", pi, j)
    if (new && isNA(gpframe$k[i])) {
        gp$x <- sdata[[paste0(Xgp_name, "_old")]]
        gp$nug <- 1e-12
        gp$yL <- .predictor_gp(gp)
        gp$x_new <- sdata[[Xgp_name]]
        gp$Igp <- sdata[[Igp_name]]
    }
    else {
        gp$x <- sdata[[Xgp_name]]
        gp$Igp <- sdata[[Igp_name]]
        if (!isNA(gpframe$k[i])) {
            gp$slambda <- sdata[[paste0("slambda", pi, j)]]
        }
    }
    gp$Jgp <- sdata[[Jgp_name]]
    gp$Cgp <- sdata[[paste0("Cgp", pi, j)]]
    gp$nug <- nug
    gp
}
