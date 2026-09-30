plot.brms_conditional_effects <- 
function (x, ncol = NULL, points = getOption("brms.plot_points", 
    FALSE), rug = getOption("brms.plot_rug", FALSE), mean = TRUE, 
    jitter_width = 0, stype = c("contour", "raster"), line_args = list(), 
    cat_args = list(), errorbar_args = list(), surface_args = list(), 
    spaghetti_args = list(), point_args = list(), rug_args = list(), 
    facet_args = list(), theme = NULL, ask = TRUE, plot = TRUE, 
    ...) 
{
    dots <- list(...)
    plot <- use_alias(plot, dots$do_plot)
    stype <- match.arg(stype)
    smooths_only <- isTRUE(attr(x, "smooths_only"))
    if (points && smooths_only) {
        stop2("Argument 'points' is invalid for objects ", "returned by 'conditional_smooths'.")
    }
    if (!is_equal(jitter_width, 0)) {
        warning2("'jitter_width' is deprecated. Please use ", 
            "'point_args = list(width = <width>)' instead.")
    }
    if (!is.null(theme) && !is.theme(theme)) {
        stop2("Argument 'theme' should be a 'theme' object.")
    }
    if (plot) {
        default_ask <- devAskNewPage()
        on.exit(devAskNewPage(default_ask))
        devAskNewPage(ask = FALSE)
    }
    dont_replace <- c("mapping", "data", "inherit.aes")
    plots <- named_list(names(x))
    for (i in seq_along(x)) {
        response <- attr(x[[i]], "response")
        effects <- attr(x[[i]], "effects")
        ncond <- length(unique(x[[i]]$cond__))
        df_points <- attr(x[[i]], "points")
        categorical <- isTRUE(attr(x[[i]], "categorical"))
        catscale <- attr(x[[i]], "catscale")
        surface <- isTRUE(attr(x[[i]], "surface"))
        ordinal <- isTRUE(attr(x[[i]], "ordinal"))
        if (surface || ordinal) {
            plots[[i]] <- ggplot(x[[i]]) + aes(.data[["effect1__"]], 
                .data[["effect2__"]]) + labs(x = effects[1], 
                y = effects[2])
            if (ordinal) {
                width <- ifelse(is_like_factor(x[[i]]$effect1__), 
                  0.9, 1)
                .surface_args <- nlist(mapping = aes(fill = .data[["estimate__"]]), 
                  height = 0.9, width = width)
                replace_args(.surface_args, dont_replace) <- surface_args
                plots[[i]] <- plots[[i]] + do_call(geom_tile, 
                  .surface_args) + scale_fill_gradientn(colors = viridis6(), 
                  name = catscale) + ylab(response)
            }
            else if (stype == "contour") {
                .surface_args <- nlist(mapping = aes(z = .data[["estimate__"]], 
                  colour = after_stat(.data[["level"]])), bins = 30, 
                  linewidth = 1.3)
                replace_args(.surface_args, dont_replace) <- surface_args
                plots[[i]] <- plots[[i]] + do_call(geom_contour, 
                  .surface_args) + scale_color_gradientn(colors = viridis6(), 
                  name = response)
            }
            else if (stype == "raster") {
                .surface_args <- nlist(mapping = aes(fill = .data[["estimate__"]]))
                replace_args(.surface_args, dont_replace) <- surface_args
                plots[[i]] <- plots[[i]] + do_call(geom_raster, 
                  .surface_args) + scale_fill_gradientn(colors = viridis6(), 
                  name = response)
            }
        }
        else {
            gvar <- if (length(effects) == 2) 
                "effect2__"
            spaghetti <- attr(x[[i]], "spaghetti")
            aes_tmp <- aes(x = .data[["effect1__"]], y = .data[["estimate__"]])
            if (!is.null(gvar)) {
                aes_tmp$colour <- aes(colour = .data[[gvar]])$colour
            }
            plots[[i]] <- ggplot(x[[i]]) + aes_tmp + labs(x = effects[1], 
                y = response, colour = effects[2])
            if (is.null(spaghetti)) {
                aes_tmp <- aes(ymin = .data[["lower__"]], ymax = .data[["upper__"]])
                if (!is.null(gvar)) {
                  aes_tmp$fill <- aes(fill = .data[[gvar]])$fill
                }
                plots[[i]] <- plots[[i]] + aes_tmp + labs(fill = effects[2])
            }
            colors <- ggplot_build(plots[[i]])
            colors <- unique(colors$data[[1]][["colour"]])
            if (points && !categorical && !surface) {
                .point_args <- list(mapping = aes(x = .data[["effect1__"]], 
                  y = .data[["resp__"]]), data = df_points, inherit.aes = FALSE, 
                  size = 2/ncond^0.25, height = 0, width = jitter_width)
                if (is_like_factor(df_points[, gvar])) {
                  .point_args$mapping[c("colour", "fill")] <- aes(colour = .data[[gvar]], 
                    fill = .data[[gvar]])
                }
                replace_args(.point_args, dont_replace) <- point_args
                plots[[i]] <- plots[[i]] + do_call(geom_jitter, 
                  .point_args)
            }
            if (!is.null(spaghetti)) {
                .spaghetti_args <- list(aes(group = .data[["sample__"]]), 
                  data = spaghetti, stat = "identity", linewidth = 0.5)
                if (!is.null(gvar)) {
                  .spaghetti_args[[1]]$colour <- aes(colour = .data[[gvar]])$colour
                }
                if (length(effects) == 1) {
                  .spaghetti_args$colour <- alpha("blue", 0.1)
                }
                else {
                  plots[[i]] <- plots[[i]] + scale_color_manual(values = alpha(colors, 
                    0.1))
                }
                replace_args(.spaghetti_args, dont_replace) <- spaghetti_args
                plots[[i]] <- plots[[i]] + do_call(geom_smooth, 
                  .spaghetti_args)
            }
            if (is.numeric(x[[i]]$effect1__)) {
                .line_args <- list(stat = "identity")
                if (!is.null(spaghetti)) {
                  if (!is.null(gvar)) {
                    .line_args$mapping <- aes(group = .data[[gvar]])
                  }
                  .line_args$colour <- alpha("white", 0.8)
                }
                replace_args(.line_args, dont_replace) <- line_args
                if (mean || is.null(spaghetti)) {
                  plots[[i]] <- plots[[i]] + do_call(geom_smooth, 
                    .line_args)
                }
                if (rug) {
                  .rug_args <- list(mapping = aes(x = .data[["effect1__"]]), 
                    sides = "b", data = df_points, inherit.aes = FALSE)
                  if (is_like_factor(df_points[, gvar])) {
                    .rug_args$mapping["colour"] <- aes(colour = .data[[gvar]])
                  }
                  replace_args(.rug_args, dont_replace) <- rug_args
                  plots[[i]] <- plots[[i]] + do_call(geom_rug, 
                    .rug_args)
                }
            }
            else {
                .cat_args <- list(position = position_dodge(width = 0.4), 
                  size = 4/ncond^0.25)
                .errorbar_args <- list(position = position_dodge(width = 0.4), 
                  width = 0.3)
                replace_args(.cat_args, dont_replace) <- cat_args
                replace_args(.errorbar_args, dont_replace) <- errorbar_args
                plots[[i]] <- plots[[i]] + do_call(geom_point, 
                  .cat_args) + do_call(geom_errorbar, .errorbar_args)
            }
            if (categorical) {
                plots[[i]] <- plots[[i]] + ylab(catscale) + labs(fill = response, 
                  color = response)
            }
        }
        if (ncond > 1) {
            if (is.null(ncol)) {
                ncol <- max(floor(sqrt(ncond)), 3)
            }
            .facet_args <- nlist(facets = "cond__", ncol)
            replace_args(.facet_args, dont_replace) <- facet_args
            plots[[i]] <- plots[[i]] + do_call(facet_wrap, .facet_args)
        }
        plots[[i]] <- plots[[i]] + theme
        if (plot) {
            plot(plots[[i]])
            if (i == 1) {
                devAskNewPage(ask = ask)
            }
        }
    }
    invisible(plots)
}
