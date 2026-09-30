plot.brmshypothesis <- 
function (x, nvariables = 5, N = NULL, ignore_prior = FALSE, 
    chars = 40, colors = NULL, theme = NULL, ask = TRUE, plot = TRUE, 
    ...) 
{
    dots <- list(...)
    nvariables <- use_alias(nvariables, N)
    if (!is.data.frame(x$samples)) {
        stop2("No posterior draws found.")
    }
    plot <- use_alias(plot, dots$do_plot)
    if (is.null(colors)) {
        colors <- bayesplot::color_scheme_get()[c(4, 2)]
        colors <- unname(unlist(colors))
    }
    if (length(colors) != 2) {
        stop2("Argument 'colors' must be of length 2.")
    }
    .plot_fun <- function(samples) {
        samples <- na.omit(samples)
        ignore_prior <- ignore_prior || length(unique(samples$Type)) == 
            1
        gg <- ggplot(samples, aes(x = .data[["values"]])) + facet_wrap("ind", 
            ncol = 1, scales = "free") + xlab("") + ylab("") + 
            theme + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
        if (ignore_prior) {
            gg <- gg + geom_density(alpha = 0.7, fill = colors[1], 
                na.rm = TRUE)
        }
        else {
            gg <- gg + geom_density(aes(fill = .data[["Type"]]), 
                alpha = 0.7, na.rm = TRUE) + scale_fill_manual(values = colors)
        }
        return(gg)
    }
    samples <- cbind(x$samples, Type = "Posterior")
    if (!ignore_prior) {
        prior_samples <- cbind(x$prior_samples, Type = "Prior")
        samples <- rbind(samples, prior_samples)
    }
    if (plot) {
        default_ask <- devAskNewPage()
        on.exit(devAskNewPage(default_ask))
        devAskNewPage(ask = FALSE)
    }
    hyps <- limit_chars(x$hypothesis$Hypothesis, chars = chars)
    if (!is.null(x$hypothesis$Group)) {
        hyps <- paste0(x$hypothesis$Group, ":  ", hyps)
    }
    names(samples)[seq_along(hyps)] <- hyps
    nplots <- ceiling(length(hyps)/nvariables)
    plots <- vector(mode = "list", length = nplots)
    for (i in seq_len(nplots)) {
        sub <- ((i - 1) * nvariables + 1):min(i * nvariables, 
            length(hyps))
        sub_hyps <- hyps[sub]
        sub_samples <- cbind(utils::stack(samples[, sub_hyps, 
            drop = FALSE]), samples[, "Type", drop = FALSE])
        sub_samples$ind <- with(sub_samples, factor(ind, levels = unique(ind)))
        plots[[i]] <- .plot_fun(sub_samples)
        if (plot) {
            plot(plots[[i]])
            if (i == 1) 
                devAskNewPage(ask = ask)
        }
    }
    invisible(plots)
}
