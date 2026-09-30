prepare_cond_data <- 
function (data, conditions, int_conditions = NULL, int_vars = NULL, 
    group_vars = NULL, surface = FALSE, resolution = 100, reorder = TRUE) 
{
    effects <- names(data)
    stopifnot(length(effects) %in% c(1, 2))
    is_factor <- ulapply(data, is_like_factor) | names(data) %in% 
        group_vars
    types <- ifelse(is_factor, "factor", "numeric")
    if (reorder) {
        new_order <- order(types, decreasing = TRUE)
        effects <- effects[new_order]
        types <- types[new_order]
    }
    if (effects[1] %in% names(int_conditions)) {
        int_cond <- int_conditions[[effects[1]]]
        if (is.function(int_cond)) {
            int_cond <- int_cond(data[[effects[1]]])
        }
        values <- int_cond
    }
    else if (types[1] == "factor") {
        values <- factor(unique(data[[effects[1]]]))
    }
    else {
        min1 <- min(data[[effects[1]]], na.rm = TRUE)
        max1 <- max(data[[effects[1]]], na.rm = TRUE)
        if (effects[1] %in% int_vars) {
            values <- seq(min1, max1, by = 1)
        }
        else {
            values <- seq(min1, max1, length.out = resolution)
        }
    }
    if (length(effects) == 2) {
        values <- setNames(list(values, NA), effects)
        if (effects[2] %in% names(int_conditions)) {
            int_cond <- int_conditions[[effects[2]]]
            if (is.function(int_cond)) {
                int_cond <- int_cond(data[[effects[2]]])
            }
            values[[2]] <- int_cond
        }
        else if (types[2] == "factor") {
            values[[2]] <- factor(unique(data[[effects[2]]]))
        }
        else {
            if (surface) {
                min2 <- min(data[[effects[2]]], na.rm = TRUE)
                max2 <- max(data[[effects[2]]], na.rm = TRUE)
                if (effects[2] %in% int_vars) {
                  values[[2]] <- seq(min2, max2, by = 1)
                }
                else {
                  values[[2]] <- seq(min2, max2, length.out = resolution)
                }
            }
            else {
                if (effects[2] %in% int_vars) {
                  median2 <- median(data[[effects[2]]])
                  mad2 <- mad(data[[effects[2]]])
                  values[[2]] <- round((-1:1) * mad2 + median2)
                }
                else {
                  mean2 <- mean(data[[effects[2]]], na.rm = TRUE)
                  sd2 <- sd(data[[effects[2]]], na.rm = TRUE)
                  values[[2]] <- (-1:1) * sd2 + mean2
                }
            }
        }
        data <- do_call(expand.grid, values)
    }
    else {
        stopifnot(length(effects) == 1)
        data <- structure(data.frame(values), names = effects)
    }
    data <- unique(data)
    data <- data[do_call(order, unname(as.list(data))), , drop = FALSE]
    data <- replicate(nrow(conditions), data, simplify = FALSE)
    cond_vars <- setdiff(names(conditions), effects)
    cond__ <- get_cond__(conditions)
    for (j in seq_rows(conditions)) {
        data[[j]] <- fill_newdata(data[[j]], cond_vars, conditions, 
            n = j)
        data[[j]]$cond__ <- cond__[j]
    }
    data <- do_call(rbind, data)
    data$cond__ <- factor(data$cond__, cond__)
    structure(data, effects = effects, types = types)
}
