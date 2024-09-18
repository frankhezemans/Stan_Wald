# define function to plot the results
plot_wald <- function(
  data, fun, facet_scales = "fixed", xlim = c(-0.5, 6.5)
) {

  require(dplyr)
  require(tidyr)
  require(purrr)
  require(forcats)
  require(scales)
  require(ggplot2)
  require(RColorBrewer)
  require(cowplot)

  y_axis_title <- c(
    "lpdf" = "value of LPDF",
    "lcdf" = "value of LCDF",
    "lccdf" = "value of LCCDF"
  )

  # in anticipation of later plotting, get a separate data frame with just the
  # unique set of parameter values (v, s, B)
  unique_param_sets <- data |>
    dplyr::select(
      dplyr::all_of(c("v", "s", "B", "method"))
    ) |>
    dplyr::distinct() |>
    dplyr::mutate(
      method = factor(
        x = method,
        levels = c("statmod", "Stan_Wald", "brms", "DMC"),
        ordered = TRUE
      ),
      mu = (B / s) / (v / s),
      lambda = (B / s)^2
    ) |>
    dplyr::arrange(
      dplyr::desc(mu), dplyr::desc(lambda)
    ) |>
    tidyr::unite(
      col = "params",
      tidyr::all_of(c("v", "s", "B")),
      remove = FALSE
    ) |>
    dplyr::mutate(
      labels = paste0("v = ", v, ", s = ", s, ", B = ", B)
    )

  data <- data |>
    # create factor column of unique parameter sets, and turn 'method' into
    # ordered factor (for plotting purposes)
    tidyr::unite(
      col = "params",
      tidyr::all_of(c("v", "s", "B"))
    ) |>
    dplyr::mutate(
      params = factor(
        x = params,
        levels = unique_param_sets[["params"]],
        labels = unique_param_sets[["labels"]],
        ordered = TRUE
      ),
      method = factor(
        x = method,
        levels = c("statmod", "Stan_Wald", "brms", "DMC"),
        ordered = TRUE
      )
    ) |>
    dplyr::arrange(method, params) |>
    dplyr::mutate(
      method_params = forcats::fct_cross(
        method, params,
        sep = ": ",
        keep_empty = TRUE
      ),
      method_params = factor(
        x = method_params,
        levels = unique(method_params),
        ordered = TRUE
      )
    )


  create_plot <- function(data, y_title) {

    get_colours <- function(
      base_colour,
      n_cols = length(levels(data[["params"]])),
      range = c(0, 0.8)
    ) {
      ramp <- scales::colour_ramp(colors = c(base_colour, "#FFFFFF"))
      ramp(seq(from = range[1], to = range[2], length = n_cols))
    }

    base_colours <- RColorBrewer::brewer.pal(
      name = "Set1",
      n = length(levels(data[["method"]]))
    )

    colour_values <- unlist(purrr::map(
      .x = base_colours,
      .f = get_colours
    ))
    names(colour_values) <- rev(levels(data[["method_params"]]))

    fill_values <- get_colours("#000000")
    names(fill_values) <- rev(levels(data[["params"]]))

    out <- ggplot2::ggplot(
      data = data,
      mapping = ggplot2::aes(
        x = x,
        y = value
      )
    ) +
      ggplot2::facet_wrap(
        facets = ggplot2::vars(method),
        nrow = 1,
        scales = facet_scales,
        drop = FALSE
      ) +
      ggplot2::geom_point(
        mapping = ggplot2::aes(
          fill = params
        ),
        shape = 21,
        alpha = 0
      ) +
      ggplot2::geom_point(
        mapping = ggplot2::aes(
          colour = method_params
        ),
        size = 3,
        show.legend = FALSE
      ) +
      ggplot2::scale_fill_manual(
        values = fill_values
      ) +
      ggplot2::scale_colour_manual(
        values = colour_values
      ) +
      ggplot2::guides(
        fill = ggplot2::guide_legend(
          override.aes = list(
            alpha = 1,
            size = 5,
            colour = "#FFFFFF"
          )
        )
      ) +
      ggplot2::coord_cartesian(
        xlim = xlim
      ) +
      ggplot2::labs(
        y = y_title,
        x = "input value"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(
          colour = "#CCCCCC",
          linewidth = 1
        ),
        axis.ticks = ggplot2::element_line(
          colour = "#CCCCCC",
          linewidth = 1
        ),
        axis.title = ggplot2::element_text(
          size = 18
        ),
        axis.text = ggplot2::element_text(
          size = 12
        ),
        strip.text = ggplot2::element_text(
          size = 25,
          family = "mono"
        ),
        legend.position = "top",
        legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(
          size = 14,
          colour = "black"
        )
      )

    return(out)

  }

  out <- vector(mode = "list", length = 2)
  names(out) <- c("finite", "non_finite")

  out[["finite"]] <- create_plot(
    data = data |>
      dplyr::filter(is.finite(value)),
    y_title = y_axis_title[fun]
  )

  neg_inf_vals <- data[["value"]] == -Inf
  nan_vals <- is.nan(data[["value"]])

  if (any(neg_inf_vals[!is.na(neg_inf_vals)]) || any(nan_vals)) {

    out[["non_finite"]] <- create_plot(
      data = data |>
        dplyr::filter(value == -Inf | is.nan(value)) |>
        dplyr::mutate(
          value = dplyr::case_when(
            is.nan(value) ~ "nan",
            value == -Inf ~ "neg_inf",
            .default = NA_character_
          ),
          value = factor(
            x = value,
            levels = c("nan", "neg_inf"),
            labels = c("NaN", "-Inf"),
            ordered = TRUE
          )
        ),
      y_title = y_axis_title[fun]
    ) +
      ggplot2::theme(
        strip.text = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        legend.position = "null"
      )

    out[["finite"]] <- out[["finite"]] +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.line.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )

    rel_heights <- c(1, 0.35)

    combined_plot <- cowplot::plot_grid(
      plotlist = out[!sapply(out, is.null)],
      ncol = 1,
      rel_heights = rel_heights[!sapply(out, is.null)],
      axis = "lr", align = "v"
    )

    return(combined_plot)

  } else {

    return(out[["finite"]])

  }

}
