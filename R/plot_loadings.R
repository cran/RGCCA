#' Bar plots of the weights or the loadings
#'
#' Bar plots of the weights or the loadings sorted in decreasing order.
#' @inheritParams plot.rgcca
#' @param df Data frame with the data to plot.
#' @param theme_RGCCA Theme of the plot.
#' @noRd
plot_loadings <- function(df, title, x, block, theme_RGCCA,
                          cex_point, var_colors, ...) {
  # Add colors depending on looking at superblock or regular blocks
  is_multiblock <- (length(block) > 1) || (block == length(x$call$blocks) + 1)
  if (is_multiblock) {
    p <- ggplot(df, aes(x = .data$x, y = .data$y, color = .data$response)) +
      ggplot2::scale_color_manual(values = var_colors) +
      ggplot2::labs(color = "Block")
  } else {
    p <- ggplot(df, aes(x = .data$x, y = .data$y))
  }
  # Construct plot
  p <- p +
    ggplot2::geom_point(size = .5 * cex_point) +
    ggplot2::geom_errorbar(aes(
      xmin = 0,
      xmax = .data$x,
      width = 0
    ), linewidth = .2 * cex_point) +
    ggplot2::geom_vline(
      xintercept = 0, lty = "longdash", linewidth = .12 * cex_point
    ) +
    theme_RGCCA +
    ggplot2::labs(title = title, x = "", y = "") +
    ggplot2::scale_y_discrete(limits = rev) +
    ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5, 0, 0, 0, "mm"),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank()
    )

  # Hide legend if not superblock
  if (!is_multiblock) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  return(p)
}
