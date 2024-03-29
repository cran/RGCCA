#' Bar plots of Average Variance Explained
#'
#' Bar plots of the model quality (based on Average Variance Explained)
#' for each blocks and sorted in decreasing order.
#' @inheritParams plot.rgcca
#' @param df Data frame with the data to plot.
#' @param theme_RGCCA Theme of the plot.
#' @noRd
plot_ave <- function(df, title, x, theme_RGCCA,
                     cex_sub, cex_point, AVE_colors, ...) {
  # Construct plot
  p <- ggplot(df, aes(x = .data$AVE, y = .data$block, fill = .data$comp)) +
    ggplot2::geom_bar(
      position = ggplot2::position_stack(reverse = TRUE), stat = "identity"
    ) +
    ggplot2::stat_identity(
      geom = "text", color = "black", size = cex_point,
      aes(label = .data$label),
      position = ggplot2::position_stack(reverse = TRUE, vjust = 0.5)
    ) +
    ggplot2::scale_fill_manual(values = AVE_colors) +
    theme_RGCCA +
    ggplot2::labs(subtitle = print_comp(x, outer = TRUE)) +
    ggplot2::labs(title = title, x = "", y = "", fill = "Component") +
    ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_text(
        size = cex_sub,
        face = "italic"
      ),
      plot.margin = ggplot2::margin(5, 0, 0, 0, "mm")
    )
  return(p)
}
