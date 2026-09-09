#' Plot trend results from BBS-style hierarchical model for point count survey data
#' 
#' Produces line and ribbon for trend estimates, and points with error bars for index estimates using \pkg{ggplot2}. Install it with: \code{install.packages("ggplot2")}
#'
#' @param type Type of plot to produce: plots of annual indices of abundance
#' ('index') or plots of the linear trend lines ('trend'), both ('both'). Defaults to 'both'.
#' @param df Summary of median, lower, and upper estimates from \code{\link{get_BBS_model_estimates}}
#' 
#' @return returns ggplot
#' 
#' @importFrom rlang .data
#' @export

plot_BBS_model <- function(df, type = 'both') {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for plot_BBS_model(). ",
      "Install it with install.packages('ggplot2').",
      call. = FALSE
    )
  }
  
  p <- ggplot2::ggplot(df) + ggplot2::theme_classic() +
    ggplot2::xlab(NULL) +
    ggplot2::ylab('Abundance index')

  if (type %in% c('both', 'trend')) {
    p <- p +
      ggplot2::geom_ribbon(
        data = df |> dplyr::filter(type == 'trend'),
        ggplot2::aes(x = .data$Year, ymin = .data$lower, ymax = .data$upper),
        alpha = 0.2, fill = 'blue') +
      ggplot2::geom_line(
        data = df |> dplyr::filter(type == 'trend'),
        ggplot2::aes(x = .data$Year, y = .data$median),
        size = 1, color = 'blue')
  }
  
  if (type %in% c('both', 'index')) {
    p <- p +
      ggplot2::geom_pointrange(
        data = df |> dplyr::filter(type == 'index'),
        ggplot2::aes(x = .data$Year,
                     y = .data$median,
                     ymin = .data$lower,
                     ymax = .data$upper)) +
      ggplot2::geom_line(
        data = df |> dplyr::filter(type == 'index'),
        ggplot2::aes(x = .data$Year, y = .data$median),
        size = 1)
  }

  print(p)
  return(p)
}
