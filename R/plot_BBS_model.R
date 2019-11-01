#' Plot trend results from BBS-style hierarchical model for point count survey data
#'
#' @param type Type of plot to produce: plots of annual indices of abundance
#' ('index') or plots of the linear trend lines ('trend'), both ('both'). Defaults to 'both'.
#' @param df Summary of median, lower, and upper estimates from \code{\link{get_BBS_model_estimates}}
#' @details Produces line and ribbon for trend estimates, and points with error bars for index estimates.
#' @return returns ggplot
#' @export
#' @import ggplot2
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#'
plot_BBS_model <- function(df, type = 'both') {
  p <- ggplot(df) +
    theme_classic() +
    xlab(NULL) +
    ylab('Abundance index')

  if (type %in% c('both', 'trend')) {
    p <- p +
      geom_ribbon(data = df %>% filter(type == 'trend'),
                  aes(x = .data$Year,
                      ymin = .data$lower,
                      ymax = .data$upper),
                  alpha = 0.2, fill = 'blue') +
      geom_line(data = df %>% filter(type == 'trend'),
                aes(x = .data$Year,
                    y = .data$median),
                size = 1, color = 'blue')
  }
  if (type %in% c('both', 'index')) {
    p <- p +
      geom_pointrange(data = df %>% filter(type == 'index'),
                      aes(x = .data$Year,
                          y = .data$median,
                          ymin = .data$lower,
                          ymax = .data$upper)) +
      geom_line(data = df %>% filter(type == 'index'),
                aes(x = .data$Year,
                    y = .data$median),
                size = 1)
  }

  # if (length(unique(df$species)) > 1) {
  #   if (length(unique(df$Project)) > 1) {
  #     p <- p + facet_wrap(~.data$species + .data$Project)
  #   } else {
  #     p <- p + facet_wrap(~.data$species)
  #   }
  # } else if (length(unique(df$Project)) > 1) {
  #   p <- p + facet_wrap(~.data$Project)
  # }
  print(p)
  return(p)
}
