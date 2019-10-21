#' Plot trend results from BBS-style hierarchical model for point count survey data
#'
#' @param modresults Model output from running \code{\link{fit_BBS_model}}
#' @param inputdat Input data from running \code{\link{setup_BBS_model}}
#' @param type Type of plot to produce: plots of annual indices of abundance
#' overall or by transect ('global_index', 'transect_index'), plots of
#' linear trend lines overall or by transect ('global_trend',
#' 'transect_trend'), or combinations of both ('global_trend+index',
#' 'transect_trend+index'). Defaults to 'global_trend+index', similar to
#' BBS output.
#' @details Runs \code{\link{get_BBS_model_estimates}} and produces ggplot
#' @return returns ggplot
#' @export
#' @import ggplot2
#' @importFrom dplyr %>%
#' @importFrom MCMCvis MCMCpstr
#' @importFrom HDInterval hdi
#' @importFrom stats median
#' @importFrom rlang .data
#'
plot_BBS_model <- function(modresults, inputdat,
                           type = 'global_trend+index') {
  if (type == 'global_trend+index') {
    res <- get_BBS_model_estimates(modresults, inputdat, type = 'global_index')
    res2 <- get_BBS_model_estimates(modresults, inputdat, type = 'global_trend')
  } else if (type == 'transect_trend+index') {
    res <- get_BBS_model_estimates(modresults, inputdat, type = 'transect_index')
    res2 <- get_BBS_model_estimates(modresults, inputdat, type = 'transect_trend')
  } else {
    res <- get_BBS_model_estimates(modresults, inputdat, type = type)
  }

  if (type %in% c('global_trend+index', 'transect_trend+index')) {
    p <- ggplot(res) +
      geom_pointrange(aes(x = .data$year.pred,
                          y = .data$median,
                          ymin = .data$lower,
                          ymax = .data$upper)) +
      geom_line(aes(x = .data$year.pred,
                    y = .data$median),
                size = 1) +
      geom_ribbon(data = res2,
                  aes(x = .data$year.pred,
                      ymin = .data$lower,
                      ymax = .data$upper),
                  alpha = 0.2, fill = 'blue') +
      geom_line(data = res2,
                aes(x = .data$year.pred,
                    y = .data$median,
                    color = .data$Transect),
                size = 1, color = 'blue') +
      theme_classic() +
      xlab(NULL) +
      ylab('Abundance index')

  } else {
    obsdat <- data.frame(Year = inputdat$dat$Year,
                         Transect = inputdat$dat$Transect,
                         Count = inputdat$observed)
    p <- ggplot(obsdat) +
      geom_jitter(aes(.data$Year, .data$Count),
                  width = 0.1, height = 0.1,
                  color = 'gray30', shape = 21, alpha = 0.5) +
      geom_ribbon(data = res,
                  aes(x = .data$year.pred,
                      ymin = .data$lower,
                      ymax = .data$upper),
                  alpha = 0.2, fill = 'blue') +
      geom_line(data = res,
                aes(x = .data$year.pred,
                    y = .data$median),
                size = 1, color = 'blue') +
      theme_classic() +
      xlab(NULL) +
      ylab('Abundance index')
  }
  if (type %in% c('transect_index', 'transect_trend', 'transect_trend+index')) {
    p <- p + facet_wrap(~.data$Transect)
  }
  print(p)
  return(p)
}
