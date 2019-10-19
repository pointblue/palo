#' Plot trend results from BBS-style hierarchical model for point count survey data
#'
#' @param modresults model output from running \code{\link[MASS]{fit_BBS_model}}
#' @param inputdat input data from running \code{\link[MASS]{setup_BBS_model}}
#' @param type type of plot to produce: plots of annual indices of abundance
#' overall or by transect ('annual_overall', 'annual_transect'), plots of
#' linear trend lines overall or by transect ('smooth_overall',
#' 'smooth_transect'), or combinations of both ('annual+smooth_overall',
#' 'annual+smooth_transect'). Defaults to 'annual+smooth_overall', similar to
#' BBS output.
#' @param plot defaults to TRUE, to display a plot of results
#'
#' @return returns ggplot
#' @export
#' @import ggplot2 MCMCvis
#' @importFrom dplyr %>%
#' @importFrom HDInterval hdi
#'
#'
plot_BBS_model <- function(modresults, inputdat,
                           type = 'annual+smooth_overall') {
  if (type == 'annual_overall') {
    res <- get_BBS_estimates(modresults, inputdat, 'abund', 'overall')
  } else if (type == 'smooth_overall') {
    res <- get_BBS_estimates(modresults, inputdat, 'count_pred', 'overall')
  } else if (type == 'annual_transect') {
    res <- get_BBS_estimates(modresults, inputdat, 'abund_transect', 'transect')
  } else if (type == 'smooth_transect') {
    res <- get_BBS_estimates(modresults, inputdat, 'count_pred_transect', 'transect')
  } else if (type == 'annual+smooth_overall') {
    res <- get_BBS_estimates(modresults, inputdat, 'abund', 'overall')
    res2 <- get_BBS_estimates(modresults, inputdat, 'count_pred', 'overall')
  } else if (type == 'annual+smooth_transect') {
    res <- get_BBS_estimates(modresults, inputdat, 'abund_transect', 'transect')
    res2 <- get_BBS_estimates(modresults, inputdat, 'count_pred_transect', 'transect')
  }

  if (type %in% c('annual+smooth_overall', 'annual+smooth_transect')) {
    p <- ggplot(res) +
      geom_pointrange(aes(x = year.pred, y = median,
                          ymin = lower, ymax = upper)) +
      geom_line(aes(x = year.pred, y = median), size = 1) +
      geom_ribbon(data = res2,
                  aes(x = year.pred, ymin = lower, ymax = upper),
                  alpha = 0.2, fill = 'blue') +
      geom_line(data = res2,
                aes(x = year.pred, y = median, color = Transect),
                size = 1, color = 'blue') +
      theme_classic() +
      xlab(NULL) +
      ylab('Abundance index')

  } else {
    obsdat <- data.frame(Year = inputdat$dat$Year,
                         Transect = inputdat$dat$Transect,
                         Count = inputdat$observed)
    p <- ggplot(obsdat) +
      geom_jitter(aes(Year, Count), width = 0.1, height = 0.1,
                  color = 'gray30', shape = 21, alpha = 0.5) +
      geom_ribbon(data = res,
                  aes(x = year.pred, ymin = lower, ymax = upper),
                  alpha = 0.2, fill = 'blue') +
      geom_line(data = res,
                aes(x = year.pred, y = median),
                size = 1, color = 'blue') +
      theme_classic() +
      xlab(NULL) +
      ylab('Abundance index')
  }
  if (type %in% c('annual_transect', 'smooth_transect', 'annual+smooth_transect')) {
    p <- p + facet_wrap(~Transect)
  }
  print(p)
  return(p)
}

get_BBS_estimates <- function(modresults, inputdat, param, type) {
  if (type == 'overall') {
    ci <- MCMCpstr(modresults, params = param,
                   func = function(x) HDInterval::hdi(x, .95))[[1]]
    med <- MCMCpstr(modresults, params = param, func = median)[[1]]
    res <- cbind(year.pred = inputdat$year.pred,
                 median = med,
                 data.frame(ci))
  } else if (type == 'transect') {
    ci <- MCMCpstr(modresults, params = param,
                   func = function(x) HDInterval::hdi(x, .95))[[1]]
    dimnames(ci)[[2]] <- levels(as.factor(as.character(inputdat$dat$Transect)))

    med <- MCMCpstr(modresults, params = param, func = median)[[1]]
    dimnames(med)[[2]] <- levels(as.factor(as.character(inputdat$dat$Transect)))

    res <- cbind(year.pred = inputdat$year.pred, data.frame(med), data.frame(ci)) %>%
      gather(-year.pred, key = 'variable', value = 'value') %>%
      separate(variable, into = c('Transect', 'variable'), fill = 'right') %>%
      mutate(variable = replace_na(variable, 'median')) %>%
      spread(key = variable, value = value) %>%
      arrange(Transect, year.pred)
  }
  return(res)
}
