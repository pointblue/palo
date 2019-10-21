#' Extract annual trend or abundance index estimates from BBS-style hierarchical model
#'
#' @param modresults Model output from running \code{\link{fit_BBS_model}}
#' @param inputdat Input data from running \code{\link{setup_BBS_model}}
#' @param type Type of estimates to extract: annual indices of abundance
#' overall or by transect ('global_index', 'transect_index'), or
#' predicted values from the linear trend overall or by transect ('global_trend',
#' 'transect_trend').
#' @details This function is automatically called by \code{\link{plot_BBS_model}}, but may be useful for generating customized plots.
#'
#' @return Dataframe containing predicted value and HDI for each year present in inputdat$year.pred (created by \code{\link{setup_BBS_model}})
#' @export
#'
get_BBS_model_estimates <- function(modresults, inputdat, type) {
  if (type %in% c('global_trend', 'global_index')) {
    if (type == 'global_trend') {param = 'count_pred'}
    if (type == 'global_index') {param = 'abund'}
    ci <- MCMCvis::MCMCpstr(modresults, params = param,
                            func = function(x) HDInterval::hdi(x, .95))[[1]]
    med <- MCMCvis::MCMCpstr(modresults, params = param, func = median)[[1]]
    res <- cbind(year.pred = inputdat$year.pred,
                 median = med,
                 data.frame(ci))
  } else if (type %in% c('transect_trend', 'transect_index')) {
    if (type == 'transect_trend') {param = 'count_pred_transect'}
    if (type == 'transect_index') {param = 'abund_transect'}
    ci <- MCMCvis::MCMCpstr(modresults, params = param,
                            func = function(x) HDInterval::hdi(x, .95))[[1]]
    dimnames(ci)[[2]] <- levels(as.factor(as.character(inputdat$dat$Transect)))

    med <- MCMCvis::MCMCpstr(modresults, params = param, func = median)[[1]]
    dimnames(med)[[2]] <- levels(as.factor(as.character(inputdat$dat$Transect)))

    res <- cbind(year.pred = inputdat$year.pred,
                 data.frame(med),
                 data.frame(ci)) %>%
      gather(-.data$year.pred, key = 'variable', value = 'value') %>%
      separate(.data$variable, into = c('Transect', 'variable'), fill = 'right') %>%
      mutate(variable = replace_na(.data$variable, 'median')) %>%
      spread(key = .data$variable, value = .data$value) %>%
      arrange(.data$Transect, .data$year.pred)
  }
  return(res)
}
