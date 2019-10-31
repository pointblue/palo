#' Extract annual trend or abundance index estimates from BBS-style hierarchical model
#'
#' @param modresults Model output from running \code{\link{fit_BBS_model}}
#' @param inputdat Input data from running \code{\link{setup_BBS_model}}
#' @param type Type of estimates to extract: annual indices of abundance
#' ('index'), or predicted values from the linear trend ('trend').
#' @details This function is automatically called by \code{\link{plot_BBS_model}},
#' but may be useful for generating customized plots.
#'
#' @return Dataframe containing predicted value and HDI for each year present
#' in inputdat$year.pred (created by \code{\link{setup_BBS_model}})
#' @export
#'
get_BBS_model_estimates <- function(modresults, inputdat, type) {

  ci <- MCMCvis::MCMCpstr(modresults, params = type,
                          func = function(x) HDInterval::hdi(x, .95))[[1]]
  dimnames(ci)[[2]] <- levels(as.factor(as.character(inputdat$dat$Project)))

  med <- MCMCvis::MCMCpstr(modresults, params = type, func = median)[[1]]
  names <- levels(as.factor(as.character(inputdat$dat$Project)))
  if (length(names) == 1) {names <- list(names)}
  dimnames(med)[[2]] <- names

  res <- cbind(year.pred = inputdat$year.pred,
               data.frame(med),
               data.frame(ci)) %>%
    gather(-.data$year.pred, key = 'variable', value = 'value') %>%
    separate(.data$variable, into = c('Project', 'variable'), fill = 'right') %>%
    mutate(variable = replace_na(.data$variable, 'median')) %>%
    spread(key = .data$variable, value = .data$value) %>%
    arrange(.data$Project, .data$year.pred)

  return(res)
}
