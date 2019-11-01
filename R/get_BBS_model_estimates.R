#' Extract annual trend or abundance index estimates from BBS-style hierarchical model
#'
#' @param modresults Model output from running \code{\link{fit_BBS_model}}
#' @param type Type of estimates to extract: annual indices of abundance
#' ('index'), or predicted values from the linear trend ('trend').
#' @param projects Character string for the 4-letter project code(s), in all caps
#' @param years Numeric vector for 4-digit years. See details.
#'
#' @details The BBS model output does not explicitly keep track of which years
#' are included in each model, but keeps them in order. For type = 'index',
#' provide a vector of years that are actually included in the original data.
#' For type = 'trend', provide the full list of years between the earliest and last.
#' @return Dataframe containing predicted value and 95% credible interval for each year
#' @export
#'
get_BBS_model_estimates <- function(modresults, type, projects = NULL, years = NULL) {

  est <- MCMCvis::MCMCpstr(modresults, params = type,
                          func = function(x) {
                            quantile(x, probs = c(0.5, 0.025, 0.975))
                          })[[1]]

  if (type == 'slope') {
    dimnames(est)[[2]] <- c('median', 'lower', 'upper')
    res <- as.data.frame(est) %>%
      mutate(type = type)

  } else if (type %in% c('index', 'trend')) {
    dimnames(est)[[2]] <- projects

    res <- as.data.frame(est) %>%
      mutate(Year = years) %>%
      pivot_longer(-Year, names_to = 'variable', values_to = 'value') %>%
      separate(.data$variable, into = c('Project', 'variable'), sep = 4, fill = 'right') %>%
      mutate(variable = case_when(variable == '.50%' ~ 'median',
                                  variable == '.2.5%' ~ 'lower',
                                  variable == '.97.5%' ~ 'upper')) %>%
      pivot_wider(names_from = 'variable', values_from = 'value') %>%
      arrange(.data$Project, .data$Year) %>%
      mutate(type = type)
  }

  return(res)
}
