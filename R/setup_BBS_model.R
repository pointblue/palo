#' Create input data for fitting BBS-style hierarchical model for estimating trends in point count survey data
#'
#' @param dat Dataframe or tibble containing point count data
#' @param count Character string giving the name of the field containing the count data to be modeled.
#'
#' @return A list containing all inputs necessary for fitting BBS-style hierarchical model.
#' @export
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#'
setup_BBS_model <- function(dat, count = 'count_max') {
  sdat <- dat %>% arrange(.data$Transect, .data$Point, .data$Year) %>%
    mutate(Transect = factor(.data$Transect, levels = unique(.data$Transect)),
           Point = factor(.data$Point, levels = unique(.data$Point)))

  year.pred <- seq(min(sdat$Year), max(sdat$Year), 1)
  prop <- sdat %>%
    group_by(.data$Transect, .data$Year) %>%
    summarize(present = sum(.data$count_mean)) %>%
    mutate(present = ifelse(.data$present > 0, 1, 0)) %>%
    group_by(.data$Year) %>%
    summarize(n = length(.data$Transect),
              n_present = sum(.data$present),
              prop = .data$n_present / n) %>%
    ungroup()

  list(
    observed = sdat %>% pull(!!count),
    zyear = sdat$Year - min(sdat$Year), #relative to baseline year
    year.pred = year.pred,
    zyear.pred = (year.pred - min(sdat$Year)), #values to predict for
    prop = prop %>% pull(prop),
    effort = sdat$n_visits,

    # transect, point, year ID numbers associated with each Count:
    transect = sdat$Transect %>% as.numeric(),
    point = sdat$Point %>% as.numeric(),
    year = sdat$Year %>% as.factor() %>% as.numeric(),

    # number of unique IDs for each:
    ntransects = sdat$Transect %>% unique() %>% length(),
    npoints = sdat$Point %>% unique() %>% length(),
    nyears = sdat$Year %>% unique() %>% length(),
    dat = sdat
  )
}
