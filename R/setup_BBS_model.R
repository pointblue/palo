#' Create input data for fitting BBS-style hierarchical model for estimating trends in point count survey data
#'
#' @param dat Dataframe or tibble containing point count data
#' @param count Option allowing users to choose between modeling the maximum
#' count on any visit to a point in a given year ('count_max'), or the mean
#' across all visits ('count_mean').
#'
#' @return A list containing all inputs necessary for fitting BBS-style hierarchical model.
#' @export
#' @importFrom dplyr %>%
#'
setup_BBS_model <- function(dat, count = 'count_max') {
  sdat <- dat %>% arrange(Transect, Point, Year)
  year.pred <- seq(min(sdat$Year), max(sdat$Year), 1)
  prop <- sdat %>%
    group_by(Transect, Year) %>%
    summarize(present = sum(count_mean)) %>%
    mutate(present = ifelse(present > 0, 1, 0)) %>%
    group_by(Year) %>%
    summarize(n = length(Transect),
              n_present = sum(present),
              prop = n_present / n) %>%
    ungroup()

  list(
    observed = sdat %>% pull(!!count),
    zyear = sdat$Year - min(sdat$Year), #relative to baseline year
    year.pred = year.pred,
    zyear.pred = (year.pred - min(sdat$Year)), #values to predict for
    prop = prop %>% pull(prop),

    # transect, point, year ID numbers associated with each Count:
    transect = sdat %>% pull(Transect) %>% as.numeric(),
    point = sdat %>% pull(Point) %>% as.numeric(),
    year = sdat %>% pull(Year) %>% as.factor() %>% as.numeric(),

    # number of unique IDs for each:
    ntransects = length(unique(sdat$Transect)),
    npoints = length(unique(sdat$Point)),
    nyears = length(unique(sdat$Year)),
    dat = sdat
  )
}
