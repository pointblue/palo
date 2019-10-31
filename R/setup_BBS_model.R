#' Create input data for fitting a BBS-style hierarchical model for estimating
#' trends in point count survey data
#'
#' @param dat Dataframe created by \code{\link{summarize_PC_dat}}
#'
#' @return A list containing all inputs necessary for fitting BBS-style
#' hierarchical model.
#' @export
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#'
setup_BBS_model <- function(dat) {
  sdat <- dat %>%
    arrange(.data$Project, .data$Transect, .data$Point, .data$Year) %>%
    mutate(Project = factor(.data$Project, levels = unique(.data$Project)),
           Transect = factor(.data$Transect, levels = unique(.data$Transect)),
           Point = factor(.data$Point, levels = unique(.data$Point)))

  year.pred <- seq(min(sdat$Year), max(sdat$Year), 1)

  # proportion of transects within a project on which the species is present
  # in each year
  prop <- sdat %>%
    group_by(.data$Project, .data$Transect, .data$Year) %>%
    summarize(present = sum(.data$Count)) %>%
    mutate(present = ifelse(.data$present > 0, 1, 0)) %>%
    group_by(.data$Project, .data$Year) %>%
    summarize(n = length(.data$Transect),
              n_present = sum(.data$present),
              prop = .data$n_present / n) %>%
    ungroup() %>%
    select(-.data$n, -.data$n_present) %>%
    spread(key = .data$Project, value = .data$prop, fill = 0) %>%
    arrange(.data$Year) %>%
    column_to_rownames('Year') %>%
    as.matrix()

  # check number of detections overall within each project
  totals <- sdat %>%
    group_by(.data$Project) %>%
    summarize(total_count = sum(.data$Count))

  cat('Total detections by project:\n\n')
  print(totals)
  if (any(totals$total_count == 0)) {warning('One or more projects have zero detections of this species. For better model fit, rerun "summarize_PC_dat" and exclude this project.')}
  if (any(totals$total_count < 100)) {warning('One or more projects have relatively few detections of this species. Model may have difficulty converging.')}

  list(
    observed = sdat %>% pull(Count),
    zyear = sdat$Year - min(sdat$Year), #relative to baseline year
    year.pred = year.pred,
    zyear.pred = (year.pred - min(sdat$Year)), #values to predict for
    prop = prop,

    # project, transect, point, year ID numbers associated with each Count:
    project = sdat$Project %>% as.numeric(),
    transect = sdat$Transect %>% as.numeric(),
    point = sdat$Point %>% as.numeric(),
    year = sdat$Year %>% as.factor() %>% as.numeric(),

    # number of unique IDs for each:
    nprojects = sdat$Project %>% unique() %>% length(),
    ntransects = sdat$Transect %>% unique() %>% length(),
    npoints = sdat$Point %>% unique() %>% length(),
    nyears = sdat$Year %>% unique() %>% length(),
    dat = sdat
  )
}
