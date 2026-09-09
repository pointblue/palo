#' Create input data for fitting a BBS-style hierarchical model for estimating
#' trends in point count survey data
#'
#' @param dat Dataframe created by \code{\link{summarize_PC_dat}}
#'
#' @return A list containing all inputs necessary for fitting BBS-style
#' hierarchical model.
#' 
#' @importFrom rlang .data
#' @export

setup_BBS_model <- function(dat) {
  sdat <- dat |> 
    dplyr::arrange(.data$Project, .data$Transect, .data$Point, .data$Year) |>
    dplyr::mutate(
      Project = factor(.data$Project, levels = unique(.data$Project)),
      Transect = factor(.data$Transect, levels = unique(.data$Transect)),
      Point = factor(.data$Point, levels = unique(.data$Point)))

  year.pred <- seq(min(sdat$Year), max(sdat$Year), 1)

  # proportion of transects within a project on which the species is present
  # in each year
  prop <- sdat |>
    dplyr::group_by(.data$Project, .data$Transect, .data$Year) |>
    dplyr::summarize(present = sum(.data$Count)) |>
    dplyr::mutate(present = ifelse(.data$present > 0, 1, 0)) |>
    dplyr::group_by(.data$Project, .data$Year) |>
    dplyr::summarize(n = length(.data$Transect),
                     n_present = sum(.data$present),
                     prop = .data$n_present / .data$n) |>
    dplyr::ungroup() |>
    dplyr::select(-'n', -'n_present') |>
    tidyr::pivot_wider(names_from = .data$Project, values_from = .data$prop,
                       fill = 0) |>
    dplyr::arrange(.data$Year) |>
    tibble::column_to_rownames('Year') |>
    as.matrix()

  # check number of detections overall within each project
  totals <- sdat |>
    dplyr::group_by(.data$Project) |>
    dplyr::summarize(total_count = sum(.data$Count))

  cat('Total detections by project:\n\n')
  print(totals)
  if (any(totals$total_count == 0)) {
    warning('One or more projects have zero detections of this species. For
            better model fit, rerun "summarize_PC_dat" and exclude this
            project.')
  } else if (any(totals$total_count < 100)) {
      warning('One or more projects have relatively few (<100) detections of
              this species. Model may have difficulty converging.')}

  list(
    observed = sdat$Count,
    zyear = sdat$Year - min(sdat$Year), #relative to baseline year
    year.pred = year.pred,
    zyear.pred = (year.pred - min(sdat$Year)), #values to predict for
    prop = prop,

    # project, transect, point, year ID numbers associated with each Count:
    project = sdat$Project |> as.numeric(),
    transect = sdat$Transect |> as.numeric(),
    point = sdat$Point |> as.numeric(),
    year = sdat$Year |> as.factor() |> as.numeric(),

    # number of unique IDs for each:
    nprojects = sdat$Project |> unique() |> length(),
    ntransects = sdat$Transect |> unique() |> length(),
    npoints = sdat$Point |> unique() |> length(),
    nyears = sdat$Year |> unique() |> length(),
    dat = sdat
  )
}
