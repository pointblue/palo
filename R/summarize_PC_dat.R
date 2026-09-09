#' Format and summarize point count survey data for a given species and project code.
#' 
#' Returns the total count (regardless of detection distance) for the selected species at each unique combination of Visit, Point, Transect, and Year in the data, including zero counts where the species was not detected.
#'
#' @param df Data frame of point count data, with the standard field names as
#' in the format downloaded from CADC Project Leader. Filter to include only
#' the detection distances that should be included in the total count.
#' @param species Character string for the 4-letter species code, in all caps
#' @param project Character string for the 4-letter project code(s), in all caps
#'
#' @return Dataframe 
#' 
#' @importFrom rlang .data
#' @export
#'

summarize_PC_dat <- function(df, species, project) {

  # format data to include surveys where count is zero,
  # summarize total Count at each unique visit across distance bins
  df |>
    format_PC_dat(species, project) |>
    dplyr::group_by(.data$Spp, .data$Project, .data$Transect, .data$Point, 
                    .data$Year, .data$Visit) |>
    dplyr::summarize(Count = sum(.data$Count)) |>
    dplyr::ungroup()
}

format_PC_dat <- function(df, species, project) {
  # filter projects, add Year field & select relevant columns
  sdf <- df |>
    dplyr::filter(.data$Project %in% project) |>
    dplyr::mutate(Year = as.numeric(format(.data$Date, '%Y')),
                  Spp = as.character(.data$Spp)) |>
    dplyr::select(
      dplyr::any_of(c('Project', 'Transect', 'Point', 'Year', 'Visit', 'Spp',
                      'Count', 'mindist', 'maxdist')))

  # generate complete list of all unique surveys in the Project (regardless of
  # species), then join to subset containing species of interest; ensure data is
  # 'complete' with all surveys listed even if count is zero
  sdf |>
    dplyr::select(
      dplyr::any_of(c('Project', 'Transect', 'Point', 'Year', 'Visit'))) |>
    dplyr::distinct() |>
    dplyr::left_join(
      sdf |> dplyr::filter(.data$Spp == species),
      by = c('Project', 'Transect', 'Point', 'Year', 'Visit')) |>
    dplyr::mutate(Spp = tidyr::replace_na(.data$Spp, species),
                  Spp = as.factor(.data$Spp),
                  Count = dplyr::if_else(is.na(.data$Count), 0, .data$Count))
}
