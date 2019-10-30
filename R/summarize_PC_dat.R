#' Format and summarize point count survey data for a given species and project code.
#'
#' @param df Data frame of point count data, with the standard field names as
#' in the format downloaded from CADC Project Leader. Filter to include only
#' the detection distances that should be included in the total count.
#' @param species Character string for the 4-letter species code, in all caps
#' @param project Character string for the 4-letter project code(s), in all caps
#'
#' @return Dataframe containing the total count (regardless of detection
#' distance) for the selected species at each unique combination of Visit,
#' Point, Transect, and Year in the data, including zero counts where the
#' species was not detected.
#' @export
#' @import dplyr
#' @importFrom rlang .data
#' @importFrom stats sd
#'
summarize_PC_dat <- function(df, species, project) {

  # format data to include surveys where count is zero,
  # summarize total Count at each unique visit across distance bins
  df %>%
    format_PC_dat(species, project) %>%
    group_by(.data$Spp, .data$Project, .data$Transect, .data$Point, .data$Year,
             .data$Visit) %>%
    summarize(Count = sum(.data$Count)) %>%
    ungroup()
}

format_PC_dat <- function(df, species, project) {
  # filter projects, add Year field & select relevant columns
  df <- df %>%
    filter(.data$Project %in% project) %>%
    mutate(Year = as.numeric(format(.data$Date, '%Y'))) %>%
    select(.data$Project, .data$Transect, .data$Point, .data$Year, .data$Visit,
           .data$Spp, .data$Count, .data$mindist, .data$maxdist)

  # generate complete list of all unique surveys in the Project (regardless of
  # species), then join to subset containing species of interest; ensure data is
  # 'complete' with all surveys listed even if count is zero
  df %>%
    select(.data$Project, .data$Transect, .data$Point, .data$Year, .data$Visit) %>%
    distinct() %>%
    full_join(df %>% filter(.data$Spp == species),
              by = c('Project', 'Transect', 'Point', 'Year', 'Visit')) %>%
    complete(.data$Spp,
             nesting(.data$Project, .data$Transect, .data$Point, .data$Year, .data$Visit)) %>%
    filter(!is.na(.data$Spp)) %>%
    mutate(Count = case_when(is.na(.data$Count) ~ 0,
                             TRUE ~ .data$Count)) %>%
    mutate_at(vars(.data$Project:.data$Point, .data$Spp), as.factor)
}
