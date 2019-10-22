#' Format and summarize point count survey data for a given species and project code.
#'
#' @param df master point count data set, with dates, transects, and distances filtered
#' @param species 4-letter species code
#' @param project 4-letter project code
#'
#' @return dataframe containing all Points, Transects, and Years in the data,
#' the number of visits, the max and mean count for the species selected,
#' including zero counts where the species was not detected, and the standard
#' error of the mean count
#' @export
#' @import dplyr
#' @importFrom rlang .data
#' @importFrom stats sd
#'
summarize_PC_dat <- function(df, species, project) {

  # format data to include surveys where count is zero,
  # summarize total Count at each unique visit in case they were in separate
  # bins, then calculate max, mean, and se of counts across visits in each
  # point/year combination
  df %>%
    format_PC_dat(species, project) %>%
    group_by(.data$Spp, .data$Project, .data$Transect, .data$Point, .data$Year,
             .data$Visit) %>%
    summarize(Count = sum(.data$Count)) %>%
    ungroup() %>%
    group_by(.data$Spp, .data$Project, .data$Transect, .data$Point, .data$Year) %>%
    summarize(n_visits = length(.data$Count),
              count_max = max(.data$Count),
              count_mean = mean(.data$Count),
              count_se = sd(.data$Count)/sqrt(.data$n_visits)) %>%
    ungroup()
}

format_PC_dat <- function(df, species, project) {
  # filter projects, add Year field & select relevant columns
  df <- df %>%
    filter(.data$Project == project) %>%
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
