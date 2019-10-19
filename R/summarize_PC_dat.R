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
#'
summarize_PC_dat <- function(df, species, project) {

  # format data to include surveys where count is zero,
  # summarize total Count at each unique visit in case they were in separate
  # bins, then calculate max, mean, and se of counts across visits in each
  # point/year combination
  df %>%
    format_PC_dat(species, project) %>%
    group_by(Spp, Project, Transect, Point, Year, Visit) %>%
    summarize(Count = sum(Count)) %>%
    ungroup() %>%
    group_by(Spp, Project, Transect, Point, Year) %>%
    summarize(n_visits = length(Count),
              count_max = max(Count),##
              count_mean = mean(Count),
              count_se = sd(Count)/sqrt(n_visits))
}

format_PC_dat <- function(df, species, project) {
  # filter projects, add Year field & select relevant columns
  df <- df %>%
    filter(Project == project) %>%
    mutate(Year = as.numeric(format(Date, '%Y'))) %>%
    select(Project, Transect, Point, Year, Visit, Spp,
           Count, mindist, maxdist)

  # generate complete list of all unique surveys in the Project (regardless of
  # species), then join to subset containing species of interest; ensure data is
  # 'complete' with all surveys listed even if count is zero
  df %>%
    select(Project, Transect, Point, Year, Visit) %>%
    distinct() %>%
    full_join(df %>% filter(Spp == species),
              by = c('Project', 'Transect', 'Point', 'Year', 'Visit')) %>%
    complete(Spp, nesting(Project, Transect, Point, Year, Visit)) %>%
    filter(!is.na(Spp)) %>%
    mutate(Count = case_when(is.na(Count) ~ 0,
                             TRUE ~ Count)) %>%
    mutate_at(vars(Project:Point, Spp), as.factor)
}
