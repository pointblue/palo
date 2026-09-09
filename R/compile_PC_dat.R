#' Compile point count survey data from a directory
#'
#' @param dir Directory containing point count survey data downloaded from CADC
#' @param pattern Regex pattern for identifying files to include (see [base::list.files()])
#'
#' @return Dataframe of compiled data, excluding any rows with all NA values.
#' @details Also tests for missing values in critical fields.
#' @export

compile_PC_dat <- function(dir, pattern = '.csv') {
  # get names of all files in directory that match pattern
  file_list <- list.files(dir, pattern = pattern, full.names = TRUE)

  # compile master data set (assuming CADC Project Leader format)
  dat <- purrr::map_dfr(
    file_list,
    ~readr::read_csv(.x, col_types = 'cccccnDtttncnccccccclcccc')) |> 
    dplyr::mutate(
      dplyr::across(
        dplyr::any_of(
          c('Project', 'Study Area', 'Transect', 'Point', 'Protocol', 
            'Time Bin', 'Spp', 'Common Name', 'Scientific Name', 
            'Detection Cue', 'Distance Bin ID', 'Distance Bin', 
            'Breeding Status', 'Researcher', 'Data Status')), as.factor))

  # check for NA values in critical fields:
  critfields <- dat |> 
    dplyr::select(
      dplyr::any_of(c('Project', 'Study Area', 'Transect', 'Point', 'Protocol', 
                    'Spp', 'Distance_Bin_ID', 'Count')))
  rows_with_na = !stats::complete.cases(critfields)
  
  if(length(rows_with_na) > 0) {
    warning(length(rows_with_na), 
            ' rows contain NA values in one or more critical fields')
  }
  return(dat)
}
