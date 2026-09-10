#' Fill missing values
#'
#' Internal utility function called by [summarize_capture_stats()] and
#' [summarize_nethrs_stats()]; not intended to be called directly
#'
#' @param df dataframe or tibble
#'
#' @details Utility function to expand a dataframe to ensure there is a row for
#' every location, year, month, and/or season present in the data, as well as
#' every species (in the case of banding data). For example, used by
#' [summarize_capture_stats()] to ensure that a lack of captures for a given
#' species in a given time frame (month, season, or year) is represented with a
#' 0 value. For net hours, to ensure each location is represented in each month,
#' season, or year, even if there was no effort during that time frame. The
#' zeroes are by default filled into a column labeled "y".
#'
#' @returns tibble
#' @export
#' @keywords internal

fill_missing = function(df) {
  candidate_cols = c('SPEC', 'LOC', 'LOCATION', 'year', 'year_adjust', 
                     'month', 'season')
  
  present_cols = dplyr::intersect(candidate_cols, names(df))
  
  df |> 
    # drop any factors already filtered out
    dplyr::mutate(
      dplyr::across(
        dplyr::any_of(c('SPEC', 'LOC', 'LOCATION')), droplevels) 
    ) |> 
    tidyr::complete(
      !!!rlang::syms(present_cols),
      fill = list(y = 0))
}

