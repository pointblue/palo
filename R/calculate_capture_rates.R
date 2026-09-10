#' Calculate capture rates from banding data
#'
#' Capture rates per 1000 net-hours, by year, season, month, location, and/or
#' species
#'
#' @param captures dataframe; see Details
#' @param effort dataframe; see Details
#' @param by_species logical; if TRUE, summarize capture rates by species
#' @param add_annual logical; if TRUE, append capture rates by year in addition
#'   to stats by season or month
#'
#' @details This function is designed to work with the output of the
#'   [summarize_capture_stats()] and [summarize_nethrs_stats()] functions, which
#'   should first be used to filter and summarize banding and nethours data by
#'   month, year, or season. The outputs of these functions become the inputs
#'   `effort` and `captures`. This function then joints the two data sets to
#'   calculate capture rates.
#'
#'   If `by_species = TRUE`, capture rates are calculated separately for each
#'   species in `captures`; if `by_species = FALSE`, the overall capture rate
#'   for all species in `captures` will instead by calculated.
#'
#'   If `add_annual = TRUE`, overall capture rates for the entire year will be
#'   appended to capture rates calculated by season or month. This argument is
#'   only relevant if `captures` and `effort` data sets contain the fields
#'   'season' or 'month'. Note that with seasonal statistics, winter seasons are
#'   often split across two calendar years; this is already handled by
#'   [summarize_capture_stats()] and [summarize_nethrs_stats()] with the
#'   argument `winter_adjustment = TRUE`. In that case, annual capture rate
#'   statistics will instead by calculated on the "year_season" field, created
#'   to keep entire winter seasons together. (For example, this means that if
#'   "winter" is considered to be December-February, then Jan and Feb of 2026
#'   will be assigned a year_season of 2025 to keep it grouped with December
#'   2025, and annual statistics will run from March through February rather
#'   than January through December.)
#'
#' @returns tibble
#' @export
#'
#' @examples
#' data(sample_band)
#' capturedat = summarize_capture_stats(
#'    df = sample_band, location = c('MUHO', 'RECR'),
#'    species = c('SWTH', 'AMGO'), by = 'season', winter_adjust = FALSE)
#'
#' data(sample_nethrs)
#' effortdat = summarize_nethrs_stats(
#'    df = sample_nethrs, location = c('MUHO', 'RECR'),
#'    by = 'season', winter_adjust = FALSE)
#'
#' # make sure capturedat and effortdat include the same locations, date ranges,
#' # and groupings (e.g. by season or month)
#' capture_stats = calculate_capture_rates(
#'    captures = capturedat, effort = effortdat, by_species = TRUE,
#'    add_annual = TRUE)
#' 

calculate_capture_rates = function(captures, effort, by_species = TRUE,
                                   add_annual = TRUE) {

  # join monthly capture totals and net hour totals by location:
  dat = dplyr::full_join(captures, 
                         effort |> dplyr::rename(LOC = .data$LOCATION))

  if (by_species) {
    
    if (add_annual & ('month' %in% names(dat) | 'season' %in% names(dat))) {
      # also add an annual total
      cols = c('SPEC', 'LOC', 'year', 'year_season')
      dat_annual = calculate_totals(dat, grouping_cols = cols)
      dat_sum = dplyr::bind_rows(dat, dat_annual)
    } else {
      dat_sum = dat
    }
    
  } else {
    # first summarize overall total annual/seasonal capture rates across all
    # species in captures data
    
    cols = c('LOC', 'year', 'year_season', 'season', 'month')
    dat_sum = calculate_totals(dat, grouping_cols = cols)
    
    if (add_annual & ('month' %in% names(dat) | 'season' %in% names(dat))) {
      # also add an annual total
      cols = c('LOC', 'year', 'year_season')
      dat_annual = calculate_totals(dat, grouping_cols = cols)
      dat_sum = dplyr::bind_rows(dat_sum, dat_annual)
      
    }

  }
  
  # combine and calculate capture rate per 1000 net hours:
  res = dat_sum |>
    dplyr::mutate(capture_rate = .data$captures / .data$nethours * 1000)
  
  res |> 
    dplyr::arrange(
      dplyr::pick(
        dplyr::any_of(c('SPEC', 'LOC', 'year', 'year_season', 'season', 
                        'month')))
    )
}

calculate_totals = function(df, grouping_cols) {
  # summarize by the columns defined above, depending on selected "stat"
  df = df |> 
    dplyr::group_by(
      dplyr::across(
        dplyr::any_of(grouping_cols))) |> 
    dplyr::summarize(
      captures = sum(.data$captures),
      nethours = sum(.data$nethours),
      .groups = 'drop')
}


