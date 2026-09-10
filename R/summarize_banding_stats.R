#' Summarize banding statistics
#'
#' @description `summarize_capture_stats()` returns total captures by location,
#' species, and either month, year, or season
#'
#' `summarize_nethrs_stats()` returns total net hours by location, and either
#' month, year or season
#'
#' @param df optional data frame representing banding or net hours database
#' @param path optional filepath to local copy of banding or net hours database
#' @param location optional character string of banding station codes to limit
#'   results
#' @param datemin optional minimum date in YYYY-mm-dd format
#' @param datemax optional maximum date in YYYY-mm-dd format
#' @param species optional character string of species codes to include (only
#'   applies to banding data)
#' @param by how statistics should be summarized, either "month", "year", or
#'   "season"; see Details
#' @param winter_adjust logical; see Details
#'
#' @details These functions are designed to work directly with standard
#'   Palomarin .dbf files, e.g. "allnumb.dbf", "band.dbf", or
#'   "allpalonthrs.dbf". If `df` is provided, it should represent a local
#'   dataframe version of the banding or net hours database already in memory,
#'   potentially already filtered to the desired locations, date range, and/or
#'   species, or with a custom "season" column added. If `df` is not provided,
#'   instead use `path` to point to a .dbf file that will be read in using
#'   [foreign::read.dbf()].
#'
#'   The `location`, `datemin`, `datemax`, and `species` arguments are optional
#'   and provide convenient ways to filter the data as needed, especially if
#'   reading in the .dbf directly from the `path` provided. Note that `species`
#'   filters apply only to the banding data.
#'
#'   After filtering, the internal [summarize_by()] function is called to either
#'   count the number of captures or sum the number of net hours as appropriate,
#'   summarized by location, species (for banding data only), and the desired
#'   time frame specified using the `by` argument. Currently supported values
#'   include "month", "season", or "year". If `by = "month"`, capture statistics
#'   will be calculated by location, species (for banding data only), year, and
#'   month. If `by = "season"`, summary statistics will be calculated by
#'   location, species (for banding data only), year, and season. If `by =
#'   "year"`, capture statistics will be calculated by location, species (for
#'   banding data only), and year.
#'
#'   If the corresponding fields "month", "season", or "year" are already
#'   present in the `df` provided, those are the fields that will be used to
#'   calculate these summary statistics, and thereby provide a means to
#'   customize how month, year, or season are defined. If those fields are not
#'   present in `df`, or if `path` is provided, the corresponding month and year
#'   will be calculated from the "DATE" field, and default seasons will be
#'   defined as: spring (March through July), fall (August through October), and
#'   winter (November through February). Note that in this case, there must be
#'   "DATE" field present in `df` or the database to which `path` points.
#'
#'   If `by = 'season'`, with `winter_adjust = TRUE` (the default), there will
#'   be an additional field created called "year_adjust". This field will be
#'   identical to "year" except that when `season == "winter"` and `month < 6`,
#'   year_adjust will be equal to year - 1, so that, for example, January will
#'   be grouped with December from the previous calendar year. In this case,
#'   summary statistics will be calculated by location, species (for banding
#'   data only), and year_adjust.
#'
#'   Finally, the resulting dataframe will be expanded to sure there is a unique
#'   row for every location, year, month, season, and/or species, including
#'   where and when there may have been zero captures or zero net hours in the
#'   data.
#'
#' @returns tibble
#' @export
#'
#' @examples
#'
#' data(sample_band)
#' capture_MUHO = summarize_capture_stats(
#'    df = sample_band, location = 'MUHO', species = c('SWTH', 'AMGO'),
#'    datemin ='2010-01-01', by = 'season', winter_adjust = FALSE)
#'
#' data(sample_nethrs)
#' effort_PGUP = summarize_nethrs_stats(
#'    df = sample_nethrs, location = 'PGUP', datemax ='2010-10-15',
#'    by = 'month') 

#'\dontrun{
#'# # summarize seasonal capture statistics for SOSP and WIWA directly from
# # allnumb, with filters for location and minimum date, and using default
# # season definitions:
# captures_season_default = summarize_capture_stats(
#    path = 'path/to/allnumb.dbf',
#    location = c("PN", "PGUP"),
#    datemin = "1979-03-01",
#    species = c("SOSP", "WIWA"),
#    by = "season",
#    winter_adjust = TRUE)
# 
# # this is an equivalent approach, but reading in the data first separately,
# # to filter the data and define custom seasons:
# band = foreign::read_dbf('path/to/allnumb.dbf') |>
#    dplyr::filter(SPEC %in% c("SOSP", "WIWA") & DATE >= "1979-03-01") |>
#    dplyr::mutate(
#      # define custom seasons based on yday (day of year); note these will
#      # all be off by one day in leap years
#      yday = format(DATE, "%j"), # day of year
#      season = case_when(
#        yday >= 74 & yday <= 212 ~ "spring", # Mar 15--Jul 31
#        yday >= 213 & yday <= 319 ~ "fall", # Aug 1--Nov 15
#        yday >= 320 | yday <= 73 ~ "winter", # Nov 16 -- Mar 14
#      )
#    )
# captures_season_custom = summarize_capture_stats(
#    df = band, by = "season", winter_adjust = TRUE)
# 
# # summarizing net hours works the same way:
# effort = summarize_nethrs_stats(
#    path = 'path/to/nethours.dbf',
#    location = c("PN", "PGUP"),
#    by = "season")
#'}


summarize_capture_stats = function(path = NULL, df = NULL, location = NULL,
                                   datemin = NULL, datemax = NULL, 
                                   species = NULL, by = 'year',
                                   winter_adjust = TRUE) {
  
  if (is.null(df) & !is.null(path)) {
    df = foreign::read.dbf(path) |> 
      dplyr::select(.data$INITIALS:.data$COM) # (drop extra columns labeled "X")
  } else if (is.null(df) & is.null(path)) {
    stop('no input dataframe or path to a .dbf provided')
  }
  
  # optional filtering
  if (!is.null(location)) {
    df = df |> dplyr::filter(.data$LOC %in% location)
  }
  if (!is.null(datemin)) {
    df = df |> dplyr::filter(.data$DATE >= datemin)
  }
  if (!is.null(datemax)) {
    df = df |> dplyr::filter(.data$DATE <= datemax)
  }
  if (!is.null(species)) {
    df = df |> dplyr::filter(.data$SPEC %in% species)
  }
  
  summarize_by(df, by = by, stat = 'captures', 
               winter_adjust = winter_adjust)
  
  
}