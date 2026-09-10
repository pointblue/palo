library(palo)
library(dplyr)

data(sample_nethrs)
data(sample_band)

# NET HOURS TESTING-------------

# net hours by season, with winter adjustment:
nethrs_season_adjust = palo:::summarize_by(
  df = sample_nethrs, stat = 'nethours', by = 'season', winter_adjust = TRUE)
nethrs_season_adjust_loc = nethrs_season_adjust |> group_by(LOCATION) |> count()
nethrs_season_adjust_seas = nethrs_season_adjust |> group_by(season) |> count()

# net hours by season, without winter adjustment (by calendar year, so for
# example this joins Jan 2026 with December 2026)
nethrs_season_calendar = palo:::summarize_by(
  df = sample_nethrs, stat = 'nethours', by = 'season', winter_adjust = FALSE)
nethrs_season_calendar_loc = nethrs_season_calendar |> 
  group_by(LOCATION) |> count()
nethrs_season_calendar_seas = nethrs_season_calendar |> 
  group_by(season) |> count()

testthat::test_that('net hours by season works as expected', {
  testthat::expect_true(
    all(c('fall', 'spring', 'winter') %in% nethrs_season_adjust$season))
  testthat::expect_true('year_adjust' %in% names(nethrs_season_adjust))
  testthat::expect_true(
    all(c(2008:2012) %in% nethrs_season_adjust$year_adjust)) 
  # because Jan-Feb 2009 count as 2008 winter season
  testthat::expect_true(all(nethrs_season_adjust_loc$n == 15))
  testthat::expect_true(all(nethrs_season_adjust_seas$n == 50)) 
  
  testthat::expect_true(
    all(c('fall', 'spring', 'winter') %in% nethrs_season_calendar$season))
  testthat::expect_true('year' %in% names(nethrs_season_calendar))
  testthat::expect_false(2008 %in% nethrs_season_calendar$year)
  testthat::expect_true(all(nethrs_season_calendar_loc$n == 12)) 
  testthat::expect_true(all(nethrs_season_calendar_seas$n == 40)) 
})


nethrs_month = palo:::summarize_by(
  df = sample_nethrs, stat = 'nethours', by = 'month')
nethrs_month_loc = nethrs_month |> group_by(LOCATION) |> count()
nethrs_month_yr = nethrs_month |> group_by(year) |> count()
testthat::test_that('net hours by month works as expected', {
  testthat::expect_true(
    all(c(1:12) %in% nethrs_month$month))
  testthat::expect_true('year' %in% names(nethrs_month))
  testthat::expect_true(all(nethrs_month_loc$n == 48)) # 4 years x 12months
  testthat::expect_true(all(nethrs_month_yr$n == 120)) # 10 locations x 12months
})

nethrs_year = palo:::summarize_by(
  df = sample_nethrs, stat = 'nethours', by = 'year')
nethrs_year_loc = nethrs_year |> group_by(LOCATION) |> count()
nethrs_year_yr = nethrs_year |> group_by(year) |> count()
testthat::test_that('net hours by year works as expected', {
  testthat::expect_true(
    all(c(2009:2012) %in% nethrs_year$year))
  testthat::expect_true(all(nethrs_year_loc$n == 4)) # 4 years of data
  testthat::expect_true(all(nethrs_year_yr$n == 10)) # 10 locations included
})

# test summarize_nethrs_stats wrapper & filters:
testthat::expect_warning(# no winter to adjust
  palo::summarize_nethrs_stats(
    df = sample_nethrs, location = 'RECR', datemin ='2010-01-01', by = 'season')) 

effort_RECR = palo::summarize_nethrs_stats(
  df = sample_nethrs, location = 'RECR', datemin ='2010-01-01', by = 'season',
  winter_adjust = FALSE)
effort_PGUP = palo::summarize_nethrs_stats(
  df = sample_nethrs, location = 'PGUP', datemax ='2010-10-15', by = 'month') 
testthat::test_that('site and mindate filters work appropriately', {
  testthat::expect_true(all(c(2010:2012) %in% effort_RECR$year))
  testthat::expect_false(2009 %in% effort_RECR$year)
  testthat::expect_true(all(effort_RECR$LOCATION == 'RECR'))
  
  testthat::expect_true(all(c(2009:2010) %in% effort_PGUP$year))
  testthat::expect_false(all(c(2011:2012) %in% effort_PGUP$year))
  testthat::expect_true(all(effort_PGUP$LOCATION == 'PGUP'))
})
# ggplot(effort_RECR, aes(year, nethours, color = season)) + geom_line()
# ggplot(effort_PGUP, aes(month, nethours, color = as.factor(year))) + 
#   geom_line() + geom_point()


# CAPTURES TESTING----------

# captures by season, without winter adjustment (by calendar year, so for
# example this joins Jan 2026 with December 2026)
captures_season_calendar = palo:::summarize_by(
  df = sample_band, stat = 'captures', by = 'season', winter_adjust = FALSE)
captures_season_calendar_loc = captures_season_calendar |> 
  group_by(LOC) |> count()
captures_season_calendar_seas = captures_season_calendar |> 
  group_by(season) |> count()

testthat::test_that(
  'captures by season works as expected', {
    
  testthat::expect_true(
    all(c('fall', 'spring', 'winter') %in% nethrs_season_calendar$season))
  testthat::expect_true('year' %in% names(nethrs_season_calendar))
  testthat::expect_false(2008 %in% captures_season_calendar$year)
  testthat::expect_true(length(unique(captures_season_calendar$SPEC)) == 123)
  testthat::expect_true(length(unique(captures_season_calendar$LOC)) == 18)
  testthat::expect_true(all(captures_season_calendar_loc$n == 1476)) 
  # >> 1476 = 3 seasons x 4 years x 123 spp
  testthat::expect_true(all(captures_season_calendar_seas$n == 8856)) 
  # >> 8856 = 4 years x 18 sites x 123 spp
})

# captures by season, with winter adjustment:
captures_season_adjust = palo:::summarize_by(
  df = sample_band, stat = 'captures', by = 'season', winter_adjust = TRUE)
captures_season_adjust_loc = captures_season_adjust |> group_by(LOC) |> count()
captures_season_adjust_seas = captures_season_adjust |> group_by(season) |>
  count()


testthat::test_that(
  'captures by season with winter adjustment works as expected', {
    
  testthat::expect_true(
    all(c('fall', 'spring', 'winter') %in% captures_season_adjust$season))
  testthat::expect_true('year_adjust' %in% names(captures_season_adjust))
  testthat::expect_true(
    all(c(2008:2012) %in% captures_season_adjust$year_adjust))
  #because Jan-Feb 2009 count as 2008 winter season
  testthat::expect_true(length(unique(captures_season_adjust$SPEC)) == 123)
  testthat::expect_true(length(unique(captures_season_adjust$LOC)) == 18)
  testthat::expect_true(all(captures_season_adjust_loc$n == 1845)) 
  # >> 1845 = 3 seasons x 5 years x 123 spp (because 2008 is filled in)
  testthat::expect_true(all(captures_season_adjust_seas$n == 11070)) 
  # >> 11070 = 5 years x 18 sites x 123 spp
  
})

capture_MUHO = palo::summarize_capture_stats(
  df = sample_band, location = 'MUHO', species = c('SWTH', 'AMGO'), 
  datemin ='2010-01-01', by = 'season', winter_adjust = FALSE)
capture_LACR = palo::summarize_capture_stats(
  df = sample_band, location = 'LACR', species = c('WIWA', 'SPTO'), 
  datemax ='2010-10-15', by = 'month') 
testthat::test_that(
  'site and mindate filters work appropriately in summarize_capture_stats ', {
    
  testthat::expect_true(all(c(2010:2012) %in% capture_MUHO$year))
  testthat::expect_false(2009 %in% capture_MUHO$year)
  testthat::expect_true(all(capture_MUHO$LOC == 'MUHO'))
  testthat::expect_true(all(capture_MUHO$SPEC %in% c('SWTH', 'AMGO')))
  
  testthat::expect_true(all(c(2009:2010) %in% capture_LACR$year))
  testthat::expect_false(all(c(2011:2012) %in% capture_LACR$year))
  testthat::expect_true(all(capture_LACR$LOC == 'LACR'))
  testthat::expect_true(all(capture_LACR$SPEC %in% c('WIWA', 'SPTO')))
})

# CAPTURE STATS-----------
capturedat = summarize_capture_stats(
   df = sample_band, location = c('MUHO', 'RECR'),
   species = c('SWTH', 'AMGO'), by = 'season', winter_adjust = FALSE)

effortdat = summarize_nethrs_stats(
   df = sample_nethrs, location = c('MUHO', 'RECR'),
   by = 'season', winter_adjust = FALSE)

capture_stats = calculate_capture_rates(
  captures = capturedat, effort = effortdat, by_species = TRUE,
  add_annual = TRUE)
