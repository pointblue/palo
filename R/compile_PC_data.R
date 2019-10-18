compile_PC_data <- function(dir, pattern = '.csv') {
  # get names of all files in directory that match pattern
  file_list <- list.files(here::here(dir), pattern = pattern,
                          full.names = TRUE)

  # compile master data set (assuming CADC Project Leader format)
  dat <- map_dfr(file_list,
                 ~read_csv(.x,
                           col_types = 'cccccnDtttncnccccccclcccc')) %>%
    filter_all(any_vars(!is.na(.))) %>%
    mutate_at(vars(Project:Protocol, `Time Bin`, Spp:`Breeding Status`,
                   Researcher:`Data Status`),
              as.factor)

  test_that('Missing values in critical fields', {
    expect_false(any(is.na(dat$Project)))
    expect_false(any(is.na(dat$`Study Area`)))
    expect_false(any(is.na(dat$Transect)))
    expect_false(any(is.na(dat$Point)))
    expect_false(any(is.na(dat$Protocol)))
    expect_false(any(is.na(dat$Spp)))
    expect_false(any(is.na(dat$`Distance Bin ID`)))
    expect_false(any(is.na(dat$Count)))
  })

  # convert distance bins to corresponding numeric values
  add_distance_key(dat)

}

add_distance_key <- function(dat) {
  key <- dat %>%
    select(Protocol, `Distance Bin ID`, `Distance Bin`) %>%
    distinct() %>%
    mutate(mindist = case_when(grepl('>', `Distance Bin`) ~
                                 gsub('>', '', `Distance Bin`),
                               grepl('<', `Distance Bin`) ~ '0',
                               grepl('\\sto\\s', `Distance Bin`) ~
                                 gsub(' to .*$', '', `Distance Bin`),
                               `Distance Bin` == 'FlyOver' ~ NA_character_),
           maxdist = case_when(grepl('>', `Distance Bin`) ~ NA_character_,
                               grepl('<', `Distance Bin`) ~
                                 gsub('<', '', `Distance Bin`),
                               grepl('\\sto\\s', `Distance Bin`) ~
                                 gsub('^.* to ', '', `Distance Bin`),
                               `Distance Bin` == 'FlyOver' ~ NA_character_)) %>%
    mutate_at(vars(mindist:maxdist), as.numeric)
  left_join(dat, key, by = c('Protocol', 'Distance Bin ID', 'Distance Bin'))
}
