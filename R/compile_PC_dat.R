#' Compile point count survey data from a directory
#'
#' @param dir Directory containing point count survey data downloaded from CADC
#' @param pattern Regex pattern for identifying files to include (see \code{\link[base]{list.files}})
#'
#' @return Dataframe of compiled data, excluding any rows with all NA values.
#' @details Also tests for missing values in critical fields.
#' @export
#' @import tidyr dplyr
#' @importFrom purrr map_dfr
#' @importFrom testthat test_that expect_false
#' @importFrom rlang .data
#' @importFrom here here

compile_PC_dat <- function(dir, pattern = '.csv') {
  # get names of all files in directory that match pattern
  file_list <- list.files(here::here(dir), pattern = pattern,
                          full.names = TRUE)

  # compile master data set (assuming CADC Project Leader format)
  dat <- map_dfr(file_list,
                 ~read_csv(.x,
                           col_types = 'cccccnDtttncnccccccclcccc'))
  dat <- dat %>%
    filter_all(any_vars(!is.na(.data))) %>%
    mutate_at(vars(.data$Project:.data$Protocol,
                   .data$`Time Bin`,
                   .data$Spp:.data$`Breeding Status`,
                   .data$Researcher:.data$`Data Status`),
              as.factor)

  testthat::test_that('Missing values in critical fields', {
    testthat::expect_false(any(is.na(dat$Project)))
    testthat::expect_false(any(is.na(dat$`Study Area`)))
    testthat::expect_false(any(is.na(dat$Transect)))
    testthat::expect_false(any(is.na(dat$Point)))
    testthat::expect_false(any(is.na(dat$Protocol)))
    testthat::expect_false(any(is.na(dat$Spp)))
    testthat::expect_false(any(is.na(dat$`Distance Bin ID`)))
    testthat::expect_false(any(is.na(dat$Count)))
  })
}
