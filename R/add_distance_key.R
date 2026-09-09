#' Convert Distance Bins to numeric mindist and maxdist fields
#'
#' @param dat Dataframe containing the field 'Distance Bin'
#'
#' @return tibble
#' @export

add_distance_key <- function(dat) {
  key <- dat |> 
    dplyr::select('Protocol', 'Distance Bin ID', 'Distance Bin') |> 
    dplyr::distinct() |>
    dplyr::mutate(
      mindist = dplyr::case_when(
        grepl('>', 'Distance Bin') ~ gsub('>', '', 'Distance Bin'),
        grepl('<', 'Distance Bin') ~ '0',
        grepl('\\sto\\s', 'Distance Bin') ~ gsub(' to .*$', '', 'Distance Bin'),
        'Distance Bin' == 'FlyOver' ~ NA_character_),
      maxdist = dplyr::case_when(
        grepl('>', 'Distance Bin') ~ NA_character_, 
        grepl('<', 'Distance Bin') ~ gsub('<', '', 'Distance Bin'),
        grepl('\\sto\\s', 'Distance Bin') ~ gsub('^.* to ', '', 'Distance Bin'),
        'Distance Bin' == 'FlyOver' ~ NA_character_)) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::any_of(c('mindist', 'maxdist')), as.numeric))
  
  dplyr::left_join(dat, key, 
                  by = c('Protocol', 'Distance Bin ID', 'Distance Bin'))
}
