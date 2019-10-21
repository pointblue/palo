#' Convert Distance Bins to numeric mindist and maxdist fields
#'
#' @param dat Dataframe containing the field 'Distance Bin'
#'
#' @return Dataframe
#' @export
#' @importFrom rlang .data

add_distance_key <- function(dat) {
  key <- dat %>%
    select(.data$Protocol, .data$`Distance Bin ID`, .data$`Distance Bin`) %>%
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
    mutate_at(vars(.data$mindist:.data$maxdist), as.numeric)
  left_join(dat, key, by = c('Protocol', 'Distance Bin ID', 'Distance Bin'))
}
