#' Summarize banding statistics (internal function)
#'
#' Internal functions called by [summarize_capture_stats()] and
#' [summarize_nethrs_stats()] and not intended to be called independently.
#'
#' @param df Input dataframe
#' @param by how statistics should be summarized, either "month", "year", or
#'   "season"
#' @param stat the type of stats to be summarized; either "nethours" or
#'   "captures"
#' @param winter_adjust logical
#'
#' @details See documentation for [summarize_capture_stats()] for details. Note
#' that the only valid options for the argument `stat` are currently: "nethours"
#' or "captures". This value determines whether the NETHOURS field is summed or
#' the number of BANDNUMBs are counted, and the name of the column containing
#' the results.
#'
#' @returns tibble
#' @export
#' @keywords internal
#' @importFrom rlang :=

summarize_by = function(df, by = 'month', stat = c('nethours', 'captures'),
                        winter_adjust = NULL) {
  
  if (by == 'month') {
    
    if (!'month' %in% names(df)) {
      if ('DATE' %in% names(df)) {
        df = df |> 
          dplyr::mutate(month = format(.data$DATE, '%m') |> as.numeric(),
                        year = format(.data$DATE, '%Y') |> as.numeric())
      } else {
        stop('df must contain either "month" or "DATE" fields')
      }
    }
    
    # define relevant columns to group by (any of these that may be present)
    cols = c("LOC", "LOCATION", "SPEC", "year", "month")
    
  } else if (by == 'year') {
    
    if (!'year' %in% names(df)) {
      
      if ('DATE' %in% names(df)) {
        df = df |> 
          dplyr::mutate(year = format(.data$DATE, '%Y') |> as.numeric())
      } else {
        stop('df must contain either "year" or "DATE" fields')
      }
    }
    
    # define relevant columns to group by
    cols = c("LOC", "LOCATION", "SPEC", "year")
    
  } else if (by == 'season') {
    
    if (!'season' %in% names(df)) {
      
      # apply default season assumptions based on groups of months
      if ('DATE' %in% names(df)) {
        df = df |> 
          dplyr::mutate(
            month = format(.data$DATE, '%m') |> as.numeric(),
            year = format(.data$DATE, '%Y') |> as.numeric(),
            season = dplyr::case_when(
              .data$month >= 3 & .data$month <= 7 ~ 'spring', # 5 months
              .data$month >= 8 & .data$month <= 10 ~ 'fall', # 3 months
              TRUE ~ 'winter'))
        
      } else {
        stop('df must contain either "season" or "DATE" fields')
      }
    }
    
    if (winter_adjust) {
      
      if (!'year' %in% names(df)) {
        
        if ('DATE' %in% names(df)) {
          df = df |> 
            dplyr::mutate(
              year = format(.data$DATE, '%Y'))
          
        } else {
          
          stop('df must contain either "year" or "DATE" fields to adjust year assignments for winter seasons' )
          
        }
      }
      
      if ('winter' %in% unique(df$season)) {
        df = df |> 
          dplyr::mutate(
            year_adjust = dplyr::case_when(
              .data$season == 'winter' & .data$month < 6 ~ .data$year - 1,
              TRUE ~ .data$year)
          )
        
        # define relevant columns to group by (any of these that may be present)
        cols = c("LOC", "LOCATION", "SPEC", "year_adjust", "season")
        
      } else {
        warning('no "winter" seasons found to adjust')
        
        cols = c("LOC", "LOCATION", "SPEC", "year", "season")
      }
      
      
      
    } else {
      # group by original "year" in df or by calendar year determined from DATE
      
      cols = c("LOC", "LOCATION", "SPEC", "year", "season")

    }
    
  }
  
  # summarize by the columns defined above, depending on selected "stat"
  df = df |> 
    dplyr::group_by(
      dplyr::across(
        dplyr::any_of(cols)))
  
  if (stat == 'nethours') {
    if('NETHOURS' %in% names(df)) {
      df = df |> 
        dplyr::summarize(
          y = sum(.data$NETHOURS),
          .groups = 'drop'
        )
    } else {
      stop('the "NETHOURS" field is missing from the input data')
    }
  } else if (stat == 'captures') {
    if ('BANDNUMB' %in% names(df)) {
      df = df |> 
        dplyr::summarize(
          y = length(.data$BANDNUMB),
          .groups = 'drop'
        )
    } else {
    stop('the "BANDNUMB" field is missing from the input data')
    }
  }

  # expand to represent any missing values (0 captures or 0 nethours during a
  # given time frame)
  df = fill_missing(df) 

  # rename results field according to "stat"
  df |> dplyr::rename(!!stat := .data$y)
  
}
