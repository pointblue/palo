#' @rdname summarize_capture_stats
#' @export
summarize_nethrs_stats = function(path = NULL, df = NULL, location = NULL, 
                                  datemin = NULL, datemax = NULL,
                                  by = 'year', winter_adjust = TRUE) {
  
  if (is.null(df) & !is.null(path)) {
    df = foreign::read.dbf(path) |> 
      dplyr::select(.data$PROJECT:.data$DUPE) # (drop extra columns labeled "X")
  } else if (is.null(df) & is.null(path)) {
    stop('no input dataframe or path to a .dbf provided')
  }
  
  # optional filtering
  if (!is.null(location)) {
    df = df |> dplyr::filter(.data$LOCATION %in% location)
  }
  if (!is.null(datemin)) {
    df = df |> dplyr::filter(.data$DATE >= datemin)
  }
  if (!is.null(datemax)) {
    df = df |> dplyr::filter(.data$DATE <= datemax)
  }
  
  summarize_by(df, by = by, stat = 'nethours', winter_adjust = winter_adjust)
  
}