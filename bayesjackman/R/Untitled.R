lint_filter <- function(ifile, encoding = "unknown") {
  x = readLines(ifile, encoding = encoding, warn = FALSE)
  n = length(x)
  if (n == 0)
    return(x)
  p = knitr:::detect_pattern(x, tolower(knitr:::file_ext(ifile)))
  if (is.null(p))
    return(x)
  p = knitr::all_patterns[[p]]
  p1 = p$chunk.begin
  p2 = p$chunk.end
  i1 = grepl(p1, x)
  i2 = knitr:::filter_chunk_end(i1, grepl(p2, x))
  m = numeric(n)
  m[i1] = 1
  m[i2] = 2
  if (m[1] == 0)
    m[1] = 2
  for (i in seq_len(n - 1)) if (m[i + 1] == 0)
    m[i + 1] = m[i]
  out <- x
  out[m == 2 | i1] = ""
  # return inline code
  # x[m == 2] = stringr::str_replace_all(x[m == 2], p$inline.code,
  #                                      "")
  x
}
