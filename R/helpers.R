.make_folds <- function(n, K = 2, seed = NULL) {
  assertthat::assert_that(is.numeric(n) && length(n) == 1 && n >= 2, msg = "n must be >= 2")
  assertthat::assert_that(is.numeric(K) && length(K) == 1 && K >= 2, msg = "K must be >= 2")
  if (!is.null(seed)) set.seed(seed)
  ids <- sample.int(n)
  folds <- split(ids, rep(seq_len(K), length.out = n))
  folds
}



.is_scalar_character <- function(x) is.character(x) && length(x) == 1 && !is.na(x)
.is_scalar_numeric   <- function(x) is.numeric(x) && length(x) == 1 && !is.na(x)


.assert_list_has_names <- function(x, nms, what = "object") {
  assert_that(is.list(x), msg = sprintf("%s must be a list", what))
  missing <- setdiff(nms, names(x))
  if (length(missing) > 0) {
    stop(sprintf("%s is missing required names: %s", what, paste(missing, collapse = ", ")), call. = FALSE)
  }
  invisible(TRUE)
}
