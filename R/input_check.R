# ...............................
# fft metrics helper functions
# ...............................

#' Check \code{for.year} argument
#' 
#' Ensure that the for.year argument is a numeric scalar value
#' @param for.year input argument to be checked
#' @return Stops execution in case of failed check
assert.for.year = function(for.year) {
  # for.year, if not NULL, must be a numeric scalar value
  assert_numeric(for.year, null.ok = TRUE)
  assert_scalar(for.year, null.ok = TRUE)
}

#' Check for numeric vector
#' 
#' Ensure that the given argument is a numeric vector
#' @param arg input argument to be checked
#' @return Stops execution in case of failed check
assert.numeric.vector = function(arg) {
  # must be a numeric vector
  assert_numeric(arg)
  assert_vector(arg)
}

#' Check for equal length inputs
#' 
#' Ensure that the given arguments are equal length vectors
#' @param arg0 At least one input argument to be checked needs to be given
#' @param ... Other inputs whose length equality needs to be checked
#' @return Stops execution in case of failed check
assert.equal.length = function(arg0, ...) {
  args = list(...)
  l = length(arg0)
  for (i in args) {
    if (length(i) != l) {
      print(l)
      print(length(i))
      stop("Inputs must be of same length!")
    }
  }
}

