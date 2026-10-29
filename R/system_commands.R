# Package-local bindings provide a stable seam for testthat scheduler mocks.
# Forward to the public base APIs instead of copying their function bodies:
# aliases expose base R's .Internal calls to R CMD check as package code.
system <- function(...) base::system(...)
system2 <- function(...) base::system2(...)
