# Package-local bindings provide a stable seam for testthat scheduler mocks.
# Keep the production behavior identical to base R's command runners.
system <- base::system
system2 <- base::system2
