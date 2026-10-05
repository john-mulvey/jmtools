# Write `lines` to a temporary .Rmd file that is deleted when the calling
# test finishes, and return its path.
local_rmd <- function(lines, env = parent.frame()) {
  path <- withr::local_tempfile(fileext = ".Rmd", .local_envir = env)
  writeLines(lines, path)
  path
}
