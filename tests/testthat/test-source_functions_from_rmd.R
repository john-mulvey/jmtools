test_that("function definitions are sourced and other code is not evaluated", {
  path <- local_rmd(c(
    "```{r}",
    "add_one <- function(x) x + 1",
    "side_effect <- 1",
    "stop('top-level code must not be evaluated')",
    "```",
    "",
    "```{r}",
    "double = function(x) x * 2",
    "\"quoted_name\" <- function() 'quoted'",
    "```"
  ))
  env <- new.env()

  source_functions_from_rmd(path, envir = env)

  expect_setequal(ls(env), c("add_one", "double", "quoted_name"))
})

test_that("top-level pkg::fn() and obj$fn() calls do not cause an error", {
  # Two arguments make each call the same length as an assignment, so the
  # call head itself must be checked rather than the call's length.
  path <- local_rmd(c(
    "```{r}",
    "add_one <- function(x) x + 1",
    "stats::median(1:3, na.rm = TRUE)",
    "base:::paste(1, 2)",
    "obj$method(1, 2)",
    "fns[[1]](1, 2)",
    "```"
  ))
  env <- new.env()

  source_functions_from_rmd(path, envir = env)

  expect_identical(ls(env), "add_one")
})

test_that("functions assigned into an object are skipped", {
  path <- local_rmd(c(
    "```{r}",
    "add_one <- function(x) x + 1",
    "obj$method <- function() 1",
    "obj[['method']] <- function() 1",
    "```"
  ))
  env <- new.env()

  source_functions_from_rmd(path, envir = env)

  expect_identical(ls(env), "add_one")
})

test_that("chunk headers with labels and options are recognised", {
  path <- local_rmd(c(
    "```{r setup, include = FALSE}",
    "from_labelled <- function() 1",
    "```",
    "",
    "```{r, echo = FALSE}",
    "from_options <- function() 2",
    "```"
  ))
  env <- new.env()

  source_functions_from_rmd(path, envir = env)

  expect_setequal(ls(env), c("from_labelled", "from_options"))
})

test_that("empty chunks are ignored", {
  path <- local_rmd(c(
    "```{r}",
    "```",
    "",
    "```{r}",
    "add_one <- function(x) x + 1",
    "```"
  ))
  env <- new.env()

  source_functions_from_rmd(path, envir = env)

  expect_identical(ls(env), "add_one")
})

test_that("chunks for other engines are skipped", {
  # The bash, python and rcpp chunks are not valid R, so the call would fail
  # to parse if their contents were included.
  path <- local_rmd(c(
    "```{r}",
    "from_r <- function() 1",
    "```",
    "",
    "```{bash}",
    "ls -l $HOME",
    "```",
    "",
    "```{python}",
    "def from_python():",
    "    return 1",
    "```",
    "",
    "```{rcpp}",
    "int from_rcpp() { return 1; }",
    "```"
  ))
  env <- new.env()

  source_functions_from_rmd(path, envir = env)

  expect_identical(ls(env), "from_r")
})

test_that("plain markdown code blocks in the prose are skipped", {
  path <- local_rmd(c(
    "Run this in a terminal:",
    "",
    "```",
    "make all",
    "```",
    "",
    "```{r}",
    "from_r <- function() 1",
    "```"
  ))
  env <- new.env()

  source_functions_from_rmd(path, envir = env)

  expect_identical(ls(env), "from_r")
})

test_that("an unclosed chunk is an error", {
  path <- local_rmd(c(
    "```{r}",
    "add_one <- function(x) x + 1",
    "```",
    "",
    "```{r}",
    "never_closed <- function() 1"
  ))

  expect_error(
    source_functions_from_rmd(path, envir = new.env()),
    "Unclosed code chunk starting at line 5"
  )
})

test_that("without `functions`, a name defined twice takes its last definition", {
  path <- local_rmd(c(
    "```{r}",
    "redefined <- function() 'first'",
    "redefined <- function() 'second'",
    "```"
  ))
  env <- new.env()

  source_functions_from_rmd(path, envir = env)

  expect_equal(env$redefined(), "second")
})

test_that("`functions` sources only the named functions", {
  path <- local_rmd(c(
    "```{r}",
    "wanted <- function() 'wanted'",
    "also_wanted <- function() 'also wanted'",
    "unwanted <- function() 'unwanted'",
    "```"
  ))
  env <- new.env()

  definitions <- source_functions_from_rmd(path, c("wanted", "also_wanted"), envir = env)

  expect_setequal(ls(env), c("wanted", "also_wanted"))
  expect_identical(definitions, c("wanted <- function() \"wanted\"",
                                  "also_wanted <- function() \"also wanted\""))
})

test_that("requesting a function that is not defined is an error", {
  path <- local_rmd(c(
    "```{r}",
    "defined <- function() 1",
    "```"
  ))
  env <- new.env()

  expect_error(
    source_functions_from_rmd(path, c("defined", "absent", "also_absent"), envir = env),
    "Not defined in .*: absent, also_absent"
  )
  expect_identical(ls(env), character(0))
})

test_that("requesting a function from a file with no definitions is an error", {
  path <- local_rmd(c(
    "```{r}",
    "x <- 1",
    "```"
  ))

  expect_error(
    source_functions_from_rmd(path, "absent", envir = new.env()),
    "Not defined in .*: absent"
  )
})

test_that("requesting a function that is defined more than once is an error", {
  path <- local_rmd(c(
    "```{r}",
    "redefined <- function() 'first'",
    "redefined <- function() 'second'",
    "```"
  ))

  expect_error(
    source_functions_from_rmd(path, "redefined", envir = new.env()),
    "Defined more than once in .*: redefined"
  )
})

test_that("an unrequested function defined more than once is not an error", {
  path <- local_rmd(c(
    "```{r}",
    "wanted <- function() 'wanted'",
    "redefined <- function() 'first'",
    "redefined <- function() 'second'",
    "```"
  ))
  env <- new.env()

  source_functions_from_rmd(path, "wanted", envir = env)

  expect_identical(ls(env), "wanted")
})

test_that("a requested function calling an unrequested helper from the file is an error", {
  path <- local_rmd(c(
    "```{r}",
    "plot_panels <- function(x) x",
    "fit_model <- function(x) plot_panels(x)",
    "```"
  ))

  expect_error(
    source_functions_from_rmd(path, "fit_model", envir = new.env()),
    "fit_model() calls plot_panels()",
    fixed = TRUE
  )
})

test_that("a requested function works when its helpers are also requested", {
  path <- local_rmd(c(
    "```{r}",
    "plot_panels <- function(x) x * 2",
    "fit_model <- function(x) plot_panels(x) + 1",
    "```"
  ))
  env <- new.env()

  source_functions_from_rmd(path, c("fit_model", "plot_panels"), envir = env)

  expect_equal(env$fit_model(1), 3)
})

test_that("calls to functions not defined in the file are not treated as helpers", {
  path <- local_rmd(c(
    "```{r}",
    "summarise_values <- function(x) stats::median(sum(x)) + defined_elsewhere(x)",
    "```"
  ))
  env <- new.env()

  source_functions_from_rmd(path, "summarise_values", envir = env)

  expect_identical(ls(env), "summarise_values")
})
