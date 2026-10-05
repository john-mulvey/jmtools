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
  expect_equal(env$add_one(1), 2)
  expect_equal(env$double(3), 6)
  expect_equal(env$quoted_name(), "quoted")
})

test_that("top-level pkg::fn() and obj$fn() calls do not cause an error", {
  path <- local_rmd(c(
    "```{r}",
    "add_one <- function(x) x + 1",
    "stats::median(1:3)",
    "base:::identity(1)",
    "obj$method(1)",
    "fns[[1]](1)",
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

test_that("braces inside strings and comments do not break parsing", {
  path <- local_rmd(c(
    "```{r}",
    "open_brace <- function() {",
    "  # a stray { in a comment",
    "  \"}\"",
    "}",
    "after <- function() 'after'",
    "```"
  ))
  env <- new.env()

  source_functions_from_rmd(path, envir = env)

  expect_equal(env$open_brace(), "}")
  expect_equal(env$after(), "after")
})

test_that("the deparsed definitions are returned invisibly", {
  path <- local_rmd(c(
    "```{r}",
    "add_one <- function(x) x + 1",
    "```"
  ))

  expect_invisible(
    definitions <- source_functions_from_rmd(path, envir = new.env())
  )
  expect_identical(definitions, "add_one <- function(x) x + 1")
})

test_that("a file with no function definitions gives a warning", {
  path <- local_rmd(c(
    "```{r}",
    "x <- 1",
    "```"
  ))

  expect_warning(
    definitions <- source_functions_from_rmd(path, envir = new.env()),
    "No function definitions found"
  )
  expect_identical(definitions, character(0))
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
  # The bash and python chunks are not valid R, so the call would fail to
  # parse if their contents were included.
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
