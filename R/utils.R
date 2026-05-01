# General Utility Functions
# jmtools package

#' Source functions from an Rmd file
#'
#' Extracts and evaluates only function definitions from an R Markdown file
#' without executing the entire document. This is useful for reusing functions
#' defined in analysis notebooks without running all the analysis code.
#'
#' @param path Path to the .Rmd file
#' @param envir Environment in which to evaluate the functions (default: parent.frame())
#'
#' @return Invisibly returns the deparsed function definitions as a character vector
#'
#' @details
#' Extracts all R code chunks (delimited by ```\{r\} and ```), parses them with
#' R's parser, and evaluates only the top-level expressions that are function
#' assignments (`name <- function(...) {...}` or `name = function(...) {...}`).
#' Because parsing is done by R's parser rather than by tracking braces in raw
#' text, comments and strings containing `\{` or `\}` are handled correctly.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Source all functions from an analysis notebook
#' source_functions_from_rmd("analysis/01_preprocessing.Rmd")
#'
#' # Source into a specific environment
#' my_env <- new.env()
#' source_functions_from_rmd("analysis/utils.Rmd", envir = my_env)
#' }
source_functions_from_rmd <- function(path, envir = parent.frame()) {

  lines <- readLines(path, warn = FALSE)

  chunk_starts <- grep("^```\\{r", lines)
  chunk_ends   <- grep("^```$", lines)
  if (length(chunk_starts) != length(chunk_ends)) {
    stop("Mismatched code chunk delimiters in ", path)
  }

  # Concatenate all chunk bodies into one parseable text block. Skip empty
  # chunks (chunk_end == chunk_start + 1) so we don't accidentally pull the
  # chunk delimiters back as content via a backwards range.
  chunks <- Map(function(s, e) {
    if (e <= s + 1) character(0) else lines[(s + 1):(e - 1)]
  }, chunk_starts, chunk_ends)
  code <- paste(unlist(chunks), collapse = "\n")

  # Parse to a list of top-level expressions; R's parser handles strings and
  # comments correctly, unlike a regex-based brace counter.
  exprs <- parse(text = code)

  is_func_def <- function(e) {
    if (!is.call(e)) return(FALSE)
    op <- as.character(e[[1L]])
    if (!op %in% c("<-", "=", "<<-")) return(FALSE)
    rhs <- e[[3L]]
    is.call(rhs) && identical(as.character(rhs[[1L]]), "function")
  }

  func_exprs <- Filter(is_func_def, as.list(exprs))

  if (length(func_exprs) == 0) {
    warning("No function definitions found in ", path)
    return(invisible(character(0)))
  }

  for (e in func_exprs) eval(e, envir = envir)

  invisible(vapply(func_exprs,
                   function(e) paste(deparse(e), collapse = "\n"),
                   character(1)))
}


#' Custom minimal theme with axes
#'
#' A professional-looking ggplot2 theme based on `theme_minimal()` with black
#' axis lines and ticks added back. Optionally adds styling suitable for
#' faceted plots or publication-ready figures.
#'
#' @param base_size Base font size. Defaults to 11 for regular plots, or 6 when
#'   `publication = TRUE`.
#' @param base_family Base font family (default: "sans")
#' @param facet Logical. If `TRUE`, applies facet styling with panel borders
#'   and rotated x-axis labels (default: `FALSE`)
#' @param publication Logical. If `TRUE`, uses smaller font sizes suitable for
#'   publication figures (default: `FALSE`)
#'
#' @return A ggplot2 theme object
#'
#' @export
#'
#' @examples
#' library(ggplot2)
#'
#' # Regular plot with default font sizes
#' ggplot(mtcars, aes(wt, mpg)) +
#'   geom_point() +
#'   theme_jm()
#'
#' # Publication-ready plot with smaller fonts
#' ggplot(mtcars, aes(wt, mpg)) +
#'   geom_point() +
#'   theme_jm(publication = TRUE)
#'
#' # Faceted plot with panel borders
#' ggplot(mtcars, aes(wt, mpg)) +
#'   geom_point() +
#'   facet_wrap(~cyl) +
#'   theme_jm(facet = TRUE)
#'
#' @importFrom ggplot2 theme_minimal theme element_line element_blank element_rect element_text %+replace%
theme_jm <- function(base_size = NULL, base_family = "sans", facet = FALSE, publication = FALSE) {
  # Set default base_size based on publication parameter
  if (is.null(base_size)) {
    base_size <- if (publication) 6 else 11
  }

  base_theme <- ggplot2::theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    ggplot2::theme(
      axis.line = ggplot2::element_line(color = "black"),
      axis.ticks = ggplot2::element_line(color = "black")
    )

  if (facet) {
    base_theme <- base_theme %+replace%
      ggplot2::theme(
        # Add panel borders for facets
        panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
        # Remove axis lines to avoid double lines with panel border
        axis.line = ggplot2::element_blank(),
        # Rotate x-axis labels for better readability
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)
      )
  }

  base_theme
}


#' Get FIt-SNE executable path
#'
#' Detects and returns the correct path to the FIt-SNE (Fast Fourier Transform
#' accelerated t-SNE) executable based on the current system. Checks common
#' installation locations for HPC, macOS, and Ubuntu systems.
#'
#' @return Character string with the path to the fast_tsne executable
#'
#' @details
#' The function checks the following paths in order:
#' \enumerate{
#'   \item HPC: `/services/tools/fit-sne/1.2.1/bin/fast_tsne`
#'   \item macOS: `/usr/local/FIt-SNE/bin/fast_tsne`
#'   \item Ubuntu: `~/FIt-SNE/bin/fast_tsne`
#' }
#'
#' If none of these paths exist, the function stops with an error.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Get the path and use with Seurat
#' fast_tsne_path <- get_fast_tsne_path()
#' seurat_obj <- RunTSNE(seurat_obj,
#'                       tsne.method = "FIt-SNE",
#'                       fast_tsne_path = fast_tsne_path)
#' }
get_fast_tsne_path <- function() {
  # Define the possible paths
  path_hpc <- "/services/tools/fit-sne/1.2.1/bin/fast_tsne"
  path_osx <- "/usr/local/FIt-SNE/bin/fast_tsne"
  path_ubuntu <- path.expand("~/FIt-SNE/bin/fast_tsne")

  # Check which path exists and return accordingly

if (file.exists(path_hpc)) {
    return(path_hpc)
  } else if (file.exists(path_osx)) {
    return(path_osx)
  } else if (file.exists(path_ubuntu)) {
    return(path_ubuntu)
  } else {
    stop("fast_tsne executable not found in expected locations:\n",
         "  - HPC: ", path_hpc, "\n",
         "  - macOS: ", path_osx, "\n",
         "  - Ubuntu: ", path_ubuntu)
  }
}


#' Save a ggplot for publication
#'
#' A wrapper around [ggplot2::ggsave()] with sensible defaults for publication
#' figures. Provides convenient column width presets and aspect ratio control,
#' whilst still allowing full customisation when needed.
#'
#' @param filename File name to create on disk.
#' @param plot Plot to save, defaults to last plot displayed.
#' @param column Column width preset: `"single"` (90mm) or `"double"` (140mm).
#'   Ignored if `width` is specified.
#' @param aspect_ratio Aspect ratio as width/height (default: 4/3, giving 67.5mm
#'   height for single column). Ignored if both `width` and `height` are specified.
#' @param width Plot width. If `NULL`, calculated from `column` preset.
#' @param height Plot height. If `NULL`, calculated from `width` and `aspect_ratio`.
#' @param units Units for width and height (default: `"mm"`).
#' @param bg Background colour (default: `"transparent"`).
#' @param ... Additional arguments passed to [ggplot2::ggsave()], such as
#'   `device`, `dpi`, `scale`, etc.
#'
#' @return Invisibly returns the `filename`.
#'
#' @details
#' Default dimensions follow common journal requirements:
#' \itemize{
#'   \item Single column: 90 x 67.5 mm (4:3 aspect ratio)
#'   \item Double column: 140 x 105 mm (4:3 aspect ratio)
#' }
#'
#' The function calculates dimensions as follows:
#' \enumerate{
#'   \item If both `width` and `height` are specified, use them directly
#'   \item If only `width` is specified, calculate `height` from `aspect_ratio`
#'   \item If only `height` is specified, calculate `width` from `aspect_ratio`
#'   \item If neither is specified, use `column` preset and `aspect_ratio`
#' }
#'
#' @export
#'
#' @examples
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(wt, mpg)) + geom_point()
#'
#' \dontrun{
#' # Save with default single-column dimensions (90 x 67.5 mm)
#' ggsave_for_publication("figure_1A.pdf", p)
#'
#' # Save as double-column width
#' ggsave_for_publication("figure_2A.pdf", p, column = "double")
#'
#' # Custom aspect ratio (16:9)
#' ggsave_for_publication("figure_3A.pdf", p, aspect_ratio = 16/9)
#'
#' # Override with custom dimensions
#' ggsave_for_publication("figure_4A.pdf", p, width = 110, height = 80)
#'
#' # Specify width only, height calculated from aspect ratio
#' ggsave_for_publication("figure_5A.pdf", p, width = 120)
#' }
#'
#' @seealso [ggplot2::ggsave()]
ggsave_for_publication <- function(filename,
                                   plot = ggplot2::last_plot(),
                                   column = c("single", "double"),
                                   aspect_ratio = 4 / 3,
                                   width = NULL,
                                   height = NULL,
                                   units = "mm",
                                   bg = "transparent",
                                   ...) {


  column <- match.arg(column)

  # Define column widths in mm
  column_widths <- c(single = 90, double = 140)

  # Calculate dimensions based on what's provided
  if (is.null(width) && is.null(height)) {
    # Neither specified: use column preset and aspect ratio
    width <- column_widths[column]
    height <- width / aspect_ratio
  } else if (!is.null(width) && is.null(height)) {
    # Width specified: calculate height from aspect ratio
    height <- width / aspect_ratio
  } else if (is.null(width) && !is.null(height)) {
    # Height specified: calculate width from aspect ratio
    width <- height * aspect_ratio
  }
  # If both width and height are specified, use them directly

  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    units = units,
    bg = bg,
    ...
  )
}
