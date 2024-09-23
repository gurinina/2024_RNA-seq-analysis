solutions <- function(solution_text) {
  output_format <- knitr::opts_knit$get("rmarkdown.pandoc.to")

  if (output_format == "html") {
    # For GitBook (HTML), output as raw HTML
    knitr::asis_output(paste0(
      '<details>\n',
      '<summary><i>Click here to see the solution</i></summary>\n',
      '<p>', solution_text, '</p>\n',
      '</details>\n'
    ))
  } else {
    # For non-HTML formats (EPUB, PDF), hide the solution
    cat("")  # No output for non-HTML formats
  }
}
