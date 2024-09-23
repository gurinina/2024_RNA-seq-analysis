# Global options for code visibility, formatting, and output width
options(width = 80)

knitr::opts_chunk$set(
  echo = TRUE,                    # Show code by default
  tidy.opts = list(width.cutoff = 80),  # Format code width
  tidy = TRUE
)
