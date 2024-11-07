# The  knit_to_gitbook.R script uses the bookdown  package to render the book in GitBook format.

# To run the  knit_to_gitbook.R  script, click on the Source  button at the top of the window.

# load the bookdown package
library(bookdown)

# removes all files from the _book  directory
clean_book(clean = TRUE)

# renders the book in GitBook format
bookdown::render_book("index.Rmd", "bookdown::gitbook")

# opens the GitBook in your browser
# rstudioapi::viewer("_book/index.html")

# To see the book in your Viewer pane, open 'index.Rmd' in your main directory and click the 'Knit' button.

