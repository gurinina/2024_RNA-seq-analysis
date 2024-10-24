# The  knit_to_gitbook.R  script is a simple R script that uses the  bookdown  package to render the book in GitBook format. The script first cleans the book by removing all files from the  _book  directory. It then renders the book in GitBook format and opens the resulting GitBook in the Viewer pane. Finally, it starts a local web server to serve the GitBook.
# To run the  knit_to_gitbook.R  script, click on the  Source  button in the script editor. This will render the book in GitBook format and open the resulting GitBook in the Viewer pane.
#
# The GitBook will be served at  http://
#

# load the bookdown package
library(bookdown)

# Clean the book; this will remove all files from the _book directory
clean_book(clean = TRUE)

# Render the book in GitBook format
bookdown::render_book("index.Rmd", "bookdown::gitbook")

# Open the resulting GitBook in the Viewer pane
viewer_path <- normalizePath("_book/index.html")
rstudioapi::viewer(viewer_path)

#
# servr::httd("_book")
# servr::daemon_stop(1)
#
