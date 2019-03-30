# provide default formats if necessary
formats <- commandArgs(trailingOnly = TRUE)
if (length(formats) == 0) {
  formats <- "bookdown::pdfbook"
}
# render the book to all formats
for (fmt in formats)
  bookdown::render_book("index.Rmd", fmt, quiet = FALSE)
