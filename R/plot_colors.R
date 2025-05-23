#' Plot Colorpalettes
#'
#' Function `plot_colors` plots colors and their labels for easy
#' visualization of a colorpalette.
#'
#' @export
#'
#' @param x A vector of colors.
#' @param labels A vector of labels for colors. If omitted, given color names are used.
#'
#' @seealso See e.g. the [colorpalette()] data and `RColorBrewer`
#' package for ready-made color palettes.
#'
#' @examples
#' plot_colors(colorpalette[[5]], labels = c("one", "two", "three", "four", "five"))
#'
#' plot_colors(colorpalette[[10]])
#'
#' plot_colors(1:7)
#'
#' plot_colors(c("yellow", "orange", "red", "purple", "blue", "green"))
#'
#' plot_colors(grDevices::rainbow(15))
plot_colors <- function(x, labels = NULL) {
  stopifnot_(
    all(isColor(x)),
    "Please provide a vector of colors."
  )
  if (is.null(labels)) {
    labels <- rev(x)
  } else if (length(x) != length(labels)) {
    warning_("The length of {.arg labels} does not match the length of {.arg x}. 
             Labels were not used.")
    labels <- rev(x)
  } else {
    labels <- rev(labels)
  }
  par(mai = c(0.1, max(graphics::strwidth(x, "inch") + 0.4, na.rm = TRUE), 0.1, 0.4))
  graphics::barplot(rep(1, length(x)),
          col = rev(x), space = 0.2, axes = FALSE,
          names.arg = labels, cex.names = 0.8, horiz = TRUE, las = 1
  )
}
