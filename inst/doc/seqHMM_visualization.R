## ----settingdata, message=FALSE, cache=FALSE, echo = FALSE, eval = TRUE----
library("seqHMM")

data("biofam", package = "TraMineR")
biofam_seq <- seqdef(biofam[, 10:25], start = 15, labels = c(
  "parent", "left", "married", "left+marr", "child", "left+child",
  "left+marr+ch", "divorced")
)

data("biofam3c")
marr_seq <- seqdef(biofam3c$married, start = 15, alphabet = c(
  "single", "married", "divorced"),
  cpal = c("violetred2", "darkgoldenrod2", "darkmagenta")
)
child_seq <- seqdef(biofam3c$children,
  start = 15,
  alphabet = c("childless", "children"),
  cpal = c("darkseagreen1", "coral3")
)
left_seq <- seqdef(biofam3c$left, start = 15, alphabet = c(
  "with parents","left home"),
  cpal = c("lightblue", "red3")
)

## ----seqplotSingle, fig.width=6.5, fig.height=3, dev.args=list(pointsize=10), echo=FALSE, fig.cap="Sequence index plot (left) and state distribution plot (right) of annual family states for 100 individuals from the \\texttt{biofam} data.", fig.align='center', fig.keep='last', message=FALSE----
library(patchwork)
p1 <- stacked_sequence_plot(biofam_seq[1:100, ], type = "i", legend_position = "none")
p2 <- stacked_sequence_plot(biofam_seq[1:100, ], type = "d", legend_position = "right")
p1 + p2 & ggplot2::xlab("Age")

## ----seqplotMulti, fig.width=6.5, fig.height=3.7, dev.args=list(pointsize=10), echo=FALSE, fig.cap="Stacked sequence plot of the first ten individuals in the \\texttt{biofam} data. The top plot shows the original sequences, and the three bottom plots show the sequences in the separate channels for the same individuals. The sequences are in the same order in each plot, i.e., the same row always matches the same individual.", fig.align='center', fig.keep='last', message=FALSE----
seq_data <- list(
  Original = biofam_seq[1:10, ], 
  Marriage = marr_seq[1:10, ], 
  Parenthood = child_seq[1:10, ],
  Residence = left_seq[1:10, ]
)
stacked_sequence_plot(
  seq_data, sort_by = "start", sort_channel = "Original"
  ) & ggplot2::xlab("Age")

## ----HMMplot, fig.width=6.5, fig.height=4.5, dev.args=list(pointsize=10), echo=FALSE, fig.cap="Illustrating a left-to-right hidden Markov model for the multichannel \\texttt{biofam} data as a directed graph. Pies represent the hidden states, with emission probabilities of combined observations as slices. Arrows illustrate transition probabilities between the hidden states. Probabilities of starting in each state are shown next to the pies.", fig.align='center', fig.keep='last', cache = TRUE, message=FALSE----
plot(hmm_biofam,
  layout = matrix(c(
    1, 2, 3, 4, 2,
    1, 1, 1, 1, 0
  ), ncol = 2),
  # varying curvature of edges
  edge.curved = c(0, -0.8, 0.6, 0, 0, -0.8, 0),
  # thinner edges and arrows
  cex.edge.width = 0.8, edge.arrow.size = 1,
  # fixing axes to the right scale
  xlim = c(0.5, 4.5), ylim = c(-0.5, 1.5), rescale = FALSE,
  # different legend properties
  with.legend = "bottom", legend.prop = 0.3, ncol.legend = 2,
  # distance of vertex labels to vertices
  vertex.label.dist = 1.1,
  # threshold for emission probabilities not shown as separate slices
  combine.slices = 0.02, combined.slice.label = "others (emission prob. < 0.02)"
)

## ----settingsequences, message = FALSE------------------------------------
library("seqHMM")

data("biofam3c")
marr_seq <- seqdef(
  biofam3c$married,
  start = 15,
  alphabet = c("single", "married", "divorced"),
  cpal = c("violetred2", "darkgoldenrod2", "darkmagenta")
)
child_seq <- seqdef(
  biofam3c$children,
  start = 15,
  alphabet = c("childless", "children"),
  cpal = c("darkseagreen1", "coral3")
)
left_seq <- seqdef(
  biofam3c$left,
  start = 15,
  alphabet = c("with parents", "left home"),
  cpal = c("lightblue", "red3")
)

## ----plottingsequences, fig.width=5, fig.height=3, dev.args=list(pointsize=10), fig.cap="Stacked sequence plot of annual state distributions in the three-channel \\texttt{biofam} data. This is the default output of the \\texttt{ssplot} function.", fig.keep='last', fig.align='center', cache=FALSE, echo = TRUE----
p <- stacked_sequence_plot(
  x = list(Marriage = marr_seq, Parenthood = child_seq, Residence = left_seq)
)
p

## ----plottingsequencesHMM, fig.width=5, fig.height=5, dev.args=list(pointsize=10), fig.cap="Stacked sequence plot of observations and hidden state paths using a hidden Markov model object.", fig.keep='last', fig.align='center', cache=FALSE, echo = TRUE----
data("hmm_biofam")
stacked_sequence_plot(x = hmm_biofam, plots = "both")

## ----seqIplot, fig.width=5, fig.height=3, dev.args=list(pointsize=10), fig.cap="Sequence index plot showing observed sequences sorted by the third channel, residence.", fig.align='center', cache=TRUE, echo = TRUE----
stacked_sequence_plot(
  hmm_biofam,
  type = "index", sort_by = "start", sort_channel = 3
)

## ----code_plottingHMMbasic, fig.height=5, fig.width=8, echo=TRUE, fig.align='center', fig.keep='last', cache = TRUE, eval = TRUE, fig.cap="A default plot of a hidden Markov model."----
plot(hmm_biofam)

## ----HMMplotCode, echo=TRUE, eval=FALSE-----------------------------------
# plot(hmm_biofam,
#   layout = matrix(c(
#     1, 2, 3, 4, 2,
#     1, 1, 1, 1, 0
#   ), ncol = 2),
#   xlim = c(0.5, 4.5), ylim = c(-0.5, 1.5), rescale = FALSE,
#   edge.curved = c(0, -0.8, 0.6, 0, 0, -0.8, 0),
#   cex.edge.width = 0.8, edge.arrow.size = 1,
#   legend.prop = 0.3, ncol.legend = 2,
#   vertex.label.dist = 1.1, combine.slices = 0.02,
#   combined.slice.label = "others (emission prob. < 0.02)"
# )

## ----HMMplotLayout, fig.width=5.5, fig.height=4, dev.args=list(pointsize=10), echo=TRUE, fig.cap="Another example of \\texttt{plot.hmm}.", fig.align='center', fig.keep='last', cache = TRUE, message=FALSE----
require("igraph")
set.seed(1234)
plot(hmm_biofam,
  layout = layout_nicely, pie = FALSE,
  vertex.size = 30, vertex.label = "names", vertex.label.dist = 0,
  edge.curved = FALSE, edge.width = 1,
  loops = TRUE, edge.loop.angle = -pi / 8,
  trim = 0.01, label.signif = 3,
  xlim = c(-1, 1.3)
)

## ----colorpalette, fig.width=5.5, fig.height=3, dev.args=list(pointsize=10), echo=TRUE, fig.cap="Helper function for plotting colour palettes with their names.", fig.align='center', fig.keep='last', cache = TRUE----
plot_colors(colorpalette[[7]])

