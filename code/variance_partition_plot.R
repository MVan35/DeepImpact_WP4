## ---------------------------------------------------------------------------
## variance_partition_plot() — graphical-abstract figure from adonis2 objects
##
##
## Usage in a .qmd chunk:
##   ```{r fig-abstract, fig.width=7, fig.height=1.9, dev='png', dev.args=list(bg='transparent')}
##   source("variance_partition_plot.R")
##   pdat <- partition_from_adonis(bacteria_sequential, fungi_sequential)
##   variance_partition_plot(pdat)
##   ```
## ---------------------------------------------------------------------------

library(ggplot2)
library(dplyr)

## Palette — keep in sync with the abstract artwork -------------------------
VP_COLS <- c(
  "Soil legacy"        = "#6B4A22",
  "Neighbour identity" = "#2E5B2E",
  "Soil x neighbour"   = "#3B6FB6",
  "Unexplained"        = "#C3C4C0"
)


## -------------------------------------------------------------------------
## 1. Pull R2 out of adonis2 objects into a tidy frame
## -------------------------------------------------------------------------
partition_from_adonis <- function(bact_fit, fung_fit,
                                  labels = c(Bacteria = "rpoB", Fungi = "ITS1")) {

  one <- function(fit, kingdom) {
    d <- as.data.frame(fit)
    d$term <- rownames(d)

    # adonis2 row names vary with the formula: "soil:species" for a * b, and
    # vegan has used both "Residual" and "Residuals" across versions.
    pick <- function(pattern, what) {
      hit <- grep(pattern, d$term, ignore.case = TRUE, value = TRUE)[1]
      if (is.na(hit)) {
        # by = "margin" on a `soil * species` 
        hint <- if (!any(grepl("^soil$", d$term)) && any(grepl(":", d$term)))
          paste0("\n  This looks like a by = \"margin\" fit: marginal tests of a ",
                 "`soil * species` model\n  report only the interaction, so there is no ",
                 "partition to plot. Refit with\n  by = \"terms\" (sequential), which is what ",
                 "the manuscript reports:\n    adonis2(d ~ soil + species + soil:species, ",
                 "data = m, by = \"terms\")")
        else ""
        stop("no '", what, "' row in the adonis2 output for ", kingdom,
             ".\n  Rows present: ", paste(d$term, collapse = ", "), hint,
             call. = FALSE)
      }
      d$R2[d$term == hit] * 100
    }

    out <- data.frame(
      kingdom   = kingdom,
      component = factor(names(VP_COLS), levels = names(VP_COLS)),
      pct       = c(pick("^soil$", "soil"),
                    pick("^species$", "species"),
                    pick("^soil:species$|^species:soil$", "soil:species"),
                    pick("^residual", "Residual")),
      stringsAsFactors = FALSE
    )

    tot <- sum(out$pct)
    if (abs(tot - 100) > 0.5)
      warning(kingdom, ": components sum to ", round(tot, 2),
              " %, not 100 %. Did you pass a by='margin' fit? ",
              "Use by='terms' (sequential) for a partition.")
    out
  }

  res <- rbind(one(bact_fit, "Bacteria"), one(fung_fit, "Fungi"))

  # facet strip text: "Bacteria (rpoB)" with the gene in italics -> plotmath
  res$facet <- factor(
    ifelse(res$kingdom == "Bacteria",
           paste0("bold('Bacteria')~~italic('", labels[["Bacteria"]], "')"),
           paste0("bold('Fungi')~~italic('",   labels[["Fungi"]],    "')")),
    levels = c(paste0("bold('Bacteria')~~italic('", labels[["Bacteria"]], "')"),
               paste0("bold('Fungi')~~italic('",   labels[["Fungi"]],    "')"))
  )
  res
}


## -------------------------------------------------------------------------
## 2. Draw it
## -------------------------------------------------------------------------
## pdat       : output of partition_from_adonis()
## inside_min : segments at least this wide (%) get their label inside;
##              narrower ones are pushed to the left of the bar
## show_key   : FALSE for the abstract (the artwork carries the key),
##              TRUE for standalone use (slide, supplement)
## base_size  : all text scales from this
## -------------------------------------------------------------------------
variance_partition_plot <- function(pdat,
                                    inside_min = 6,
                                    show_key   = FALSE,
                                    base_size  = 9,
                                    label_fmt  = "%.1f %%") {

  d <- pdat %>%
    group_by(facet) %>%
    arrange(component, .by_group = TRUE) %>%
    mutate(
      # stack left-to-right in component order
      right  = cumsum(pct),
      left   = right - pct,
      centre = (left + right) / 2,
      inside = pct >= inside_min,
      lab    = sprintf(label_fmt, pct)
    ) %>%
    ungroup()

  p <- ggplot(d) +
    geom_rect(aes(xmin = left, xmax = right, ymin = -0.5, ymax = 0.5,
                  fill = component),
              colour = "white", linewidth = 0.4) +
    # labels that fit inside their segment
    geom_text(data = filter(d, inside),
              aes(x = centre, y = 0, label = lab,
                  colour = component == "Unexplained"),
              size = base_size * 0.30, fontface = "bold", show.legend = FALSE) +
    # labels too narrow to sit inside: park them left of the bar
    geom_text(data = filter(d, !inside),
              aes(x = -1.5, y = 0, label = lab, colour = TRUE),
              hjust = 1, size = base_size * 0.30, fontface = "bold",
              show.legend = FALSE) +
    scale_fill_manual(values = VP_COLS, name = NULL) +
    scale_colour_manual(values = c("TRUE" = "#3A3B38", "FALSE" = "white")) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    facet_wrap(~ facet, ncol = 1, labeller = label_parsed) +
    coord_cartesian(xlim = c(-9, 100), ylim = c(-0.62, 0.62), clip = "off") +
    theme_void(base_size = base_size) +
    theme(
      strip.text      = element_text(hjust = 0, size = base_size,
                                     margin = margin(b = 3, t = 5)),
      panel.spacing.y = unit(0.9, "lines"),
      legend.position = if (show_key) "bottom" else "none",
      legend.key.size = unit(0.42, "lines"),
      legend.text     = element_text(size = base_size * 0.8),
      legend.margin   = margin(t = 2),
      plot.margin     = margin(2, 4, 2, 2),
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA)
    ) +
    guides(fill = guide_legend(nrow = 1))

  p
}


## -------------------------------------------------------------------------
## 3. Saving with a genuinely transparent background
## -------------------------------------------------------------------------
save_variance_partition <- function(p, stem = "FigAbstract_variance_partition",
                                    width = 7, height = 1.9, dpi = 400,
                                    show = TRUE) {
  ggsave(paste0(stem, ".png"), p, width = width, height = height,
         dpi = dpi, bg = "transparent")
  ggsave(paste0(stem, ".pdf"), p, width = width, height = height,
         bg = "transparent", device = cairo_pdf)
  if (show) print(p)
  invisible(p)
}


## -------------------------------------------------------------------------
## 4. Displaying it
## -------------------------------------------------------------------------
## variance_partition_plot() returns a ggplot object, which auto-prints ONLY at
## the top level of the console or a chunk:
##
##   variance_partition_plot(pdat)          # draws
##   p <- variance_partition_plot(pdat)     # assigns, draws nothing
##   p                                      # now draws
##   print(p)                               # always draws — use this inside
##                                          # for(), sapply(), or a function,
##                                          # where auto-printing is suppressed
## -------------------------------------------------------------------------

