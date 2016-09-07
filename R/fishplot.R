# plot functions for


# pops off everything after the last instance of pattern
# in a vector
pop <- function(string, pattern = ">"){
  string <- unlist(strsplit(string, pattern))
  string <- string[1:length(string)-1]
  paste(string, collapse = pattern)
}

# matches clone index to their respective
# parent clone index
match.parents <- function(ids){
  parents <- unlist(lapply(ids, pop, ">"))
  parents <- match(parents, ids)
  ifelse(is.na(parents), 0, parents)

}

fishPlotBD <- function(tdat, allele.freq, shape = "polygon", vlines = 0, col.vline = "#FFFFFF99",
                      vlab = NULL, border = 0.5, col.border = "#777777", pad.left = 0.0,
                      ramp.angle = 0.5, title = NULL, title.btm = NULL, cex.title = NULL,
                      cex.vlab = 0.7){
  require(dplyr)
  require(tidyr)
  require(fishplot)

  ancestors <- (tdat %>% filter(time == 0))$unique_id
  topallele <- tdat %>% filter(unique_id %in% ancestors) %>% group_by(time) %>%
    summarize(timefreq = sum(allelefreq)) %>%
    summarize(maxfreq = max(timefreq))

  maxalleles <- (tdat %>% group_by(unique_id) %>%
    summarize(maxallelefreq = max(allelefreq)) %>%
    filter(maxallelefreq >= allele.freq * topallele$maxfreq))$unique_id

  subdat <- tdat %>% filter(unique_id %in% maxalleles)
  subdat <- subdat %>% mutate(allele_frac = allelefreq / topallele$maxfreq)

#  subdat <- subdat %>% filter((time * 100) %% 10 == 0)
  subdat_wide <- subdat %>% select(time, unique_id, allele_frac) %>%
    spread(time, allele_frac, fill = 0)

  frac.table <- as.matrix(subdat_wide[,2:ncol(subdat_wide)]) * 100
  timepoints <- as.numeric(colnames(subdat_wide[,2:ncol(subdat_wide)]))

  parents <- match.parents(maxalleles)

  fish <- createFishObject(frac.table,parents,timepoints=timepoints)
  fish <- setCol(fish, rainbow(length(unique(maxalleles))))
  #calculate the layout of the drawing
  fish <- layoutClones(fish)

  #draw the plot, using the splining method (recommended)
  #and providing both timepoints to label and a plot title
  fishPlot(fish, shape = "polygon", vlines, col.vline, vlab, border,
           col.border, pad.left, ramp.angle, title, title.btm,
           cex.title, cex.vlab)
}
