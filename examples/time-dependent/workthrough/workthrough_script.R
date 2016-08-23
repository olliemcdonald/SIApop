# use install.packages('package-name') to install
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)

# Example from code to import and visualize plots


setwd('./examples/inhomogeneous/workthrough/')

tdat <- read_delim('tdata.txt', delim = "\t")
cdat <- read_delim('clonedata.txt', delim = "\t")

tdat %>% group_by(run, time) %>% summarize(totcells = sum(numcells)) %>%
  ggplot(aes(x = time, y = totcells, colour = as.factor(run))) + geom_line() +
  scale_y_log10()

tdat %>% filter(unique_id %in% c("a1.a", "a2.a", "a3.a", "a4.a", "a5.a")) %>%
  ggplot(aes(x = time, y = numcells, colour = as.factor(unique_id)))  +
  geom_line() + facet_wrap(~run, scales = "free") +
  theme(legend.position = "none")

subdat <- timedat %>% filter(numcells >= 5) %>% select(unique_id)
tdat %>% filter(unique_id %in% subdat$unique_id) %>%
  ggplot(aes(x = time, y = numcells, colour = as.factor(run)))  +
  geom_line() + facet_wrap(~unique_id, scales = "free") +
  theme(legend.position = "none")




tdat <- tdat %>% filter(numcells > 0)
clonesize <- 0 : (round(log10(max(tdat$numcells))) + 1)
tottime <- max(tdat$initialtime)
h <- ggplot(tdat , aes(x = initialtime, y = growthrate, group = unique_id, frame = time)) +
  geom_segment(aes(xend = parent_initialtime, yend = parent_growth_rate,
                   size = 0.5 * (tottime - initialtime) / tottime),
               col = "grey80", alpha = (0.5)) +
  geom_point(aes(size = log10(numcells), colour = log10(numcells))) +
  theme_bw(base_size = 24) + labs(x = "Time of appearance of subclone", y = "Growth rate") +
  xlim(0, max(tdat$time)) +
  scale_size_area(breaks = clonesize,
                  labels = as.character(10^clonesize)) +
  scale_color_gradientn(colours = rainbow(12)[1:9],
                        breaks = clonesize,
                        labels = as.character(10^clonesize)) +
  guides(color = guide_legend(), size = guide_legend()) +
  labs(colour = "Number\nof Cells", size = "Number\nof Cells") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~run)

h

############################################################
# For animations of evolution - WARNING: Takes a long time
library(gganimate)
gg_animate(h, filename = "bpanimation.mp4", saver = "mp4", interval = 0.2, ani.width = 1500, ani.height = 1000)
