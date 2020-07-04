library(nls.multstart)
library(broom)
library(tidyverse)
library(nlstools)
library(pdftools)
library(magick)
library(cowplot)



# To plot functions without data, specify range of x-axis
base <- ggplot(data.frame(x = c(-5, 7)), aes(x))
(ConflictwithF <- base + geom_rect(xmin =0, xmax = 1.5, ymin = 0, ymax = .42, fill = "grey", alpha = 0.3)+ stat_function(fun = dnorm, colour = I(rgb(.5,0,.5)), args = list(mean = 1.5, sd =1), size = 1.5) + stat_function(fun = dnorm, colour = I(rgb(0,0.5,0)), size = 1.5) +theme_cowplot() + coord_cartesian() + theme(axis.title = element_text(size = 20, face = "bold"), axis.text = element_blank(), axis.ticks = element_blank()) + labs(x = " ", y = "Fitness") +annotate("text", x = 0.75, y = 0.5, label =  "Conflict", size = 6, fontface = "italic")  +annotate("text", x = 0.75, y = 0.45, label =  "With fitness feedback", size = 4) + stat_function(fun = dnorm, colour = I(rgb(.5,0,.5)), args = list(mean = 3.5, sd =1), size = 1.5, alpha = 0.4)  + stat_function(fun = dnorm, colour = I(rgb(0,0.5,0)), args = list(mean = -1.75, sd =1), size = 1.5, alpha = 0.5) +  geom_segment(x = c(-1.25, 2.75), y =0.39, xend = c(-.3, 1.8), yend = 0.39, lineend = "round", linejoin = "mitre", size = 1, arrow = arrow(length = unit(0.13, "inches"))) + coord_cartesian(clip = "off"))

#+ geom_vline(xintercept = c(-1.75, 3.5), linetype = 3, size = .8) + geom_vline(xintercept = c(0, 1.5), linetype = 2, size = .8)
ggsave(file = "FitnessAlignmentFig.pdf")
(Conflict <- base + geom_rect(xmin =-1.75, xmax = 3.5, ymin = 0, ymax = .42, fill = "grey", alpha = 0.3)+theme_cowplot() + coord_cartesian() + theme(axis.title = element_text(size = 20, face = "bold"), axis.text = element_blank(), axis.ticks = element_blank()) + labs(x = " ", y = "Fitness") +annotate("text", x = 0.75, y = 0.5, label =  "Conflict", size = 6, fontface = "italic") +annotate("text", x = 0.75, y = 0.45, label =  "No fitness feedback", size = 4)  + stat_function(fun = dnorm, colour = I(rgb(.5,0,.5)), args = list(mean = 3.5, sd =1), size = 1.5)  + stat_function(fun = dnorm, colour = I(rgb(0,0.5,0)), args = list(mean = -1.75, sd =1), size = 1.5) + coord_cartesian(clip = "off"))
#pg <- ggplot_build(p) # If I wanted to access the values that were plotted
(AlignFig <-base + stat_function(fun = dnorm, colour = I(rgb(0,0.5,0)), args = list(mean = 0.75, sd =1), size = 1.5) + stat_function(fun = dnorm, args = list(mean = 0.75, sd =1.1), colour = I(rgb(.5,0,.5)), size = 1.5)  + theme_cowplot() + coord_cartesian() + theme(axis.title = element_text(size = 20, face = "bold"), axis.text = element_blank(), axis.ticks = element_blank()) + labs(x = " ", y = "Fitness") +annotate("text", x = 0.75, y = 0.5, label = c("Alignment"), size = 6, fontface = "italic")+annotate("text", x = c(-4,6), y = 0.2, label = c("Host", "Microbe"), colour = c(I(rgb(0,0.5,0)), I(rgb(.5,0,.5))), size = 8) + coord_cartesian(clip = "off"))#+annotate("text", x = c(-3,5.5), y = 0.2, label = c("Plant", "Microbe"), colour = c("blue", "red"), size = 10)

Concept <- image_read_pdf("~/Downloads/concept_fig_whose(1).pdf")# %>% image_resize("570x380")

p1 <- ggdraw() + draw_image(Concept, scale = 0.9)


# Haney data
FitDF <- read.csv("~/Documents/Projects/Azotobacter_Competition/fit_align.csv")

(FitDF <- FitDF %>% mutate(PlantWtFit = PlantWt/ max(PlantWt), LatRootsFit = LatRoots/ max(LatRoots), CFUFit = CFU/ max(CFU)))


# 
# use nls.multstart to fit optima

FitDF2 <- FitDF %>% select(PlantWtFit, LatRoots, CFUFit) %>% pivot_longer(- LatRoots, names_to = "Species", values_to = "RelativeFit") %>% mutate(Species= recode(Species, PlantWtFit = "Plant", CFUFit = "Microbe"))
ggplot(FitDF2, aes(x = LatRoots, y = RelativeFit, colour = Species)) + geom_point(size =3)

# Define the function:
GausFit <- function(LatRoots, optimum, sd.fit){
  RelativeFit <- exp(-1*((LatRoots-optimum)^2)/(2*(sd.fit^2)))
  return(RelativeFit)
}
# 
# # try only plant first
# fit <- nls_multstart(RelativeFit ~ GausFit(LatRoots, optimum, sd.fit),
#                      data = subset(FitDF2, Species == "Microbe"),
#                      iter = 250,
#                      start_lower = c(optimum = 40, sd.fit = 30),
#                      start_upper = c(optimum = 50, sd.fit = 35),
#                      supp_errors = 'Y',
#                      convergence_count = 100,
#                      na.action = na.omit)
# 
# (info <- glance(fit))
# (params <- tidy(fit))
# (preds <- augment(fit))
# 
# ggplot() +
#   geom_point(aes(LatRoots, RelativeFit), subset(FitDF2, Species == "Plant")) +
#   geom_line(aes(LatRoots, .fitted), preds)
# 
# 
# Lrs <- seq(10,100)
# new_preds <- data.frame(LatRoots = Lrs, re2 = GausFit(Lrs, 53.9, 32.9))
# ggplot(new_preds, aes(x = LatRoots, y = re2)) + geom_point(size =3)
                          

# lets try with both plant and microbe

fits <- FitDF2 %>%
  group_by(., Species) %>%
  nest()%>%
  mutate(fit= purrr::map(data, ~ nls_multstart(RelativeFit ~ GausFit(LatRoots, optimum, sd.fit),
                                               data = .x,
                                               iter = 1000,
                                               start_lower = c(optimum = 45, sd.fit = 30),
                                               start_upper = c(optimum = 100, sd.fit = 100),
                                               supp_errors = 'N',
                                               na.action = na.omit)))

info <- fits %>%
  mutate(summary = map(fit, glance)) %>%
  unnest(summary)

params <- fits %>%
  mutate(., p = map(fit, tidy)) %>%
  unnest(p)

CI <- fits %>%
  mutate(., cis = map(fit, confint2),
         cis = map(cis, data.frame)) %>%
  unnest(cis) %>%
  rename(., conf.low = X2.5.., conf.high = X97.5..) %>%
  group_by(., Species) %>%
  mutate(., term = c('optimum', 'sd.fit')) %>%
  ungroup() %>%
  select(., -data, -fit)

params <- merge(params, CI, by = intersect(names(params), names(CI)))

# get predictions
preds <- fits %>%
  mutate(., p = map(fit, augment)) %>%
  unnest(p)

new_preds <- data.frame(LatRoots = seq(0, 200))

preds2 <- fits %>%
  mutate(., p = map(fit, augment, newdata = new_preds)) %>%
  unnest(p) %>%
  rename(., RelativeFit = .fitted)

(Haney <-ggplot()+  geom_rect(aes(xmin =53.8549, xmax = 67.5758, ymin = 0, ymax = 1.01), fill = "grey", alpha = 0.5)+geom_point(aes(LatRoots, RelativeFit, colour = Species), size = 1, FitDF2) +geom_line(aes(LatRoots, RelativeFit, colour = Species), size = 2, preds2)  + scale_colour_manual(values = c( "#800080","#008000")) +theme_cowplot() + labs(x = "Trait", y = "Fitness") + theme(axis.title = element_text(size = 20, face = "bold"), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "none"))

#  # ggplot(params, aes(col = Species)) +
#  #  geom_point(aes(Species, estimate), size = 1) +
#  #  facet_wrap(~ term, scale = 'free_x', ncol = 4) +
#  #  geom_linerange(aes(Species, ymin = conf.low, ymax = conf.high)) + coord_flip()
#      
# RightSide <- plot_grid(AlignFig,Conflict, Haney, ncol = 1)         
# plot_grid(p1, RightSide, ncol = 2)
# ggsave("Fig1Draft3.pdf")
# 
# # Draft 4 pattern on R: 1,2,1
# Mid2panel <- plot_grid(Conflict, ConflictwithF, ncol = 2, labels = c("C", "D"))
# RightSide2 <- plot_grid(AlignFig, Mid2panel, Haney, ncol = 1, labels = c("B", "", "E"))
# plot_grid(p1, RightSide2, ncol = 2, labels = c("A", ""))
# ggsave("Fig1Draft4.pdf")
# 
# # Draft 5 pattern on R: 3,1
# Toppanel <- plot_grid(AlignFig,Conflict, ConflictwithF, ncol = 3, labels = c("B", "C", "D"))
# RightSide3 <- plot_grid(Toppanel, Haney, ncol = 1, labels = c("", "E"))
# plot_grid(p1, RightSide3, ncol = 2, labels = c("A", ""))
# ggsave("Fig1Draft5.pdf")
# 
# # Draft 6 pattern on R: 2,2
# Toppanel <- plot_grid(AlignFig, Conflict, ncol = 2, labels = c("B","C" ))
# Bottom <- plot_grid(ConflictwithF, Haney, ncol = 2, labels = c("D", "E"))
# RightSide4 <- plot_grid(Toppanel, Bottom, ncol = 1)
# plot_grid(p1, RightSide4, ncol = 2, labels = c("A", ""))
# ggsave("Fig1Draft6.pdf")

# Draft 7 pattern on R: 1
RightSide5 <- plot_grid(AlignFig,Conflict, ConflictwithF, Haney, ncol = 1, labels = c("B", "C", "D", "E"))
FullPlot <-plot_grid(p1, RightSide5, ncol = 2, labels = c("A", ""))
ggsave("Fig1Draft8.pdf")




