#####################################################################
###   Targeted gene flow under a changing environment.            ###
###   Modified from generic TGF model by Kelly & Phillips 2018    ###  
###   Mods by Adam Smart                                          ###
###   For publication in Smart & Phillips 2020. "X"               ###
###                                                               ###
###   Updated 25.02.2020                                          ###
#####################################################################

### Other scripts required for full analysis: 'spartan_TGF_parameters.R; spartan_TGF_dataClean.R'

### Load in required packages
.libPaths("/home/asmart1/R/lib/3.4")
lib<- .libPaths()[1]
repo<- "https://cran.ms.unimelb.edu.au/"
pkgs = c("tidyverse", "viridis", "purrr", "ggplot2", "gridExtra", "wesanderson")
new.pkgs <- pkgs[!(pkgs %in% installed.packages(lib=lib)[,"Package"])]
if(length(new.pkgs)) install.packages(new.pkgs, lib=lib, repos=repo)
inst = lapply(pkgs, library, character.only = TRUE)

### Set path and read in output files from Spartan
path = "out/output/"
output.RDS <- list.files(path = path, pattern = "\\.rds$", full.names = TRUE)
summaryStack <- lapply(output.RDS, readRDS)

### Create reference grid for title
initFreq_options = c(0.05, 0.3)
recomb_options = c(10, 50)
h2_options = c(0.1, 0.2, 0.3)
Rmax_options = c(2, 3, 5)
Nstar_options = c(500, 1000, 2000)
m_options = c(0.15,0.3,1)
decay_options = 1
alpha_options = c(0, 0.04, 0.277)

# Create grid of all possible parameter combinations
testGrid<-expand.grid(initFreq_options,
                      recomb_options,
                      h2_options, 
                      Rmax_options,
                      Nstar_options,
                      m_options, 
                      decay_options,
                      alpha_options)

# rename & print testGrid
colnames(testGrid)<-c("initFreq_options", 
                      "recomb_options",
                      "h2_options", 
                      "Rmax_options",
                      "Nstar_options",
                      "m_options", 
                      "decay_options",
                      "alpha_options"
                      )

### test case input for manuscript figures c(1,1,2,2,2,2,1)
init_sel<-initFreq_options[1]
recomb_sel<-recomb_options[1]
h2_sel<-h2_options[1]
Rmax_sel<-Rmax_options[2]
Nstar_sel<-Nstar_options[1]
m_sel<- m_options[2]
alpha_sel<-alpha_options[1]

testCase_ID<-subset(testGrid, Rmax_options==Rmax_sel & recomb_options==recomb_sel & h2_options==h2_sel & alpha_options==alpha_sel & initFreq_options==init_sel & m_options == m_sel & Nstar_options==Nstar_sel)
print(testCase_ID)

### Convert list of matricies to list of data.frames
TGFstack<- map(summaryStack, as.data.frame)

#----------------------- Manuscript Plots -------------------------------

#----------------------- Figure 1 -------------------------------
{
tile1<- ggplot((TGFstack[[301]]), aes(x=props, y=k, scale=props)) +
  geom_tile(aes(fill =extinctRisk)) +
  ylim(1, 51) +
  scale_fill_viridis(option = "C", begin = 0 , end = 1, na.value = "black", limits = c(0, 1), direction = -1) +
  theme_classic() +
  labs(title = "(a)", fill= "Probability of extinction (x)", x="Proportion of pre-adapted individuals introduced", y=expression(atop("Timing of TGf action (years) in relation", paste("to environmental shift")))) +
  theme(legend.position=c("bottom"), legend.title = element_text(size=10), legend.text = element_text(size=7))

tile2<- ggplot((TGFstack[[301]]), aes(x=props, y=k, scale=props)) +
  geom_tile(aes(fill =simpsonsDI)) +
  ylim(1, 51) +
  scale_fill_viridis(option = "D", begin = 0, end = 1, na.value = "black", limits = c(0, 1), direction = 1) +
  theme_classic()+
  labs(title= "(b)", fill = expression(Gini-Simpson~Diversity~Index~(D[GS])), x="Proportion of pre-adapted individuals introduced", y=expression(atop("Timing of TGf action (years) in relation", paste("to environmental shift")))) +
  theme(legend.position=c("bottom"), legend.title = element_text(size=10), legend.text = element_text(size=7))

tile3<- ggplot((TGFstack[[301]]), aes(x=props, y=k, scale=props)) +
  geom_tile(aes(fill =(1-extinctRisk)*simpsonsDI)) +
  ylim(1, 51) +
  scale_fill_viridis(option = "A", begin = 0 , end = 1, na.value = "black", limits = c(0, 1)) +
  theme_classic() +
  stat_contour(breaks= c(0.5, 0.9), aes(z=(1-extinctRisk)*simpsonsDI), color = "white", size = 0.6) +
  labs(title = "(c)", fill= expression(Expected~Return~(E[Y])), x="Proportion of pre-adapted individuals introduced", y=expression(atop("Timing of TGf action (years) in relation", paste("to environmental shift")))) +
  theme(legend.position=c("bottom"), legend.title = element_text(size=10), legend.text = element_text(size=7))

grob<- arrangeGrob(tile1, tile2, tile3, ncol=1)
ggsave(file="man/images/a4_figure1.pdf", grob, width = 110, height = 297, units = "mm")
}
#----------------------- Figure 3 -------------------------------
{
tile1<- ggplot((TGFstack[[301]]), aes(x=props, y=k, scale=props)) +
  geom_tile(aes(fill =(1-extinctRisk)*simpsonsDI)) +
  ylim(1, 51) +
  scale_fill_viridis(option = "A", begin = 0 , end = 1, na.value = "black", limits = c(0, 1), direction = 1) +
  theme_classic() +
  stat_contour(breaks= c(0.5, 0.9), aes(z=(1-extinctRisk)*simpsonsDI), color = "white", size = 0.6) +
  labs(title = "(a)", fill= expression(Expected~Return~(E[Y])), x="Proportion of pre-adapted individuals introduced", y=expression(atop("Timing of TGf action (years) in relation", paste("to environmental shift")))) +
  theme(legend.position=c("bottom"), legend.title = element_text(size=10), legend.text = element_text(size=7))

tile2<- ggplot((TGFstack[[625]]), aes(x=props, y=k, scale=props)) +
  geom_tile(aes(fill =(1-extinctRisk)*simpsonsDI)) +
  ylim(1, 51) +
  scale_fill_viridis(option = "A", begin = 0, end = 1, na.value = "black", limits = c(0, 1), direction = 1) +
  theme_classic()+
  stat_contour(breaks= c(0.5, 0.9), aes(z=(1-extinctRisk)*simpsonsDI), color = "white", size = 0.6) +
  
  labs(title= "(b)", fill= expression(Expected~Return~(E[Y])), x="Proportion of pre-adapted individuals introduced", y=expression(atop("Timing of TGf action (years) in relation", paste("to environmental shift")))) +
  theme(legend.position=c("bottom"), legend.title = element_text(size=10), legend.text = element_text(size=7))

tile3<- ggplot((TGFstack[[229]]), aes(x=props, y=k, scale=props)) +
  geom_tile(aes(fill =(1-extinctRisk)*simpsonsDI)) +
  ylim(1, 51) +
  scale_fill_viridis(option = "A", begin = 0 , end = 1, na.value = "black", limits = c(0, 1), direction = 1) +
  theme_classic() +
  stat_contour(breaks= c(0.5, 0.9), aes(z=(1-extinctRisk)*simpsonsDI), color = "white", size = 0.6) +
  labs(title = "(c)", fill= expression(Expected~Return~(E[Y])), x="Proportion of pre-adapted individuals introduced", y=expression(atop("Timing of TGf action (years) in relation", paste("to environmental shift")))) +
  theme(legend.position=c("bottom"), legend.title = element_text(size=10), legend.text = element_text(size=7))

grob<- arrangeGrob(tile1, tile2, tile3, ncol=1)
ggsave(file="man/images/a4_figure3.pdf", grob, width = 110, height = 297, units = "mm")
}

#----------------------- Figure 6 -------------------------------

{
  tile1<- ggplot((TGFstack[[13]]), aes(x=props, y=k, scale=props)) +
    geom_tile(aes(fill =(1-extinctRisk)*simpsonsDI)) +
    ylim(1, 51) +
    scale_fill_viridis(option = "A", begin = 0 , end = 1, na.value = "black", limits = c(0, 1), direction = 1) +
    theme_classic() +
    stat_contour(breaks= c(0.5, 0.9), aes(z=(1-extinctRisk)*simpsonsDI), color = "white", size = 0.6) +
    labs(title = "(a)", fill= expression(Expected~Return~(E[Y])), x="Proportion of pre-adapted individuals introduced", y=expression(atop("Timing of TGf action (years) in relation", paste("to environmental shift")))) +
    theme(legend.position=c("bottom"), legend.title = element_text(size=10), legend.text = element_text(size=7))
  
  tile2<- ggplot((TGFstack[[121]]), aes(x=props, y=k, scale=props)) +
    geom_tile(aes(fill =(1-extinctRisk)*simpsonsDI)) +
    ylim(1, 51) +
    scale_fill_viridis(option = "A", begin = 0, end = 1, na.value = "black", limits = c(0, 1), direction = 1) +
    theme_classic()+
    stat_contour(breaks= c(0.5, 0.9), aes(z=(1-extinctRisk)*simpsonsDI), color = "white", size = 0.6) +
    
    labs(title= "(b)", fill= expression(Expected~Return~(E[Y])), x="Proportion of pre-adapted individuals introduced", y=expression(atop("Timing of TGf action (years) in relation", paste("to environmental shift")))) +
    theme(legend.position=c("bottom"), legend.title = element_text(size=10), legend.text = element_text(size=7))
  
  tile3<- ggplot((TGFstack[[229]]), aes(x=props, y=k, scale=props)) +
    geom_tile(aes(fill =(1-extinctRisk)*simpsonsDI)) +
    ylim(1, 51) +
    scale_fill_viridis(option = "A", begin = 0 , end = 1, na.value = "black", limits = c(0, 1), direction = 1) +
    theme_classic() +
    stat_contour(breaks= c(0.5, 0.9), aes(z=(1-extinctRisk)*simpsonsDI), color = "white", size = 0.6) +
    labs(title = "(c)", fill= expression(Expected~Return~(E[Y])), x="Proportion of pre-adapted individuals introduced", y=expression(atop("Timing of TGf action (years) in relation", paste("to environmental shift")))) +
    theme(legend.position=c("bottom"), legend.title = element_text(size=10), legend.text = element_text(size=7))
  
  grob<- arrangeGrob(tile1, tile2, tile3, ncol=1)
  ggsave(file="man/images/a4_figure6.pdf", grob, width = 110, height = 297, units = "mm")
}
#----------------------- Figure 2, 4 & 5 -------------------------------
testGrid$extinct<- rep(NA, times= length(TGFstack))
testGrid$GSI<- rep(NA, times= length(TGFstack))
testGrid$EY<- rep(NA, times= length(TGFstack))
testGrid$k<- rep(NA, times= length(TGFstack))
testGrid$props<- rep(NA, times= length(TGFstack))
testGrid$EY_percent <- rep(NA, times= length(TGFstack))
testGrid$k_percent<- rep(NA, times= length(TGFstack))
testGrid$props_percent<- rep(NA, times= length(TGFstack))

for (i in 1:nrow(testGrid)){
  testGrid[i,]$extinct<- apply(TGFstack[[i]][2], 2, min)
}

for (i in 1:nrow(testGrid)){
  testGrid[i,]$GSI<- apply(TGFstack[[i]][5], 2, max)
}

maxFind<- function(TGF){
  ER<-(1-TGF[[2]])*TGF[[5]]
  key <- which(ER == max(ER))
  
  if (length(key) > 1) {
    col_vals <- TGF[[4]][key]
    key <- key[which(col_vals==min(col_vals))]
  }
  
  if (length(key) > 1) {
    col_vals <- TGF[[1]][key]
    key <- key[which(col_vals==min(col_vals))]
  }
  key
}

maxFind_percent<- function(TGF){
  ER<-(1-TGF[[2]])*TGF[[5]]
  key <- which(ER >= max(ER)*0.9)
  
  if (length(key) > 1) {
    col_vals <- TGF[[4]][key]
    key <- key[which(col_vals==min(col_vals))]
  }
  
  if (length(key) > 1) {
    col_vals <- TGF[[1]][key]
    key <- key[which(col_vals==min(col_vals))]
  }
  key
}

for (i in 1:nrow(testGrid)){
  TGF <- TGFstack[[i]]
  max<- maxFind(TGF)
  max_percent<-maxFind_percent(TGF)
  temp<- TGFstack[[i]][max,]
  temp_percent<- TGFstack[[i]][max_percent,]
  testGrid[i,]$k<- TGFstack[[i]][max,1]
  testGrid[i,]$props <- TGFstack[[i]][max, 4]
  testGrid[i,]$EY<- (1-temp[1,2])*temp[1,5]
  testGrid[i,]$k_percent<- TGFstack[[i]][max_percent,1]
  testGrid[i,]$props_percent <- TGFstack[[i]][max_percent, 4]
  testGrid[i,]$EY_percent<- (1-temp_percent[1,2])*temp_percent[1,5]
}

head(testGrid)

# Set plotting paramters
plotWidth = 210
plotHeight= 180

## -------------------- NSTAR ----------------------------------
{
pal<- wes_palette("FantasticFox1", n=3)

Rmax_sel<-Rmax_options[1]
recomb_sel<-recomb_options[1]
alpha_sel<-alpha_options[1]
init_sel<-initFreq_options[1]



testGrid_select<-subset(testGrid,
                          Rmax_options==Rmax_sel & 
                          recomb_options==recomb_sel &
                          alpha_options==alpha_sel & 
                          initFreq_options==init_sel)



panel1<- ggplot(testGrid_select, aes(x=factor(h2_options), y=EY)) +
  geom_point(aes(colour=factor(Nstar_options), shape = factor(h2_options)), size = 4) +
  geom_line(aes(colour= factor(Nstar_options), group = factor(Nstar_options))) +
  scale_shape_manual(values = c(16, 17, 15)) +
  scale_color_manual(values = pal) +
  ylim(0, 1.01) +
  labs(title = "(a)", y= expression(Expected~Return~(E[Y])), x=expression(Heritability~(h^2)), color = expression(atop("Carrying", paste("capacity"))), shape=expression(Heritability~(h^2))) +
  theme(panel.background = element_rect(fill = "white"), 
        panel.border = element_rect(fill=NA), 
        axis.line = element_blank())
  
grob3<- panel1 + facet_wrap(~ m_options, nrow = 1) + 
  theme(legend.position = "right",
        legend.box = "vertical",
        legend.box.background = element_blank(),
        strip.background =element_rect(colour = "black", fill="lightgrey"),
        legend.key = element_rect(colour = "white", fill = "white", size = 0.25)) 

panel2<- ggplot(testGrid_select, aes(x=(props_percent), y=(k_percent),)) +
  geom_point(aes(colour=factor(Nstar_options), shape = factor(h2_options)), size = 5) +
  xlim(-0.1, 0.51)+
  scale_shape_manual(values = c(1, 2, 22)) +
  scale_color_manual(values = pal) +
  ylim(0, 51) +
  labs(title = "(b)", y=expression(atop("Timing of TGf action in relation", paste("to environmental shift (years)"))), x="Proportion of pre-adapted individuals introduced", shape=expression(Heritability~(h^2)), color = expression(atop("Carrying", paste("capacity")))) +
  theme(panel.background = element_rect(fill = "white"), 
        panel.border = element_rect(fill=NA), 
        axis.line = element_blank())

grob4<- panel2 + facet_wrap(~ m_options, nrow = 1) + 
  theme(legend.position = "right",
        legend.box = "vertical",
        legend.box.background = element_blank(),
        strip.background =element_rect(colour = "black", fill="lightgrey"),
        legend.key = element_rect(colour = "white", fill = "white", size = 0.25)) 

grob5<- arrangeGrob(grob3, grob4, ncol =1, nrow=2)
ggsave(file="man/images/a4_figure2.pdf", grob5, width = plotWidth, height = plotHeight, units = "mm")
}
## -------------------- Outbreeding depression plot ----------------------------------
{
pal<- wes_palette("FantasticFox1", n=3)

#subset data to contain only m_options c(1, 3, 6)

Rmax_sel<-Rmax_options[1] 
Nstar_sel<-Nstar_options[2] 
init_sel<-initFreq_options[1] 
recomb_sel<-recomb_options[1]

testGrid_select<-subset(testGrid, 
                        Rmax_options==Rmax_sel & 
                        Nstar_options==Nstar_sel &
                        initFreq_options==init_sel & 
                        recomb_options==recomb_sel)

panel3<- ggplot(testGrid_select, aes(x=factor(h2_options), y=EY)) +
  geom_point(aes(colour=factor(alpha_options), shape = factor(h2_options)), size = 4) +
  geom_line(aes(colour= factor(alpha_options), group = factor(alpha_options))) +
  scale_shape_manual(values = c(16, 17, 15)) +
  scale_color_manual(values = pal, labels = c("0", "10", "50")) +
  ylim(0, 1.05) +
  labs(title = "(a)", y= expression(Expected~Return~(E[Y])), x=expression(Heritability~(h^2)), shape=expression(Heritability~(h^2)), color=expression(atop("Outbreeding", paste("depression")))) +
  theme(panel.background = element_rect(fill = "white"), 
        panel.border = element_rect(fill=NA), 
        axis.line = element_blank())

grob6<- panel3 + facet_wrap(~ m_options, nrow = 1) + 
  theme(legend.position = "right",
        legend.box = "vertical",
        legend.box.background = element_blank(),
        strip.background =element_rect(colour = "black", fill="lightgrey"),
        legend.key = element_rect(colour = "white", fill = "white", size = 0.25)) 

panel4<- ggplot(testGrid_select, aes(x=(props_percent), y=(k_percent))) +
  geom_point(aes(colour=factor(alpha_options), shape = factor(h2_options)), size = 4) +
  xlim(-0.1, 0.51)+
  scale_shape_manual(values = c(1, 2, 22)) +
  scale_color_manual(values = pal, labels = c("0", "10", "50")) +
  ylim(0, 51) +
  labs(title = "(b)", y=expression(atop("Timing of TGf action in relation", paste("to environmental shift (years)"))), x="Proportion of pre-adapted individuals introduced", shape=expression(Heritability~(h^2)), color=expression(atop("Outbreeding", paste("depression")))) +
  theme(panel.background = element_rect(fill = "white"), 
        panel.border = element_rect(fill=NA), 
        axis.line = element_blank())

grob7<- panel4 + facet_wrap(~ m_options, nrow = 1) + 
  theme(legend.position = "right",
        legend.box = "vertical",
        legend.box.background = element_blank(),
        strip.background =element_rect(colour = "black", fill="lightgrey"),
        legend.key = element_rect(colour = "white", fill = "white", size = 0.25)) 

grob8<- arrangeGrob(grob6, grob7, ncol =1, nrow=2)
ggsave(file="man/images/a4_figure4.pdf", grob8, width = plotWidth, height = plotHeight, units = "mm")
}
## -------------------- recomb ----------------------------------
{
pal<- wes_palette("FantasticFox1", n=2)

Rmax_sel<-Rmax_options[1] 
alpha_sel<-alpha_options[1] 
initFreq_sel<-initFreq_options[1] 
Nstar_sel<-Nstar_options[2]

testGrid_select<-subset(testGrid, 
                        Rmax_options==Rmax_sel &
                        alpha_options==alpha_sel &
                        initFreq_options==initFreq_sel &
                        Nstar_options==Nstar_sel)

panel1<- ggplot(testGrid_select, aes(x=factor(h2_options), y=EY)) +
  geom_point(aes(colour=factor(recomb_options), shape = factor(h2_options)), size = 4) +
  geom_line(aes(colour= factor(recomb_options), group = factor(recomb_options))) +
  scale_shape_manual(values = c(16, 17, 15)) +
  scale_color_manual(values = pal) +
  ylim(0, 1.01) +
  labs(title = "(a)", y= expression(Expected~Return~(E[Y])), x=expression(Heritability~(h^2)), color = expression(atop("Recombination", paste("rate"))), shape=expression(Heritability~(h^2))) +
  theme(panel.background = element_rect(fill = "white"), 
        panel.border = element_rect(fill=NA), 
        axis.line = element_blank())

grob3<- panel1 + facet_wrap(~ m_options, nrow = 1) + 
  theme(legend.position = "right",
        legend.box = "vertical",
        legend.box.background = element_blank(),
        strip.background =element_rect(colour = "black", fill="lightgrey"),
        legend.key = element_rect(colour = "white", fill = "white", size = 0.25)) 

panel2<- ggplot(testGrid_select, aes(x=(props_percent), y=(k_percent),)) +
  geom_point(aes(colour=factor(recomb_options), shape = factor(h2_options)), size = 4) +
  xlim(-0.1, 0.51)+
  scale_shape_manual(values = c(1, 2, 22)) +
  scale_color_manual(values = pal) +
  ylim(0, 51) +
  labs(title = "(b)", y=expression(atop("Timing of TGf action in relation", paste("to environmental shift (years)"))), x="Proportion of pre-adapted individuals introduced", shape=expression(Heritability~(h^2)), color = expression(atop("Recombination", paste("rate")))) +
  theme(panel.background = element_rect(fill = "white"), 
        panel.border = element_rect(fill=NA), 
        axis.line = element_blank())

grob4<- panel2 + facet_wrap(~ m_options, nrow = 1) + 
  theme(legend.position = "right",
        legend.box = "vertical",
        legend.box.background = element_blank(),
        strip.background =element_rect(colour = "black", fill="lightgrey"),
        legend.key = element_rect(colour = "white", fill = "white", size = 0.25)) 

grob5<- arrangeGrob(grob3, grob4, ncol =1, nrow=2)
ggsave(file="man/images/a4_figure5.pdf", grob5, width = plotWidth, height = plotHeight, units = "mm")
}

## -------------------- Rmax ----------------------------------
{
  pal<- wes_palette("FantasticFox1", n=3)
  
  init_sel<-initFreq_options[1]  
  alpha_sel<-alpha_options[1] 
  recomb_sel<-recomb_options[1] 
  Nstar_sel<-Nstar_options[2]
  
  testGrid_select<-subset(testGrid, 
                          initFreq_options==init_sel & 
                            alpha_options==alpha_sel & 
                            recomb_options==recomb_sel & 
                            Nstar_options==Nstar_sel)
  
  panel1<- ggplot(testGrid_select, aes(x=factor(h2_options), y=EY)) +
    geom_point(aes(colour=factor(Rmax_options), shape = factor(h2_options)), size = 4) +
    geom_line(aes(colour= factor(Rmax_options), group = factor(Rmax_options))) +
    scale_shape_manual(values = c(16, 17, 15)) +
    scale_color_manual(values = pal) +
    ylim(0, 1.01) +
    labs(title = "(a)", y= expression(Expected~Return~(E[Y])), x=expression(Heritability~(h^2)), color = expression(atop("Maximum reproductive", paste("rate"))), shape=expression(Heritability~(h^2))) +
    theme(panel.background = element_rect(fill = "white"), 
          panel.border = element_rect(fill=NA), 
          axis.line = element_blank())
  
  grob3<- panel1 + facet_wrap(~ m_options, nrow = 1) + 
    theme(legend.position = "right",
          legend.box = "vertical",
          legend.box.background = element_blank(),
          strip.background =element_rect(colour = "black", fill="lightgrey"),
          legend.key = element_rect(colour = "white", fill = "white", size = 0.25)) 
  
  panel2<- ggplot(testGrid_select, aes(x=(props_percent), y=(k_percent),)) +
    geom_point(aes(colour=factor(Rmax_options), shape = factor(h2_options)), size = 4) +
    xlim(-0.1, 0.51)+
    scale_shape_manual(values = c(1, 2, 22)) +
    scale_color_manual(values = pal) +
    ylim(0, 51) +
    labs(title = "(b)", y=expression(atop("Timing of TGf action in relation", paste("to environmental shift (years)"))), x="Proportion of pre-adapted individuals introduced", shape=expression(Heritability~(h^2)), color = expression(atop("Maximum reproductive", paste("rate")))) +
    theme(panel.background = element_rect(fill = "white"), 
          panel.border = element_rect(fill=NA), 
          axis.line = element_blank())
  
  grob4<- panel2 + facet_wrap(~ m_options, nrow = 1) + 
    theme(legend.position = "right",
          legend.box = "vertical",
          legend.box.background = element_blank(),
          strip.background =element_rect(colour = "black", fill="lightgrey"),
          legend.key = element_rect(colour = "white", fill = "white", size = 0.25)) 
  
  grob5<- arrangeGrob(grob3, grob4, ncol =1, nrow=2)
  ggsave(file="man/images/a4_figure7.pdf", grob5, width = plotWidth, height = plotHeight, units = "mm")
}

{print("MANUSCRIPT PLOTS FINISHED")}

#----------------------- Supplementary Figures  -------------------------------

#----------------------- Supp 1. shape of m  -------------------------------
x<-seq(from = -50, to = 50-1, by = 2)
initSD=1
m=0.15
plot((initSD*2)/(1+exp(-m*(x-(x/2))))+0.00,
     xlab = "Generation",
     ylab= "Environmental trait optimum")

m=0.3
lines((initSD*2)/(1+exp(-m*(x-(x/2))))+0.00)

m=1
lines((initSD*2)/(1+exp(-m*(x-(x/2))))+0.00, lty="dashed")

lines(seq(from = 0, to = initSD*2, by = initSD*2/26.75), lty="dotted", col="red", lwd=2)

#----------- Supp 2. comparison of managment OB  -------------------------------
{
  tile1<- ggplot((TGFstack[[303]]), aes(x=props, y=k, scale=props)) +
    geom_tile(aes(fill =(1-extinctRisk)*simpsonsDI)) +
    ylim(1, 51) +
    scale_fill_viridis(option = "A", begin = 0 , end = 1, na.value = "black", limits = c(0, 1), direction = 1) +
    theme_classic() +
    stat_contour(breaks= c(0.5, 0.9), aes(z=(1-extinctRisk)*simpsonsDI), color = "white", size = 0.6) +
    labs(title = "(a)", fill= expression(Expected~Return~(E[Y])), x="Proportion of pre-adapted individuals introduced", y=expression(atop("Timing of TGf action (years) in relation", paste("to environmental shift")))) +
    theme(legend.position=c("bottom"), legend.title = element_text(size=10), legend.text = element_text(size=7))
  
  tile2<- ggplot((TGFstack[[303]]), aes(x=props, y=k, scale=props)) +
    geom_tile(aes(fill =(1-extinctRisk)*r.pop)) +
    ylim(1, 51) +
    scale_fill_viridis(option = "A", begin = 0 , end = 1, na.value = "black", limits = c(0, 1), direction = 1) +
    theme_classic() +
    stat_contour(breaks= c(0.5, 0.9), aes(z=(1-extinctRisk)*r.pop), color = "white", size = 0.6)+
    labs(title = "(b)", fill="Expected Return", x="Proportion of pre-adapted individuals introduced", y="Timing of TGF action")+
    theme(legend.position=c("bottom"), legend.title = element_text(size=10), legend.text = element_text(size=7))
  
grob<- arrangeGrob(tile1, tile2, ncol=1)
ggsave(file="man/images/a4_S2.pdf", grob, width = 110, height = 297, units = "mm")
}

#----------- Supp 3. total managemnt space (E(Y_GSI))  -------------------------------

for (i in 1:length(TGFstack)){
  
  #create parameter specific title string  
  title<-paste0(testGrid[i,], collapse = " | ")
  
  plot<- ggplot((TGFstack[[i]]), aes(x=props, y=k, scale=props)) +
    geom_tile(aes(fill =extinctRisk) , colour = "white") +
    ylim(1, 51) +
    scale_fill_viridis(option = "A", begin = 0 , end = 1, na.value = "black", limits = c(0, 1)) +
    theme_classic() +
    labs(title= title, fill= "Probability of extinction (x)", x="Proportion of pre-adapted individuals introduced", y="Timing of TGf action") +
    theme(legend.position=c("bottom"), legend.title = element_text(size=10), legend.text = element_text(size=7))
  
  ggsave(paste0("extinctPlot",i,".pdf"), plot, width = 10, height = 8, dpi = 72, path ="out/plot/extinct" )
  
  #------------------------------------------------------------------------------------
  
  plot<- ggplot((TGFstack[[i]]), aes(x=props, y=k, scale=props)) +
    geom_tile(aes(fill =r.pop), colour = "white") +
    ylim(1, 51) +
    scale_fill_viridis(option = "D", begin = 0.5 , end = 1, na.value = "black", limits = c(0, 1)) +
    theme_classic()+
    labs(title= title,fill= "Proportion of recipient genome intact", x="Proportion of pre-adapted individuals introduced", y="Timing of TGF action")+
    theme(legend.position=c("bottom"), legend.title = element_text(size=10), legend.text = element_text(size=7))
  
  ggsave(paste0("r.popPlot",i,".png"), plot, width = 10, height = 8, dpi = 72, path ="out/plot/r.pop")
  
  #------------------------------------------------------------------------------------
  
  plot<- ggplot((TGFstack[[i]]), aes(x=props, y=k, scale=props)) +
    geom_tile(aes(fill =(1-extinctRisk)*r.pop), colour = "white") +
    ylim(1, 51) +
    scale_fill_viridis(option = "D", begin = 0 , end = 1, na.value = "black", limits = c(0, 1)) +
    theme_classic() +
    stat_contour(breaks= c(0.5, 0.9), aes(z=(1-extinctRisk)*r.pop), color = "white", size = 0.6)+
    labs(title = title, fill="Expected Return", x="Proportion of pre-adapted individuals introduced", y="Timing of TGF action")+
    theme(legend.position=c("bottom"), legend.title = element_text(size=10), legend.text = element_text(size=7))
  
  ggsave(paste0("expectPlot",i,".png"), plot, width = 10, height = 8, dpi = 72, path ="out/plot/expect")
  
  #------------------------------------------------------------------------------------
  
  plot<- ggplot((TGFstack[[i]]), aes(x=props, y=k, scale=props)) +
    geom_tile(aes(fill =simpsonsDI), colour = "white") +
    ylim(1, 51) +
    scale_fill_viridis(option = "D", begin = 0 , end = 1, na.value = "black", limits = c(0, 1)) +
    theme_classic()+
    labs(title= title,fill= "Gini-Simpson Index", x="Proportion of pre-adapted individuals introduced", y="Timing of TGF action") +
    theme(legend.position=c("bottom"), legend.title = element_text(size=10), legend.text = element_text(size=7))
  
  ggsave(paste0("gini-simpsonPlot",i,".png"), plot, width = 10, height = 8, dpi = 72, path ="out/plot/gini-simpson")
  
  #------------------------------------------------------------------------------------
  
  plot<- ggplot((TGFstack[[i]]), aes(x=props, y=k, scale=props)) +
    geom_tile(aes(fill = meanPheno), colour = "white") +
    ylim(1, 51) +
    scale_fill_viridis(option = "D", limits = c(-1, 3)) +
    theme_classic() +
    labs(title= title,fill= "Mean Phenotype", x="Proportion of pre-adapted individuals introduced", y="Timing of TGF action") +
    theme(legend.position=c("bottom"), legend.title = element_text(size=10), legend.text = element_text(size=7))
  
  ggsave(paste0("phenoPlot",i,".png"), plot, width = 10, height = 8, dpi = 72, path ="out/plot/pheno")
  
  #------------------------------------------------------------------------------------
  
  plot<- ggplot((TGFstack[[i]]), aes(x=props, y=k, scale=props)) +
    geom_tile(aes(fill=(1-extinctRisk)*simpsonsDI), colour = "white") +
    ylim(1, 51) +
    scale_fill_viridis(option = "D", begin = 0 , end = 1, na.value = "black", limits = c(0, 1)) +
    theme_classic() +
    stat_contour(breaks= c(0.5, 0.9), aes(z=(1-extinctRisk)*simpsonsDI), color = "white", size = 0.6) +
    labs(title = title, fill="Expected Return (Gini-Simpson)", x="Proportion of pre-adapted individuals introduced", y="Timing of TGF action")+
    theme(legend.position=c("bottom"), legend.title = element_text(size=10), legend.text = element_text(size=7))
  
  ggsave(paste0("expect_GSPlot",i,".png"), plot, width = 10, height = 8, dpi = 72, path ="out/plot/expect_GS")
  
  print(i)
}
