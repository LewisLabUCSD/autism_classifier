#Fig 4 ============
data %>%
  ggbetweenstats(
    x     = AgeBins.DJx.Comb,
    y     = weightedEstimate,
    p.adjust.method  = "BH",
    palette = "Dark2",
    title = "Distribution of Classification Score",
    ggplot.component = list(theme(text = element_text(size = 15)
    )))


#Fig 5 ============
theme_set(theme_gray(base_size =15))
data %>%
ggplot( aes(y = EnsembleScore, 
                                        x = Train.Test, fill = Train.Test)) + 
  scale_fill_manual(name = "",values = c("orchid1","pink1")) + 
  geom_boxplot (alpha =  0.5, width = 0.7, outlier.shape = NA) + 
  geom_jitter(  aes( colour = Mutation.color, 
                     shape = SFARI.Level.1.OR.2 ), size=3, alpha=0.9, 
                position = position_jitterdodge(jitter.width = .3, dodge.width = 0.1))   +
  scale_color_manual(name = "", values = c("blue" = "blue", "grey" = "gray", "red" = "red", "orange" = "orange"),  
                     label = c("No Mutation", "No MIPs data", "Missense mutation", "LGD")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        legend.title=element_text(size=0))+
  facet_wrap(diagnosis_binary~.) 