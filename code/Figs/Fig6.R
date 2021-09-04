#***********************************************
#*Fig.6. example of lung cancer case study
#***********************************************

#***************
#*file path
#***************
code_file = "code/landmark_detection/"
fig.output = "figs"
input = "data/nlst"
raw_dat = file.path(input, "raw")
data.loc =  file.path(input, "processed")
rdatloc = file.path(input, "BayesLASA")


##good 10374

sg = 10374
f = paste0(sg, ".Rdata")
load(file.path(raw_dat, f), verbose = T)

pdf(file.path(fig.output, "Fig1_GoodPrognosis_rawimage.pdf"), width = 4, height = 4)
image(image, col = c(  "#230077","#2fbc83","#f7f70f"), xaxt = "n", yaxt = "n")
dev.off()

f = paste0(sg, "_outline.Rdata")
f2 = paste0(sg, "Sigma_bs500_normal_4chain_largest.rdata")
load(file.path(data.loc, f))
load(file.path(rdatloc, f2), verbose = T)
L = L_ppm[[1]]
pc  = outline_polygon[[1]]
colnames(pc) = c("x", "y")

p_g1 = ggplot(as.data.frame(pc), aes(x = x, y = y)) +
  geom_polygon(fill = NA, color= "black") +
  xlim(0, 1400) +
  ylim(0, 900)+
  theme_pubr()+
  ggtitle("Good prognosis") +
  theme(axis.title=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(colour = "black", size=0.8))

p_g2 = p_g1 +
  geom_point(data = as.data.frame(pc[L,]), aes(x = x, y = y), col = "red") +
  geom_polygon(data = as.data.frame(pc[L,]), aes(x = x, y = y), fill = NA, color= "red", linetype = "dashed") 

p_g_density = Roughness %>% filter(sample == sg) %>%
  ggplot(aes(x = Ra)) +
  geom_density(fill="lightgrey") +
  ggtitle("Good prognosis") +
  theme_pubr()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.line=element_blank(),
        panel.background = element_rect(colour = "black", size=0.8)) +
  xlim(c(0, 8))

ggsave(file.path(fig.output, "Fig6_GoodPrognosis_TumorRegion.pdf"), p_g1, width = 4, height = 4)
ggsave(file.path(fig.output, "Fig6_GoodPrognosis_TumorLandmak.pdf"), p_g2, width = 4, height = 4)
ggsave(file.path(fig.output, "Fig6_GoodPrognosis_density.pdf"), p_g_density, width = 4, height = 4)
##poor 10127
sp = 10127
f = paste0(sp, ".Rdata")
load(file.path(raw_dat, f), verbose = T)
pdf(file.path(fig.output, "Fig1_PoorPrognosis_rawimage.pdf"), width = 4, height = 4)
image(image, col = c(  "#230077","#2fbc83","#f7f70f"), xaxt = "n", yaxt = "n")
dev.off()
f = paste0(sp, "_outline.Rdata")
f2 = paste0(sp, "Sigma_bs500_normal_4chain_largest.rdata")
load(file.path(data.loc, f))
load(file.path(rdatloc, f2), verbose = T)
L = L_ppm[[1]]
pc  = outline_polygon[[1]]
colnames(pc) = c("x", "y")

p_p1 = ggplot(as.data.frame(pc), aes(x = x, y = y)) +
  geom_polygon(fill = NA, color= "black") +
  theme_pubr()+
  xlim(0, 1400) +
  ylim(0, 900) +
  ggtitle("Poor prognosis") +
  theme(axis.title=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(colour = "black", size=0.8))

p_p2 = p_p1 +
  geom_point(data = as.data.frame(pc[L,]), aes(x = x, y = y), col = "red") +
  geom_polygon(data = as.data.frame(pc[L,]), aes(x = x, y = y), fill = NA, color= "red", linetype = "dashed") 

p_p_density = Roughness %>% filter(sample == sp) %>%
  ggplot(aes(x = Ra)) +
  geom_density(fill = "lightgrey") +
  ggtitle("Poor prognosis") +
  theme_pubr()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.line=element_blank(),
        panel.background = element_rect(colour = "black", size=0.8)) +
  xlim(c(0, 8))

ggsave(file.path(fig.output, "Fig6_PoorPrognosis_TumorRegion.pdf"), p_p1, width = 4, height = 4)
ggsave(file.path(fig.output, "Fig6_PoorPrognosis_TumorLandmak.pdf"), p_p2, width = 4, height = 4)
ggsave(file.path(fig.output, "Fig6_PoorPrognosis_density.pdf"), p_p_density, width = 4, height = 4)