#######################################################################
##*  Fig.S10-S15 Application case of United States maps
##* 1. read the raw data
##* 2. run BayesLASA
##* 3. plot identified landmarks based on BayesLASA output
#######################################################################
library(latex2exp)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(mcclust)
#***************
#* file path
#***************
input <- "manuscript_reproducibility/results/real_data_us_state/"
input_raw <- "data/real_data_us_state/"
fig.output <- "manuscript_reproducibility/figures_and_tables/"
code_path <- "code/"

#* raw data
load(file.path(input_raw, "states_outline.Rdata"), verbose = T)
state_list <- names(outline_polygon)
state_pns <- unlist(lapply(outline_polygon, nrow))

# Identify landmarks using BayesLASA --------------------------------------

####################################################################################
##** BayesLASA was on using the Rscript 'manuscript_reproducibility/results/real_data_us_state/run_BayesLASA.R'
##** Results are already saved in the 'manuscript_reproducibility/results/real_data_us_state/BayesLASA'
##** Landmark points and credible interval were identified using below code section
##** ready-to-use results saved in 'manuscript_reproducibility/results/real_data_us_state'
##** uncomment the below code section to re-run the codes
####################################################################################

# Get credible interval --------------------------------------
# df_cred_ppm_list = list()
# df_cred_ppi_list = list()
# for(mm in 1:length(state_list)){
#   state_name = state_list[mm]
#   original_dat = outline_polygon[[mm]]
#   normalized = pc_normalizor(original_dat)
#
#   dat = normalized$pc
#   dat = dat[-nrow(dat),]
#   n = nrow(dat)
#   beta_sigma = c(0.01, 0.001, 0.0001, 0.00001)
#   names(beta_sigma) = c("0.01","0.001",  "0.0001","0.00001")
#
#   for(k in 1:4){
#     r = readRDS(file = file.path(input, "BayesLASA",
#                                  paste0("state_", mm,"_bs_", names(beta_sigma)[k],"_MCMCres.rds")))
#     ##ppm
#     ppms = r$ppm
#     z_ppm <- minbinder(r$ppm, method = "comp")$cl
#     L_ppm = which(diff(c(z_ppm[n], z_ppm)) != 0)
#     ##ppi
#     ppis = colSums(r$PPI)/(r$iter-r$burn)
#     L_ppi = which(r$gamma_map > 0)
#
#     df_cred_ppm = data.frame(L_ppm)
#     df_cred_ppm = df_cred_ppm %>%
#       mutate(Method = "PPM", bs = names(beta_sigma)[k], state = state_name, state_ind = mm)
#     df_cred_ppi = data.frame(L_ppi)
#     df_cred_ppi = df_cred_ppi %>%
#       mutate(Method = "MAP", bs = names(beta_sigma)[k], state = state_name, state_ind = mm)
#     df_cred_ppm_list = c(df_cred_ppm_list,
#                          list(df_cred_ppm))
#     df_cred_ppi_list = c(df_cred_ppi_list,
#                          list(df_cred_ppi))
#     rm(r)
#   }
# }
#
# do.call('rbind', df_cred_ppi_list) %>%
#   write.csv(file.path(input, "map_df_cred_ppi.csv"))
#
# do.call('rbind', df_cred_ppm_list) %>%
#   write.csv(file.path(input, "map_df_cred_ppm.csv"))

# Load data and plot --------------------------------------

df_cred_ppi <- read.csv(file.path(input, "map_df_cred_ppi.csv"),
  row.names = 1
) %>% rename(L = L_ppi)
df_cred_ppm <- read.csv(file.path(input, "map_df_cred_ppm.csv"),
  row.names = 1
) %>% rename(L = L_ppm)

landmarks <- rbind(df_cred_ppi, df_cred_ppm)

#* list store landmarks
landmarks_list <- list()
#* list store credible interval
landmarks_cred_list <- list()
pcs <- list()
for (ind in unique(landmarks$state_ind)) {
  landmarks2 <- landmarks %>% filter(state_ind == ind)

  normalized <- pc_normalizor(outline_polygon[[ind]])
  dat <- normalized$pc
  centroid <- st_centroid(st_polygon(list(dat)))

  dat_c <- dat - rep(centroid, each = nrow(dat))

  pc <- as.data.frame(dat_c) %>%
    rename(x = long, y = lat) %>%
    mutate(state = unique(landmarks2$state), state_ind = unique(landmarks2$state_ind))
  pns <- nrow(outline_polygon[[ind]]) - 1
  pcs <- c(pcs, list(cbind(pc)))
  landmarks_list <- c(
    landmarks_list,
    list(cbind.data.frame(pc[landmarks2$L, c(1, 2)], cloc = "Landmarks", loc = "Landmarks", landmarks2))
  )
}
pcs <- do.call("rbind", pcs)
landmarks_list <- do.call("rbind", landmarks_list)
rownames(landmarks_list) <- NULL
landmarks_cred_list <- do.call("rbind", landmarks_cred_list)
landmarks_plot <- rbind(landmarks_list, landmarks_cred_list)

landmarks_plot <- landmarks_list
rownames(landmarks_plot) <- NULL

# split states to different plot
stat_ind_list <- split(1:48, ceiling(1:48 / 8))

#* change label
landmarks_plot <- landmarks_plot %>%
  mutate(
    bs = factor(bs,
      levels = c("0.01", "0.001", "0.0001", "0.00001"),
      labels = c(
        expression(b[sigma] == "0.01"),
        expression(b[sigma] == "0.001"),
        expression(b[sigma] == "0.0001"),
        expression(b[sigma] == "0.00001")
      )
    ),
    state = factor(state)
  )

# list of state names
landmarks_plot <- landmarks_plot %>%
  mutate(
    state = stringr::str_to_title(state),
    state_pn = paste0("(n = ", state_pns[state_ind], ")"),
    state_label = paste0(state, "\n", state_pn)
  ) %>%
  group_by(state, bs, Method) %>%
  mutate(landmark_n = n(), landmark_n = paste0("K = ", landmark_n)) # number of identified landmarks
#* state outline
pcs <- pcs %>%
  mutate(
    state = stringr::str_to_title(state),
    state_pn = paste0("(n = ", state_pns[state_ind], ")"),
    state_label = paste0(state, "\n", state_pn)
  )

pdf(file.path(fig.output, paste0("figure_s10-s15.pdf")), height = 14, width = 10)
for (i in 1:length(stat_ind_list)) {
  ind <- stat_ind_list[[i]]
  landmarks_plot2 <- landmarks_plot %>% filter(state_ind %in% ind, Method == "MAP")
  landmark_legend <- subset(pcs, state_ind %in% ind) %>%
    ungroup() %>%
    group_by(state) %>%
    mutate(legend_x = max(x), legend_y = max(y)) %>%
    slice(1)

  landmark_legend <- landmarks_plot2 %>%
    group_by(state, bs) %>%
    slice(1) %>%
    dplyr::select(bs, state, landmark_n) %>%
    left_join(landmark_legend, by = c("state"))

  p_map <- ggscatter(
    data = subset(landmarks_plot2, cloc == "Landmarks"),
    x = "x", y = "y",
    group = "state",
    color = "#FC4E07",
    size = 1.5, shape = 15
  ) +
    geom_point(
      data = subset(pcs, state_ind %in% ind),
      aes(x = x, y = y, group = state), color = "grey20",
      size = 0.6
    ) +
    geom_polygon(
      data = subset(pcs, state_ind %in% ind),
      aes(x = x, y = y, group = state), fill = NA, linetype = "solid",
      size = 0.6, color = "grey20"
    ) +
    geom_point(
      data = subset(landmarks_plot2, cloc == "Landmarks"),
      aes(x = x, y = y),
      colour = "#FC4E07", size = 1.5, shape = 15
    ) +
    geom_polygon(
      data = subset(landmarks_plot2, cloc == "Landmarks"),
      aes(x = x, y = y),
      colour = "#FC4E07", fill = NA,
      size = 0.7, linetype = "longdash", alpha = 0.2
    ) +
    geom_text(
      data = landmark_legend,
      aes(x = 0, y = legend_y + 0.02, label = landmark_n)
    ) +
    theme_void() +
    theme(
      strip.text = element_text(size = 12, colour = "grey30"),
      strip.background = element_rect(
        colour = "grey",
        fill = "grey"
      ),
      strip.text.y = element_text(angle = 90, colour = "grey30"),
      panel.border = element_rect(fill = NA)
    ) +
    # theme(aspect.ratio=1)
    coord_quickmap()
  p_maps <- p_map +
    facet_grid(state_label ~ bs, scales = "free", labeller = labeller(bs = label_parsed))
  print(p_maps)
}
dev.off()
