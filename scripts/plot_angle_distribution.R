library(ggplot2)
library(dplyr)
library(readr)

angle_data <- read_tsv("results/angles.tsv")

angle_data$size_class <- factor(angle_data$size_class,
                                levels = c("Tiny", "Small", "Intermediate", "Large", "Bulky"))

angle_data <- angle_data %>%
  mutate(angle = ((angle + 180) %% 360) - 180)

size_class_colors <- c(
  "Tiny"         = "#74c0e0",
  "Small"        = "#3a9fd1",
  "Intermediate" = "#1971b8",
  "Large"        = "#0c4a8c",
  "Bulky"        = "#08306b"
)

total_count <- nrow(angle_data)

angle_plot <- ggplot(angle_data, aes(x = angle, color = size_class)) +

  geom_density(size = 1.2, adjust = 1.2) +

  scale_color_manual(values = size_class_colors) +

  scale_x_continuous(
    limits = c(-180, 180),
    breaks = seq(-180, 180, by = 50)
  ) +

  labs(
    title = paste0("Tripeptide (XRX) in Helix (n = ", total_count, ")"),
    x = "Angle between adjacent C-alpha → Centroid vectors [°]",
    y = "Norm. Freq. [A.U.]",
    color = NULL
  ) +

  theme_minimal(base_size = 14) +

  theme(
    panel.background = element_rect(fill = "#bdbdbd", color = NA),
    plot.background  = element_rect(fill = "#bdbdbd", color = NA),

    panel.grid.major.x = element_line(color = "#dddddd", linetype = "dotted"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),

    legend.position   = c(0.2, 0.8),
    legend.background = element_rect(fill = "white", color = "black"),

    plot.title = element_text(hjust = 0.5)
  )

ggsave("results/angle_plot.png", plot = angle_plot, width = 10, height = 6, dpi = 300)

print(angle_plot)
