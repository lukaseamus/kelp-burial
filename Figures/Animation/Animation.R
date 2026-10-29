# 1. Proportion ####
# 1.1 Load data ####
sed <- here("Carbon", "RDS", "sed.rds") %>%
  read_rds() %T>%
  print()

images <- sed %>%
  distinct(Species, Sediment, Load) %>%
  mutate(
    image = case_when(
      Species == "Ecklonia radiata" & Sediment == "Coarse sand" ~
        here("Figures", "Illustrations", "Ecklonia.pdf"),
      Species == "Laminaria hyperborea" & Load == 1 ~
        here("Figures", "Illustrations", "Laminaria hyperborea.pdf"),
      Species == "Laminaria ochroleuca" & Load == 1 ~
        here("Figures", "Illustrations", "Laminaria ochroleuca.pdf")
    ),
    x = 8.5, y = .8
  ) %T>%
  print()

# 1.2 Animated figures ####
# Define theme
mytheme <- theme(
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  plot.margin = margin(0.2, 0.5, 0.2, 0.2, unit = "cm"),
  axis.line = element_line(),
  axis.title = element_text(size = 12, hjust = 0),
  axis.text = element_text(size = 10, colour = "black"),
  axis.ticks.length = unit(.25, "cm"),
  axis.ticks = element_line(colour = "black", lineend = "square"),
  legend.key = element_blank(),
  legend.key.width = unit(.25, "cm"),
  legend.key.height = unit(.45, "cm"),
  legend.key.spacing.x = unit(.5, "cm"),
  legend.key.spacing.y = unit(.05, "cm"),
  legend.background = element_blank(),
  legend.position = "top",
  legend.justification = 0,
  legend.text = element_text(size = 12, hjust = 0),
  legend.title = element_blank(),
  legend.margin = margin(0, 0, 0, 0, unit = "cm"),
  strip.background = element_blank(),
  strip.text = element_text(size = 12, hjust = 0),
  panel.spacing = unit(.8, "cm"),
  text = element_text(family = "Futura")
)

# For this simple animation task I can use gganimate directly
# rather than tweening manually
require(gganimate)
require(ggh4x)
require(ggimage)
( sed %>%
    filter(
      Depth < 10 & # filter for depths represented across all data
        # remove additional fine controls because they stack in geom_col()
        !ID_core %in% c("E_C_2_1", "E_C_2_2", "E_C_2_3") 
    ) %>%
    ggplot() +
      geom_col(aes(Depth, Proportion, fill = Species)) +
      # Since the images are just static, they can be overlaid afterwards
      # geom_image(
      #   data = images,
      #   aes(x, y, image = image),
      #   size = .3, by = "height"
      # ) +
      scale_fill_manual(values = c("#c3b300", "#627d0e", "#f1c700"),
                        guide = "none") +
      scale_x_reverse(breaks = 0:10) +
      scale_y_continuous(
        breaks = seq(0, 1, .5),
        labels = scales::label_number(accuracy = c(1, .1, 1))
      ) +
      facet_nested(
        ~ Species + Sediment + Load,
        nest_line = TRUE, solo_line = TRUE,
        scales = "free",
        strip = strip_nested(
          clip = "off", by_layer_x = TRUE,
          text_x = list(element_text(face = "italic"), NULL, NULL)
        ),
        labeller = labeller(
          Load = c(
            "1" = "1 g",
            "7.5" = "7.5 g",
            "10" = "10 g",
            "15" = "15 g"
          ),
          Sediment = c( # Abbreviate sand
            "Coarse sand" = "Coarse",
            "Fine sand" = "Fine"
          )
        )
      ) +
      labs(
        y = "Remaining proportion of kelp carbon after day {floor(frame_time +.5)}",
        x = "Sediment depth (cm)"
      ) +
      coord_flip(
        xlim = c(0, 10), ylim = c(0, 1), 
        expand = FALSE, clip = "off"
      ) +
      transition_time(Day) +
      mytheme ) %>%
  animate(
    fps = 25, duration = 20, width = 20, height = 10, units = "cm", 
    res = 300, renderer = gifski_renderer(), device = "ragg_png"
  ) %>%
  anim_save(
    filename = "prop.gif", 
    path = here("Figures", "Animation")
  )

# Alternative version on log scale
( sed %>%
    mutate(
      Proportion = if_else(
        Proportion == 0, # replace zeros with minima (presumed detection limit)
        min(Proportion[Proportion > 0]),
        Proportion
      )
    ) %>%
    filter(
      Depth < 10 & # filter for depths represented across all data
        # remove additional fine controls to avoid overplotting
        !ID_core %in% c("E_C_2_1", "E_C_2_2", "E_C_2_3") 
    ) %>%
    ggplot() +
      geom_rect(
        aes(xmin = Depth - 0.45, xmax = Depth + 0.45,
            ymin = 1e-6, ymax = Proportion, fill = Species)
      ) +
      scale_fill_manual(values = c("#c3b300", "#627d0e", "#f1c700"),
                        guide = "none") +
      scale_x_reverse(breaks = 0:10) +
      scale_y_log10(
        breaks = 10^seq(-6, 0, 3),
        labels = scales::label_log()
      ) +
      facet_nested(
        ~ Species + Sediment + Load,
        nest_line = TRUE, solo_line = TRUE,
        scales = "free",
        strip = strip_nested(
          clip = "off", by_layer_x = TRUE,
          text_x = list(element_text(face = "italic"), NULL, NULL)
        ),
        labeller = labeller(
          Load = c(
            "1" = "1 g",
            "7.5" = "7.5 g",
            "10" = "10 g",
            "15" = "15 g"
          ),
          Sediment = c( # Abbreviate sand
            "Coarse sand" = "Coarse",
            "Fine sand" = "Fine"
          )
        )
      ) +
      labs(
        y = "Remaining proportion of kelp carbon after day {floor(frame_time +.5)}",
        x = "Sediment depth (cm)"
      ) +
      coord_flip(
        xlim = c(0, 10), ylim = 10^c(-6, 0), 
        expand = FALSE, clip = "off"
      ) +
      transition_time(Day) +
      mytheme ) %>%
  animate(
    fps = 25, duration = 30, width = 20, height = 10, units = "cm", 
    res = 300, renderer = gifski_renderer(), device = "ragg_png"
  ) %>%
  anim_save(
    filename = "prop_log.gif", 
    path = here("Figures", "Animation")
  )