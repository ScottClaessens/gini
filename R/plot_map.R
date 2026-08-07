# library(targets)
# library(rnaturalearth)
# library(gganimate)
# 
# tar_load(data)
# draws <- tar_read(fit_draws_model)
# 
# 
# # calculate average posterior trajectories
# pred <- tibble()
# for (r in 1:86) {
#   for (t in 1:121) {
#     pred <-
#       bind_rows(
#         pred,
#         tibble(
#           region = r,
#           date = t,
#           pop_size = 
#             median(
#               plogis(draws[[paste0("regional_latent_rep[", r, ",", t, ",1]")]])
#               * exp(draws[[paste0("regional_latent_rep[", r, ",", t, ",2]")]])
#             ),
#           cropland =
#             median(
#               plogis(draws[[paste0("regional_latent_rep[", r, ",", t, ",3]")]])
#               * exp(draws[[paste0("regional_latent_rep[", r, ",", t, ",4]")]])
#             ),
#           gini = 
#             median(
#               plogis(draws[[paste0("regional_latent_rep[", r, ",", t, ",5]")]])
#             )
#         )
#       )
#   }
# }
# 
# # edit predictions dataset
# pred2 <-
#   pred |>
#   rowwise() |>
#   mutate(
#     region = unique(data$subregion)[region],
#     longitude = unique(data$longitude[data$subregion == region]),
#     latitude = unique(data$latitude[data$subregion == region]),
#     date = seq(-100, 20)[date],
#     pop_size = log(pop_size + 1),
#     cropland = log(cropland + 1)
#   ) |>
#   filter(date <= 16)
# 
# world <-
#   ne_countries(
#     scale = "small",
#     returnclass = "sf"
#   ) |>
#   filter(continent != "Antarctica")
# 
# p <-
#   pred2 |>
#   ggplot() +
#   geom_sf(
#     data = world,
#     fill = "grey80",
#     colour = NA
#   ) +
#   geom_point(
#     mapping = aes(
#       x = longitude,
#       y = latitude,
#       colour = gini,
#       size = cropland
#     )
#   ) +
#   scale_colour_viridis_c() +
#   scale_size(range = c(0.1, 4)) +
#   theme_void()
# 
# animation <- 
#   p +
#   transition_time(date) +
#   labs(subtitle = "Date: {frame_time}") +
#   shadow_wake(wake_length = 0)
# 
# animate(
#   animation,
#   height = 500,
#   width = 800,
#   fps = 100,
#   duration = 10,
#   end_pause = 200,
#   res = 200
# )
# 
# anim_save("test.gif")
