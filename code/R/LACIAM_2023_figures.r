library(tidyverse)
library(outbreaks)

data(influenza_england_1978_school)

boarding_data <- ggplot(data = influenza_england_1978_school,
       aes(x = date, y = in_bed)) + 
  geom_point(shape = 21, colour = "black", fill = "white", size = 3) + 
  scale_x_date("", expand = c(0, 0)) + 
  scale_y_continuous("Number of students in bed",
                     expand = c(0, 0), limits = c(0, 310)) +
  theme_bw(base_size = 16)

boarding_data
ggsave(plot = boarding_data,
       filename = "../../figures/boarding_school_data.pdf",
       scale = 1,
       width = 297,
       height = 210,
       units = "mm",
       dpi = 300)

#############
## pushforward

prior1 <- rnorm(1E6, sd = 1)
prior2 <- rnorm(1E6, sd = 5)

p1 <- arm::invlogit(x = prior1)
p2 <- arm::invlogit(x = prior2)

par(mfrow = c(1, 2))
hist(prior2, density = TRUE, main = "",
     xlab = expression(beta[1]), col = "grey50")
hist(prior1, density = TRUE, add = TRUE, col = 2)
legend(x = "topright",
       legend = c("Informative", "Non-informative"),
       col = c("red", "grey50"),
       pch = 16,
       bty = 'n')

hist(p2, density = TRUE, main = "",
     xlab = expression(p(x)), col = "grey50")
hist(p1, density = TRUE, add = TRUE, col = 2)



