library(AccelStab)
library(ggplot2)

## Data ----
data(antigenicity)
plot(antigenicity$time, antigenicity$conc, col = antigenicity$K)

antigenicity$Validate <- as.factor(ifelse(antigenicity$time <= 0.5, 0, 1))

head(antigenicity)
table(antigenicity$Celsius)
unique(antigenicity$time)

## AccelStab fit ----
res <- step1_down(
  data = antigenicity, y = "conc",
  .time = "time", C = "Celsius",
  validation = "Validate",
  parms = list(k1 = 50, k2 = 10000, k3 = 3, c0 = 100),
  max_time_pred = 3
)
summary(res$fit)
confint(res$fit)


img_path <- "images/"

pdf(paste0(img_path, "accel_stab_fit.pdf"), width = 10, height = 6.5)
step1_plot_pred(res, yname = "Antigenicity")
dev.off()

# WO starting point
## AccelStab fit ----
res_wo <- step1_down(
  data = antigenicity, y = "conc",
  .time = "time", C = "Celsius",
  validation = "Validate",
  max_time_pred = 3
)
summary(res_wo$fit)
confint(res_wo$fit)


step1_plot_pred(res_wo, yname = "Antigenicity")
