#Probability of activity of mammals in composition, configuration and quality of native vegetation.
#10 different species of medium-large- size mammals in brazilian cerrado
#30/09/2026
#Sofia
rm(list=ls())

install.packages("TMB")
install.packages("glmmTMB")

library(GLMMadaptive)
library(readxl)
library(glmmTMB)
library(ggplot2)
source("https://raw.githubusercontent.com/MarcusRowcliffe/make_chm_data/refs/heads/main/chm_helpers.R")
library(dplyr)
library(lubridate)
library(openxlsx)

# Load data
obs <- read.xlsx("Data/all_species_data_solartime.xlsx") %>%
  droplevels() %>% 
  mutate(
    Date = as.Date(Date, origin = "1899-12-30"),
    timestamp = ymd_hm(
      paste(format(Date, "%Y-%m-%d"), substr(Solartime, 1, 5))
    ),
    deploymentID = LocationName
  )

# deployments
dep <- read.xlsx("data/cameratrapdata_native_veg.xlsx") %>%
  mutate(
    # Datas do Excel
    Date_setup = as.Date(data_instalacao, origin = "1899-12-30"),
    Date_retr  = as.Date(data_retirada, origin = "1899-12-30"),
    
    # Junta data + hora exatamente como registradas
    deploymentStart = as.POSIXct(
      Date_setup,
      tz = "UTC"
    ) + hora_instalacao * 86400,
    
    deploymentEnd = as.POSIXct(
      Date_retr,
      tz = "UTC"
    ) + hora_retirada * 86400,
    
    # Estrutura espacial
    deploymentID = paste0(toupper(ID), camera, "R"),
    locationName = deploymentID,
    cluster = ID
  )

View(
  dep %>%
    select(
      ID,
      camera,
      deploymentStart,
      deploymentEnd,
      deploymentID,
      locationName,
      cluster
    )
)


# ============================================================
# CHM DATA - Veado-catingueiro
# Melhores escalas selecionadas:
# Native: 1680 m
# ENN: 2120 m
# NDVI: 1680 m
# ============================================================

dat_veado_comp <- make_chm_data(
  dep,
  subset(obs, Species == "Veado-catingueiro"),
  covs = c(
    "locationName",
    "cluster",
    
    "native_1680",
    "enn_mn_2120",
    "ndvi_mean_1680"
  )
)

View(dat_veado_comp)

dat_veado_comp %>%
  select(
    locationName,
    cluster,
    timeHour,
    timeRadian,
    success,
    failure,
    native_1680,
    enn_mn_2120,
    ndvi_mean_1680
  ) %>%
  head(20)

# Número de registros independentes do Veado-catingueiro
nrow(subset(obs, Species == "Veado-catingueiro"))

sum(dat_veado_comp$success)


# Fit basic models...
# ...using mixed_model
uniform <- mixed_model(fixed = cbind(success, failure) ~ 1,
                       random = ~ 1 | locationName,
                       family = binomial(),
                       data = dat_veado_comp)

unimodal <- mixed_model(fixed = cbind(success, failure) ~ 
                          cos(timeRadian) + sin(timeRadian),
                        random = ~ 1 | locationName,
                        family = binomial(),
                        data = dat_veado_comp)

bimodal <- mixed_model(fixed = cbind(success, failure) ~ 
                         cos(timeRadian) + sin(timeRadian) +
                         cos(2*timeRadian) + sin(2*timeRadian), 
                       random = ~ 1 | locationName,
                       family = binomial(),
                       data = dat_veado_comp)

# AIC comparison
AIC(uniform, unimodal, bimodal)

# Generate predictions 
newdat <- data.frame(timeRadian = seq(0, 24, len=100) * pi / 12)

predict_uniform <- effectPlotData(uniform, newdat, marginal = FALSE) %>% 
  mutate(Model = "Uniform")
predict_unimodal <- effectPlotData(unimodal, newdat, marginal = FALSE) %>% 
  mutate(Model = "Unimodal")
predict_bimodal <- effectPlotData(bimodal, newdat, marginal = FALSE) %>% 
  mutate(Model = "Bimodal")

# join and plot results
rbind(predict_unimodal, predict_bimodal, predict_uniform) %>% 
  ggplot(., aes(x = timeRadian, y = plogis(pred), group = Model, fill = Model)) +
  geom_line(aes(colour = Model)) +
  geom_ribbon(aes(ymin = plogis(low), ymax = plogis(upp), colour = NULL), alpha = 0.3) +
  labs(x = "Time of Day (Hour)", y = "Predicted Activity Pattern \n (probability)") +
  theme_minimal() +
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=10, face="bold"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,-10,-10,-10),
        plot.title = element_blank(),
        axis.line = element_line(colour = 'black', linetype = 'solid'),
        axis.ticks = element_line(colour = 'black', linetype = 'solid'),
        axis.title = element_text(size=9, face="bold"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = 'lightgrey', linetype = 'dashed', linewidth=0.5),
        panel.grid.minor.x = element_blank(),
        strip.text = element_text(size = 9, colour = "black", face = "bold", hjust = 0)
  ) +
  scale_x_continuous(breaks=seq(0, 2*pi, len=5), labels=seq(0, 24, len=5))


# ============================================================
# COMPARISON OF LANDSCAPE VARIABLES (BEST SCALES)
# Veado-catingueiro
# ============================================================

# Null model
mod_null <- chmTMB(
  cbind(success, failure) ~
    time +
    (1 | cluster) +
    (1 | locationName),
  data = dat_veado_comp
)

# Native vegetation (1680 m scale)
mod_native_1680 <- chmTMB(
  cbind(success, failure) ~
    time * native_1680 +
    (1 | cluster) +
    (1 | locationName),
  data = dat_veado_comp
)

# Native vegetation (1680 m scale) and fragments (2120 m scale)
mod_native_1680_fragments_2120 <- chmTMB(
  cbind(success, failure) ~
    time * native_1680 + enn_mn_2120 +
    (1 | cluster) +
    (1 | locationName),
  data = dat_veado_comp
)

# Native vegetation (1680 m scale) and ndvi (1680 m scale)
mod_native_1680_ndvi_1680 <- chmTMB(
  cbind(success, failure) ~
    time * native_1680 + ndvi_mean_1680 +
  (1 | cluster) +
    (1 | locationName),
  data = dat_veado_comp
)
# Distance between fragments (2120 m scale)
mod_fragments_2120 <- chmTMB(
  cbind(success, failure) ~
    time * enn_mn_2120 +
    (1 | cluster) +
    (1 | locationName),
  data = dat_veado_comp
)

# NDVI (1680 m scale)
mod_ndvi_1680 <- chmTMB(
  cbind(success, failure) ~
    time * ndvi_mean_1680 +
    (1 | cluster) +
    (1 | locationName),
  data = dat_veado_comp
)

# Distance between fragments (2120 m scale) and ndvi (1680 m scale)
mod_fragments_2120_ndvi_1680 <- chmTMB(
  cbind(success, failure) ~
    time * enn_mn_2120 + ndvi_mean_1680 +
    (1 | cluster) +
    (1 | locationName),
  data = dat_veado_comp
)


# ============================================================
# CHECK MODEL CONVERGENCE
# ============================================================

c(
  null = mod_null$sdr$pdHess,
  native_1680 = mod_native_1680$sdr$pdHess,
  native_1680_fragments_2120 = mod_native_1680_fragments_2120$sdr$pdHess,
  native_1680_ndvi_1680 = mod_native_1680_ndvi_1680$sdr$pdHess,
  fragments_2120 = mod_fragments_2120$sdr$pdHess,
  ndvi_1680 = mod_ndvi_1680$sdr$pdHess,
  fragments_2120_ndvi_1680 = mod_fragments_2120_ndvi_1680$sdr$pdHess
)

# ============================================================
# MODEL SELECTION BASED ON AIC
# ============================================================

tab_aic <- AIC(
  mod_null,
  mod_native_1680,
  mod_native_1680_fragments_2120,
  mod_native_1680_ndvi_1680,
  mod_fragments_2120,
  mod_ndvi_1680,
  mod_fragments_2120_ndvi_1680
) %>%
  mutate(
    dAIC = AIC - min(AIC, na.rm = TRUE)
  ) %>%
  arrange(dAIC)

tab_aic


## ============================================================
# PREDICTED ACTIVITY PATTERN - CORRIGIDO E LIMPO
# Native vegetation (1680 m scale)
# ============================================================

predict(
  mod_native_1680,
  list(
    native_1680 = quantile(
      dat_veado_comp$native_1680,
      c(0.10, 0.90), # Usando 10% e 90% para evitar explosão nos extremos
      na.rm = TRUE
    )
  )
) %>%
  mutate(
    native_level = factor(
      ifelse(native_1680 == min(native_1680), "Low", "High"),
      levels = c("Low", "High") # Força a ordem correta: Low à esquerda, High à direita
    )
  ) %>%
  ggplot(
    aes(
      time,
      response,
      col = native_level,
      group = native_level
    )
  ) +
  geom_ribbon(
    aes(
      ymin = lcl.response,
      ymax = ucl.response,
      fill = native_level
    ),
    col = NA,
    alpha = 0.2
  ) +
  geom_line(linewidth = 0.8) +
  scale_x_continuous(
    breaks = seq(0, 24, len = 5)
  ) +
  facet_grid(~native_level) +
  theme_classic()


# Todos os números devem ser rigorosamente IGUAIS (ex: 3888)
sapply(list(Null = mod_null, Native = mod_native_1680, Fragments = mod_fragments_2120, NDVI = mod_ndvi_1680), nobs)
# Todos devem retornar TRUE
c(null = mod_null$sdr$pdHess, native = mod_native_1680$sdr$pdHess, enn = mod_fragments_2120$sdr$pdHess, ndvi = mod_ndvi_1680$sdr$pdHess)
summary(mod_native_1680)
