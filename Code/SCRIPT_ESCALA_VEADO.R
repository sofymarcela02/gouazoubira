#Scale of effect of composition (native vegetation),configuration ( enn_mn distance
#between fragments) and NDVi of native vegetation for subulo gouazobira

#10 different species of medium-large- size mammals in brazilian cerrado
#30/09/2026
#Sofia
#try

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
# CHM DATA - VEADO-CATINGUEIRO
# ============================================================

dat_veado_comp <- make_chm_data(
  dep,
  subset(obs, Species == "Veado-catingueiro"), 
  covs = c(
    "locationName",
    "cluster",
    
    # Distância entre fragmentos
    "enn_mn_800",
    "enn_mn_1240",
    "enn_mn_1680",
    "enn_mn_2120",
    "enn_mn_2560",
    "enn_mn_3000",
    
    # Vegetação Nativa (%)
    "native_800",
    "native_1240",
    "native_1680",
    "native_2120",
    "native_2560",
    "native_3000",
    
    # Média do NDVI
    "ndvi_mean_800",
    "ndvi_mean_1240",
    "ndvi_mean_1680",
    "ndvi_mean_2120",
    "ndvi_mean_2560",
    "ndvi_mean_3000"
  )
)




# Número de registros independentes do Veado-catingueiro
nrow(subset(obs, Species == "Veado-catingueiro"))

sum(dat_veado_comp$success)
# ============================================================
# MODELO NULO
# ============================================================

mod_null <- chmTMB(
  cbind(success, failure) ~
    time +
    (1 | cluster) +
    (1 | locationName),
  data = dat_veado_comp
)
# ============================================================
# EFFECT-SCALE MODELS - NATIVE VEGETATION (native)
# Veado-catingueiro
# ============================================================

# 800 m scale
mod_native_800 <- chmTMB(
  cbind(success, failure) ~
    time * native_800 +
    (1 | cluster) +
    (1 | locationName),
  data = dat_veado_comp
)


# 1240 m scale
mod_native_1240 <- chmTMB(
  cbind(success, failure) ~
    time * native_1240 +
    (1 | cluster) +
    (1 | locationName),
  data = dat_veado_comp
)


# 1680 m scale
mod_native_1680 <- chmTMB(
  cbind(success, failure) ~
    time * native_1680 +
    (1 | cluster) +
    (1 | locationName),
  data = dat_veado_comp
)


# 2120 m scale
mod_native_2120 <- chmTMB(
  cbind(success, failure) ~
    time * native_2120 +
    (1 | cluster) +
    (1 | locationName),
  data = dat_veado_comp
)


# 2560 m scale
mod_native_2560 <- chmTMB(
  cbind(success, failure) ~
    time * native_2560 +
    (1 | cluster) +
    (1 | locationName),
  data = dat_veado_comp
)


# 3000 m scale
mod_native_3000 <- chmTMB(
  cbind(success, failure) ~
    time * native_3000 +
    (1 | cluster) +
    (1 | locationName),
  data = dat_veado_comp
)


# ============================================================
# CHECK MODEL CONVERGENCE - NATIVE VEGETATION
# ============================================================

c(
  null = mod_null$sdr$pdHess,
  native_800 = mod_native_800$sdr$pdHess,
  native_1240 = mod_native_1240$sdr$pdHess,
  native_1680 = mod_native_1680$sdr$pdHess,
  native_2120 = mod_native_2120$sdr$pdHess,
  native_2560 = mod_native_2560$sdr$pdHess,
  native_3000 = mod_native_3000$sdr$pdHess
)


# ============================================================
# MODEL SELECTION BASED ON AIC - NATIVE VEGETATION
# ============================================================

tab_aic_native <- AIC(
  mod_null,
  mod_native_800,
  mod_native_1240,
  mod_native_1680,
  mod_native_2120,
  mod_native_2560,
  mod_native_3000
) %>%
  mutate(
    dAIC = AIC - min(AIC, na.rm = TRUE)
  ) %>%
  arrange(dAIC)

tab_aic_native


# ============================================================
# EFFECT-SCALE MODELS - NDVI (ndvi_mean)
# Veado-catingueiro
# ============================================================

# 800 m scale
mod_ndvi_800 <- chmTMB(
  cbind(success, failure) ~
    time * ndvi_mean_800 +
    (1 | cluster) +
    (1 | locationName),
  data = dat_veado_comp
)


# 1240 m scale
mod_ndvi_1240 <- chmTMB(
  cbind(success, failure) ~
    time * ndvi_mean_1240 +
    (1 | cluster) +
    (1 | locationName),
  data = dat_veado_comp
)


# 1680 m scale
mod_ndvi_1680 <- chmTMB(
  cbind(success, failure) ~
    time * ndvi_mean_1680 +
    (1 | cluster) +
    (1 | locationName),
  data = dat_veado_comp
)


# 2120 m scale
mod_ndvi_2120 <- chmTMB(
  cbind(success, failure) ~
    time * ndvi_mean_2120 +
    (1 | cluster) +
    (1 | locationName),
  data = dat_veado_comp
)


# 2560 m scale
mod_ndvi_2560 <- chmTMB(
  cbind(success, failure) ~
    time * ndvi_mean_2560 +
    (1 | cluster) +
    (1 | locationName),
  data = dat_veado_comp
)


# 3000 m scale
mod_ndvi_3000 <- chmTMB(
  cbind(success, failure) ~
    time * ndvi_mean_3000 +
    (1 | cluster) +
    (1 | locationName),
  data = dat_veado_comp
)


# ============================================================
# CHECK MODEL CONVERGENCE - NDVI
# ============================================================

c(
  null = mod_null$sdr$pdHess,
  ndvi_mean_800 = mod_ndvi_800$sdr$pdHess,
  ndvi_mean_1240 = mod_ndvi_1240$sdr$pdHess,
  ndvi_mean_1680 = mod_ndvi_1680$sdr$pdHess,
  ndvi_mean_2120 = mod_ndvi_2120$sdr$pdHess,
  ndvi_mean_2560 = mod_ndvi_2560$sdr$pdHess,
  ndvi_mean_3000 = mod_ndvi_3000$sdr$pdHess
)


# ============================================================
# MODEL SELECTION BASED ON AIC - NDVI
# ============================================================

tab_aic_ndvi <- AIC(
  mod_null,
  mod_ndvi_800,
  mod_ndvi_1240,
  mod_ndvi_1680,
  mod_ndvi_2120,
  mod_ndvi_2560,
  mod_ndvi_3000
) %>%
  mutate(
    dAIC = AIC - min(AIC, na.rm = TRUE)
  ) %>%
  arrange(dAIC)

tab_aic_ndvi



# ============================================================
# EFFECT-SCALE MODELS - DISTÂNCIA ENTRE FRAGMENTOS (enn_mn)
# Veado-catingueiro
# ============================================================

# 800 m scale
mod_fragments_800 <- chmTMB(
  cbind(success, failure) ~
    time * enn_mn_800 +
    (1 | cluster) +
    (1 | locationName),
  data = dat_veado_comp
)


# 1240 m scale
mod_fragments_1240 <- chmTMB(
  cbind(success, failure) ~
    time * enn_mn_1240 +
    (1 | cluster) +
    (1 | locationName),
  data = dat_veado_comp
)


# 1680 m scale
mod_fragments_1680 <- chmTMB(
  cbind(success, failure) ~
    time * enn_mn_1680 +
    (1 | cluster) +
    (1 | locationName),
  data = dat_veado_comp
)


# 2120 m scale
mod_fragments_2120 <- chmTMB(
  cbind(success, failure) ~
    time * enn_mn_2120 +
    (1 | cluster) +
    (1 | locationName),
  data = dat_veado_comp
)


# 2560 m scale
mod_fragments_2560 <- chmTMB(
  cbind(success, failure) ~
    time * enn_mn_2560 +
    (1 | cluster) +
    (1 | locationName),
  data = dat_veado_comp
)


# 3000 m scale
mod_fragments_3000 <- chmTMB(
  cbind(success, failure) ~
    time * enn_mn_3000 +
    (1 | cluster) +
    (1 | locationName),
  data = dat_veado_comp
)


# ============================================================
# CHECK MODEL CONVERGENCE - DISTÂNCIA ENTRE FRAGMENTOS
# ============================================================

c(
  null = mod_null$sdr$pdHess,
  enn_mn_800 = mod_fragments_800$sdr$pdHess,
  enn_mn_1240 = mod_fragments_1240$sdr$pdHess,
  enn_mn_1680 = mod_fragments_1680$sdr$pdHess,
  enn_mn_2120 = mod_fragments_2120$sdr$pdHess,
  enn_mn_2560 = mod_fragments_2560$sdr$pdHess,
  enn_mn_3000 = mod_fragments_3000$sdr$pdHess
)


# ============================================================
# MODEL SELECTION BASED ON AIC - DISTÂNCIA ENTRE FRAGMENTOS
# ============================================================

tab_aic <- AIC(
  mod_null,
  mod_fragments_800,
  mod_fragments_1240,
  mod_fragments_1680,
  mod_fragments_2120,
  mod_fragments_2560,
  mod_fragments_3000
) %>%
  mutate(
    dAIC = AIC - min(AIC, na.rm = TRUE)
  ) %>%
  arrange(dAIC)

tab_aic




