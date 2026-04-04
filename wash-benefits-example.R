


library(devtools)
load_all('/home/salbalkus/Research/RieszCML/')

library(nadir)
library(readr)
library(tidyverse)
# read in data via readr::read_csv
wash <- read_csv(
  paste0(
    "https://raw.githubusercontent.com/tlverse/tlverse-data/master/",
    "wash-benefits/washb_data.csv"
  )
)
normalize <- function(v){(v - min(v)) / (max(v) - min(v))}

df = data.frame(Y = wash$whz,
                A1 = wash$momage,
                L1 = wash$Nlt18,
                A2 = ifelse(wash$tr == "Control", 1, 0),
                L2 = wash$momedu) %>% drop_na()

#---------------------------------------------
# Stage 1 RieszCurve:
#   parameter is E[Y | L1, A1] under intervention A2 = 1
#---------------------------------------------
rc1 <- RieszCurve$new(
  nuis = function(data) {
    # regression for Y given current history
    m_fit <- glm(Y ~ L1 + A1 + L2 + A2,
                 family = gaussian(),
                 data = data)

    data_A2eq1 <- data
    data_A2eq1$A2 <- 1

    # treatment mechanism for A2
    g_fit <- glm(A2 ~ L1 + A1 + L2,
                 family = binomial(),
                 data = data)

    list(
      m = predict(m_fit, newdata = data, type = "response"),
      m1 = predict(m_fit, newdata = data_A2eq1, type = "response"),
      g = predict(g_fit, newdata = data, type = "response")
    )
  },
  h = ~ m1,
  f = ~ m,
  alpha = ~ A2 / g,
  ic_expr = ~ h + alpha * (Y - f)
)

#---------------------------------------------
# Stage 2 RieszCurve:
#   parameter is E[h[j=1] | L1] under intervention A1 = 1
#---------------------------------------------
rc2 <- function(delta){return(RieszCurve$new(
  nuis = function(data) {
    m_fit <- glm(`h[j=1]` ~ L1 + A1,
                 family = gaussian(),
                 data = data)

    data_shift <- data
    data_shift$A1 <- data_shift$A1 + delta

    # Fit the probabilistic classifier from Diaz et al. (2021)
    data_aug <- rbind(cbind(data, I = 0), cbind(data_shift, I = 1))
    r_fit <- glm(I ~ L1 + A1,
                 family = binomial(),
                 data = data_aug)

    list(
      m = predict(m_fit, newdata = data, type = "response"),
      m1 = predict(m_fit, newdata = data_shift, type = "response"),
      r = predict(r_fit, newdata = data_shift, type = "response")
    )
  },
  h = ~ m1,
  f = ~ m,
  alpha = ~ r / (1 - r),
  ic_expr = ~ h + alpha * (`h[j=1]` - f)
))}

#---------------------------------------------
# Compose the two stages
#---------------------------------------------
rc_composed <- function(delta){return(ComposedRieszCurve$new(
  rc_list = list(rc1, rc2(delta))
))}

#---------------------------------------------
# One-step estimate using composed object
#---------------------------------------------
est_onestep <- riesz_estimate(
  data = df,
  rc = rc_composed(3)
)

est_onestep

#---------------------------------------------
# TMLE using the composed object
#---------------------------------------------
est_tmle <- riesz_tmle(
  data = df,
  rc = rc_composed(3),
  fluctuation = "logistic",
  bounds = c(min(df$Y), max(df$Y)),
  outcome_col = "Y"
)

est_tmle
