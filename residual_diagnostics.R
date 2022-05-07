library(tidyverse)

setwd("D:/Academic/Projects/2020_Covid 19 Prediction Model_Anirban Ghatak_IIMV/")
estList <- readRDS('./Datasets/estimate.Rds')

# C1, C2, C3 plots with 95\% CI
tibble("Date" = as.Date(names(estList$c1)), "C1" = estList$c1, 
       "Lower" = estList$c1 - 1.96 * sqrt(estList$varc1) / sqrt(1000), 
       "Upper" = pmin(estList$c1 + 1.96 * sqrt(estList$varc1) / sqrt(1000), 1)  ) %>%
    ggplot(aes(x = Date)) +
    geom_line(aes(y = C1), color = "darkgreen", size = 1.5) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "darkgreen", alpha = 0.25) +
    ylab(expression(c[1~t])) +
    ggtitle(expression(Time~Varying~Proportion~of~New~`In-Migrants`~Not~Infected~frac(dM[nt],dt)/M[t] )) +
    scale_x_date(breaks = scales::pretty_breaks(n = 5)) +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 12)) +
    ggsave("Chhattisgarh-mig-c1.pdf", device = "pdf", width = 12, height = 5)
    

tibble("Date" = as.Date(names(estList$c2)), "C2" = estList$c2, 
       "Lower" = pmax(estList$c2 - 1.96 * sqrt(estList$varc2) / sqrt(1000), 0), 
       "Upper" = pmin(estList$c2 + 1.96 * sqrt(estList$varc2) / sqrt(1000), 1)  ) %>%
    ggplot(aes(x = Date)) +
    geom_line(aes(y = C2), color = "orange", size = 1.5) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "orange", alpha = 0.25) +
    ylab(expression(c[2~t])) +
    ggtitle(expression(Time~Varying~Proportion~of~New~`In-Migrants`~Infected~and~Pre-Symptomatic~frac(dM[pt],dt)/M[t] )) +
    scale_x_date(breaks = scales::pretty_breaks(n = 5)) +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 12)) +
    ggsave("Chhattisgarh-mig-c2.pdf", device = "pdf", width = 12, height = 5)


tibble("Date" = as.Date(names(estList$c1)), "C3" = estList$c3, 
       "Lower" = pmax(0, estList$c3 - 1.96 * sqrt(estList$varc3) / sqrt(1000)), 
       "Upper" = pmin(estList$c3 + 1.96 * sqrt(estList$varc3) / sqrt(1000), 1)  ) %>%
    ggplot(aes(x = Date)) +
    geom_line(aes(y = C3), color = "red", size = 1.5) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "red", alpha = 0.25) +
    ylab(expression(c[3~t])) +
    ggtitle(expression(Time~Varying~Proportion~of~New~`In-Migrants`~Infected~and~Symptomatic~frac(dM[it],dt)/M[t] )) +
    scale_x_date(breaks = scales::pretty_breaks(n = 5)) +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 12)) +
    ggsave("Chhattisgarh-mig-c3.pdf", device = "pdf", width = 12, height = 5)





##############################
# check the residuals

# Residuals (It)
residIt <- estList$Prediction[, "TrueIt"] / estList$Prediction[, "It"]
plot(log(residIt[20:160]))

p1 <- forecast::ggtaperedacf(log(residIt[20:160]), type = "correlation")
p1$layers[[2]]$aes_params$fill <- "grey60"
p1$layers[[2]]$aes_params$alpha <- 0.5
p1$layers[[3]]$aes_params$size <- 1
p1$layers[[4]]$data$sig <- !p1$layers[[4]]$data$sig
p1$layers[[4]]$aes_params$size <- 3
p1 <- p1 + ggtitle("") + theme_bw() +
    guides(color = guide_legend(title = "Significant")) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + 
    theme(legend.position = "top")


p2 <- forecast::ggtaperedpacf(log(residIt[20:160]))
p2$layers[[2]]$aes_params$fill <- "grey60"
p2$layers[[2]]$aes_params$alpha <- 0.5
p2$layers[[3]]$aes_params$size <- 1
p2$layers[[4]]$data$sig <- !p2$layers[[4]]$data$sig
p2$layers[[4]]$aes_params$size <- 3
p2 <- p2 + ggtitle("") + theme_bw() +
    guides(color = guide_legend(title = "Significant")) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    theme(legend.position = "top")

ggpubr::ggarrange(p1, p2)
ggsave('Resid_It.pdf', device = "pdf", width = 12, height = 4)

# Residuals (R)
residR <- estList$Prediction[, "TrueR"] / estList$Prediction[, "R"]
plot(log(residR[25:160]))

p1 <- forecast::ggtaperedacf(log(residR[25:160]), type = "correlation")
p1$layers[[2]]$aes_params$fill <- "grey60"
p1$layers[[2]]$aes_params$alpha <- 0.5
p1$layers[[3]]$aes_params$size <- 1
p1$layers[[4]]$data$sig <- !p1$layers[[4]]$data$sig
p1$layers[[4]]$aes_params$size <- 3
p1 <- p1 + ggtitle("") + theme_bw() +
    guides(color = guide_legend(title = "Significant")) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + 
    theme(legend.position = "top")


p2 <- forecast::ggtaperedpacf(log(residR[25:160]))
p2$layers[[2]]$aes_params$fill <- "grey60"
p2$layers[[2]]$aes_params$alpha <- 0.5
p2$layers[[3]]$aes_params$size <- 1
p2$layers[[4]]$data$sig <- !p2$layers[[4]]$data$sig
p2$layers[[4]]$aes_params$size <- 3
p2 <- p2 + ggtitle("") + theme_bw() +
    guides(color = guide_legend(title = "Significant")) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    theme(legend.position = "top")

ggpubr::ggarrange(p1, p2)
ggsave('Resid_R.pdf', device = "pdf", width = 12, height = 4)


###############################
# Residual parameter estimation

# T_t
W1t <- log(estList$Prediction[, "TrueIt"]) - log(estList$Prediction[, "It"])
W2t <- log(estList$Prediction[, "TrueR"]) - log(estList$Prediction[, "R"])
W3t <- log(estList$Prediction[, "TrueE"]) - log(estList$Prediction[, "E"])

fit1 <- arima(W1t[80:160], order = c(3, 1, 0))
fit2 <- arima(W2t[80:160], order = c(3, 1, 0))
fit3 <- arima(W3t[80:160], order = c(3, 1, 0))


# Simulate to find theta_t
get_theta <- function(confirmed) {
    N_urban <- 5937237  # Total urban population of Chattisgarh (http://www.chtenvis.nic.in/pdf/demography.pdf)
    avg_hh_size <- 4.8    # average household size (https://en.wikipedia.org/wiki/Indian_states_ranking_by_household_size)
    low_income_group <- 0.4   # low income percentage (http://documents1.worldbank.org/curated/en/166551468194958356/pdf/105848-BRI-P157572-PUBLIC-Chhattisgarh-Proverty.pdf)
    
    N_adj <- ((N_urban / avg_hh_size) * low_income_group * 2) + ((N_urban / avg_hh_size) * (1 - low_income_group) * 1)
    
    N_total <- 3.23e7
    N_c <- (23007 / 3257)
    
    # Estimate of theta_t
    theta_t <- diff(confirmed) * ( (N_c / N_adj) + (avg_hh_size-2) / (N_total - N_adj))/2
    theta_t <- c(theta_t[1], theta_t)
    return(theta_t)
}


B <- 1000
thetas <- matrix(NA, nrow = B, ncol = 160)
confirmed <- matrix(NA, nrow = B, ncol = 160)

set.seed(12345)
for (b in 1:B) {
    W1t_new <- arima.sim(list("ar" = fit1$coef, "sd" = fit1$sigma2), n = 160)
    W2t_new <- arima.sim(list("ar" = fit2$coef, "sd" = fit2$sigma2), n = 160)
    W3t_new <- arima.sim(list("ar" = fit3$coef, "sd" = fit3$sigma2), n = 160)
    
    It_new <- estList$Prediction[-161, "It"] * exp(W1t_new)
    R_new <- estList$Prediction[-161, "R"] * exp(W2t_new)
    E_new <- estList$Prediction[-161, "E"] * exp(W3t_new)
    confirmed_new <- It_new + R_new + E_new
    
    confirmed[b, ] <- sort(confirmed_new)
    thetas[b, ] <- get_theta(confirmed[b, ])
}

plot(confirmed[1, ])
plot(estList$Prediction[, "It"] + estList$Prediction[, "R"] + estList$Prediction[, "E"])


theta_df <- tibble(
    "Lower" = apply(thetas, 2, FUN = function(x) { quantile(x, probs = 0.025) }),
    "Upper" = apply(thetas, 2, FUN = function(x) { min(2.95e-3, quantile(x, probs = 0.975)) }),
    "Estimate" = estList$theta[-161],
    "Date" = estList$Prediction[-161, "time"]
)

ggplot(theta_df, aes(x = Date)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "blue", alpha = 0.25) +
    geom_line(aes(y = Estimate), color = "blue", size = 1) + 
    ylim(0, 3e-3) + theme_bw() +
    xlab("Date") + ylab(expression(theta[t])) + 
    theme(axis.title = element_text(size = 14)) +
    ggsave("theta_t_prediction.pdf", device = "pdf", width = 6, height = 4)


# CI of delta_t
a <- 0.22315656
N_urban <- 5937237  # Total urban population of Chattisgarh (http://www.chtenvis.nic.in/pdf/demography.pdf)
avg_hh_size <- 4.8    # average household size (https://en.wikipedia.org/wiki/Indian_states_ranking_by_household_size)
low_income_group <- 0.4   # low income percentage (http://documents1.worldbank.org/curated/en/166551468194958356/pdf/105848-BRI-P157572-PUBLIC-Chhattisgarh-Proverty.pdf)

N_adj <- ((N_urban / avg_hh_size) * low_income_group * 2) + ((N_urban / avg_hh_size) * (1 - low_income_group) * 1)
b_t <- 0.1 * (estList$Prediction[, "It"] + estList$Prediction[, "R"] + estList$Prediction[, "E"]) / N_adj

k_hat <- 0.0113
e_hat <- 0.01955
k_var <- 9.8137e-4
e_var <- 3.1419e-3
t <- 2:160
vardelta <- k_var * ((t * b_t[t] / (a * e_hat - b_t[t])) - ((t-1) * b_t[t-1] / (a * e_hat - b_t[t-1]) ) )^2

delta_df <- tibble(
    "Lower" = pmax(estList$delta[t-1] - 1.96 * sqrt(vardelta), 0),
    "Upper" = pmin(estList$delta[t-1] + 1.96 * sqrt(vardelta), 1.5e-2),
    "Estimate" = estList$delta[t-1],
    "Date" = estList$Prediction[1:159, "time"]
)

ggplot(delta_df, aes(x = Date)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "red", alpha = 0.25) +
    geom_line(aes(y = Estimate), color = "red", size = 1) +
    ylim(0, 1.5e-2) + theme_bw() +
    xlab("Date") + ylab(expression(delta[t])) + 
    theme(axis.title = element_text(size = 14)) +
    ggsave("delta_t_prediction.pdf", device = "pdf", width = 6, height = 4)












