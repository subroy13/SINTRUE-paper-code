library(readr)
library(dplyr)
library(tidyr)
library(deSolve)
library(ggplot2)

pred_diag <- function(preds, N, start_date = as.Date('2020-02-20')) {
    preds[, -1] <- round(preds[, -1] * N)
    preds <- as_tibble(preds)
    preds$time <- start_date + preds$time - 1
    
    return(preds)
}

# logit function
logit <- function(x) { log(x/(1 - x)) }

# function for Numerical derivatives
numDeriv <- function(x, order = 1) {
    if (order == 1) {
        diffx <- (lead(x) - lag(x))/2
    } 
    if (order == 2) {
        diffx <- (-lead(x, 2) + 8 * lead(x, 1) - 8 * lag(x, 1) + lag(x, 2)) / 12
    }
    if (order > 2) {
        stop("Not Implemented Yet...")
    }
    
    return(diffx)
}



# Take data and model
dat <- read_csv('./datasets/time_series_covid19_confirmed_global.csv') %>% 
    dplyr::filter(`Country/Region` == "Italy") %>% 
    select(-c(`Province/State`, Lat, Long)) %>%
    tidyr::pivot_longer(-`Country/Region`, names_to = "Date", values_to = "Confirmed") %>%
    mutate(Date = lubridate::mdy(Date)) %>%
    left_join(
        read_csv('./datasets/time_series_covid19_deaths_global.csv') %>% 
        dplyr::filter(`Country/Region` == "Italy") %>% 
        select(-c(`Province/State`, Lat, Long)) %>%
        tidyr::pivot_longer(-`Country/Region`, names_to = "Date", values_to = "Deaths") %>%
        mutate(Date = lubridate::mdy(Date))
    ) %>%
    left_join(
        read_csv('./datasets/time_series_covid19_recovered_global.csv') %>% 
            dplyr::filter(`Country/Region` == "Italy") %>% 
            select(-c(`Province/State`, Lat, Long)) %>%
            tidyr::pivot_longer(-`Country/Region`, names_to = "Date", values_to = "Recovered") %>%
            mutate(Date = lubridate::mdy(Date))
    ) %>%
    mutate(Active = Confirmed - Recovered - Deaths) %>%
    dplyr::filter(Date > as.Date('2020-02-19'))


subdat <- jsonlite::fromJSON('./datasets/italy_covid.json')
subdat$totalCases <- as.numeric(subdat$totalCases)

ggplot(dat, aes(x = Date)) +
    geom_line(aes(y = Active), color = "red") +
    geom_line(aes(y = Recovered), color = "darkgreen")

train_dat <- dat %>% dplyr::filter(Date < as.Date('2020-03-13'))



# SIDARTHE MODEL (For Italy Feb 20 to March 12)
states <- c(S = (60e6 - (200 + 20 + 3))/60e6, I = 200/60e6, D = 20/60e6, 
            A = 1/60e6, R = 2/60e6, T_ = 0, H = 0, E = 0)

parameters <- c(theta = 0.371, zeta = 0.125, eta = 0.125, mu = 0.012, nu = 0.027,
                tau = 0.003, lambda = 0.034, rho = 0.034, kappa = 0.017, xi = 0.017, sigma = 0.017)

alpha_fun <- approxfun(c(rep(0.57, 4), rep(0.422, 18), rep(0.2, 176) ))
gamma_fun <- approxfun(c(rep(0.456, 4), rep(0.285, 18), rep(0.136, 176) ))
betadel_fun <- approxfun(c(rep(0.011, 4), rep(0.0057, 18), rep(0.0057, 176) ))
epsilon_fun <- approxfun(c(rep(0.171, 12), rep(0.143, 10 + 186)))

SIDARTHE <- function(t, states, parameters, alpha_fun, gamma_fun, betadel_fun, epsilon_fun) {
    with(as.list(c(states, parameters)), {
        beta_fun <- betadel_fun
        delta_fun <- betadel_fun
        
        dI <- S * (alpha_fun(t) * I + beta_fun(t) * D + gamma_fun(t) * A + delta_fun(t) * R) - (epsilon_fun(t) + zeta + lambda) * I
        dD <- epsilon_fun(t) * I - (eta + rho) * D
        dA <- zeta * I - (theta + mu + kappa) * A
        dR <- eta * D + theta * A - (nu + xi) * R
        dT <- mu * A + nu * R - (sigma + tau) * T_
        dH <- lambda * I + rho * D + kappa * A + xi * R + sigma * T_
        dE <- tau * T_
        dS <- -(dI + dD + dA + dR + dT + dH + dE)
        
        return(list(c(dS, dI, dD, dA, dR, dT, dH, dE)))
    })
}

preds <- pred_diag(ode(states, times = 1:198, func = SIDARTHE, parms = parameters, 
             alpha_fun = alpha_fun, gamma_fun = gamma_fun, betadel_fun = betadel_fun, 
             epsilon_fun = epsilon_fun), N = 60e6)
#View(preds)

preds$Regime <- c(rep("Observed", 22), rep("Predicted", nrow(preds) - 22 ) )


ggplot(preds, aes(x = time)) + 
    geom_line(aes(y = D + R + T_, color = Regime), size = 1) +
    geom_point(aes(y = dat$Active[1:nrow(preds)], color = Regime)) + 
    scale_x_date(breaks = scales::pretty_breaks(n = 10)) +
    theme_bw() +
    xlab('Date') + ylab('Count of Active Cases')

ggsave("Italy-SIDARTHE.pdf", device = "pdf")

sqrt(mean(((preds$D + preds$R + preds$T_ -  dat$Active)[23:nrow(preds)])^2, na.rm = T))
# 33219.58

################
# Our Model

SINTRUE <- function(t, states, parameters, 
                  beta1_fun, theta_fun, delta_fun) {
    with(as.list(c(states, parameters)), {
        
        beta2_fun <- function(x) { beta1_fun(x) * 0.617 }
        
        dIp <- beta1_fun(t) * S * Ip + beta2_fun(t) * S * Isn - (alpha + lambda + theta_fun(t)) * Ip
        dIsn <- alpha * Ip - (kappa + delta_fun(t)) * Isn
        dU <- kappa * Isn + lambda * Ip    # U = unrecognized people who becomes good from symptomatic
        dIt <- theta_fun(t) * Ip + delta_fun(t) * Isn - (tau + gamma) * It
        dR <- gamma * It
        dE <- tau * It
        
        dS <-  - (dIp + dIsn + dIt + dR + dE + dU)
        
        return(list(c(dS, dIp, dIsn, dU, dIt, dE, dR)))
        
    })
}

# function to estimate Tau
est.TauGamma <- function(It, R, E) {
    dE <- numDeriv(E)
    fit <- MASS::rlm(dE ~ It - 1, maxit = 100)
    
    plot(It, dE, pch = 19)
    abline(a = 0, b = fit$coefficients, col = "blue")
    
    tau <- fit$coefficients
    
    dR <- numDeriv(R)
    fit <- MASS::rlm(dR ~ It - 1, maxit = 100)
    
    plot(It, dR, pch = 19)
    abline(a = 0, b = fit$coefficients, col = "blue")
    
    gam <- fit$coefficients
    
    return(c("Tau" = tau, "Gamma" = gam))
}

est.TauGamma(train_dat$Active, train_dat$Recovered, train_dat$Deaths)
#est.TauGamma(dat$Active, dat$Recovered, dat$Deaths)

N_c <- 14
covid_prevalence <- 0.02

theta_t <- c(subdat$totalCases[1], diff(subdat$totalCases)) * (N_c / 60e6)
theta_t <- pmax(theta_t, 0)
plot(theta_t)

p_testSymp <- (subdat$hospitalizedWithSymptoms / subdat$totalHospitalized) * 
    (as.numeric(subdat$totalCases)/60e6) / (0.822 * 0.617 * 0.44 * covid_prevalence)
plot(p_testSymp)

plot_our_model <- function(parameters, beta1_fun) {
    states <- c(S = (60e6 - (200 + 20 + 3 + 50))/60e6, 
                Ip = 200/60e6, Isn = 1/60e6,  U = 50/60e6, It = 20/60e6, R = 2/60e6, E = 0)
    
    
    #delta_t_int <- parameters["kappa"] * c(1:length(p_testSymp)) * p_testSymp / (1 - p_testSymp)
    #delta_t <- pmax(diff(delta_t_int), 0)
    #delta_t <- c(delta_t[1:2], delta_t)
    delta_t <- parameters["kappa"] * p_testSymp / (1 - p_testSymp)
    plot(delta_t)
    
    theta_fun = approxfun(c(theta_t[1], theta_t)) 
    delta_fun = approxfun(delta_t)
    preds <- ode(states, times = 1:198, func = SINTRUE, parms = parameters, 
        beta1_fun = beta1_fun, theta_fun = theta_fun, delta_fun = delta_fun)
    
    #View(preds)
    preds <- pred_diag(preds, N = 60e6)
    preds$Regime <- c(rep("Observed", 22), rep("Predicted", nrow(preds) - 22 ) )
    
    p <- ggplot(preds, aes(x = time)) + 
        geom_line(aes(y = It, color = Regime), size = 1) +
        geom_point(aes(y = dat$Active[1:nrow(preds)], color = Regime)) + 
        scale_x_date(breaks = scales::pretty_breaks(n = 10)) +
        theme_bw() +
        xlab('Date') + ylab('Count of Active Cases')
    
    ggsave("Italy-SINTRUE.pdf", device = "pdf")
    
    cat('RMSE: ', sqrt(mean(((preds$It -  dat$Active[1:198])[23:nrow(preds)])^2, na.rm = T)),'\n')
    
    p
}

parameters <- c(tau = 0.01719181, gamma = 0.01567709,
                kappa = 0.15, alpha = 0.012, lambda = 0.081)
beta1_fun <- function(x) { return(exp(-0.0715 * x)) }

plot_our_model(parameters, beta1_fun)




















