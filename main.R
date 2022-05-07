library(dplyr)
library(readr)
library(tidyr)
library(EpiEstim)
library(deSolve)
library(ggplot2)
library(forecast)

setwd("D:/Academic/Projects/2020_Covid 19 Prediction Model_Anirban Ghatak_IIMV/")


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


# Extract State level data
# updated Date
upto_date <- as.Date('2020-08-22')

statedf <- read_csv('./datasets/india-state-data.csv')
statedf <- statedf %>% dplyr::filter(Date < upto_date)   # today's data is NA

CGdf <- statedf %>% dplyr::filter(State == "Chhattisgarh")   # consider only Chattishgarh's data

# Extract District level data to subset migrant data
distdf <- read_csv('./datasets/india-district-level.csv')

migrant_terms<- c( "Others","Other", "Other Region", "Other State", "Other States",  "Unknown","Unassigned", "BSF Camp", "Airport Quarantine", "Railway Quarantine")

migdf <- distdf %>% 
    dplyr::filter(District %in% migrant_terms) %>% 
    group_by(Date, State) %>%
    summarise(Confirmed = sum(Confirmed), Active = sum(Active), Recovered = sum(Recovered), Deaths = sum(Deaths)) %>%
    dplyr::filter(Active >= 0)  # remove all rows where active becomes < 0 (because of late publication of state bulletin)

times <- seq(statedf$Date[1], max(statedf$Date), by = "day")

###############################
# ESTIMATION of R0_t for a specified state

estimate.R0 <- function(state_name) {
    
    tmpdf <- statedf %>% dplyr::filter(State == state_name) %>% arrange(Date)
    tmpdf$New_Confirmed <- c(0, diff(tmpdf$Confirmed))
    
    epicurve <- tmpdf$New_Confirmed
    names(epicurve) <- tmpdf$Date
    
    w <- 14   # window length
    t_start <- seq(2, nrow(tmpdf)-w)
    t_end <- t_start + w
    
    a <- estimate_R(tmpdf$New_Confirmed, method = "parametric_si", 
                    config = make_config(list(mean_si = 4.7, std_si = 2.9, t_start = t_start, t_end = t_end)))
    
    R0 <- a$R[, "Mean(R)"]
    names(R0) <- tmpdf$Date[a$R[, "t_end"]]
    
    R0_all <- rep(NA, nrow(tmpdf))
    names(R0_all) <- as.character(tmpdf$Date)
    
    R0_all[names(R0)] <- R0
    R0_all[1:(w+1)] <- R0_all[w+2]
    return(R0_all[as.character(tmpdf$Date)])   # return the required portion
}


#####################################
# ESTIMATION of c1, c2, c3

estimate.migrant <- function() {
    # percentage of migrants from each state in Chattisgarh is taken same as
    # percentage of migrants from each state into Chattisgarh who have been infected.
    # If we get data about the former we can replace it.
    
    # origin states in order:
    
    # MH, UP, Delhi, TG, Guj, MP, TN, Haryana, Odisha, AP, JK, Ktka, Bihar, J&K, RJ,UK
    # Punjab, Assam, Goa, HP, Kerala
    
    migrant_incoming_rate <- c(466,127,58,63,43,34,28,22,23,16,18,10,7,8,6,3,3,2,1,1,1,2)/944
    migrant_incoming_rate<- migrant_incoming_rate[1:12] / sum(migrant_incoming_rate[1:12])
    
    # we consider first 12 origin states
    # for each significant state take its capital city's population+
    # some other important city (if applicable) 
    # as most migrants would work in those cities
    
    origin_states<- c("Maharashtra","Uttar Pradesh","Delhi","Telangana","Gujarat","Madhya Pradesh",
                      "Tamil Nadu","Haryana","Odisha","Andhra Pradesh","Jharkhand","Karnataka")
    names(migrant_incoming_rate) <- origin_states
    
    # mumbai+pune, kanpur+lakhnow, delhi, hyderabad, ahmedabad+vadodara,bhopal,
    # chennai, chandigarh, bhubaneswar, vizag+other cities, ranchi, bangalore
    total_popn <- c(1.84e7+31.2e5,31.24e5+28.2e5, 1.9e7, 1e7, 55.7e5+21.8e5, 18e5 ,70.9e5, 10.6e5, 8.38e5,30e5,10e5, 1.2e7)
    names(total_popn) <- origin_states
    
    # assuming migrants take 3 days to come on average
    migrant_popn <- c(7018877, 18000, 6352546, 101146, 4439145, 2000, 2.9e5, 3507625, 26621, 27548, 1000, 5439221)
    names(migrant_popn) <- origin_states
    
    
    # This function calculates infected proportion among migrant in an origin state
    infected.origin <- function(origin_state) {
        # extract number of active time series
        tmp_state <- statedf %>% dplyr::filter(State == origin_state)
        
        mig_count <- rep(NA, length(times))     # this holds the time series of migration counts
        names(mig_count) <- times
        
        tmp_mig <- migdf %>% dplyr::filter(State == origin_state)
        
        # if there is enough data on migrant, use migrant based formula, otherwise overall formula
        if (nrow(tmp_mig) > 30) {
            migrant.inf.rate <- tmp_mig$Active/migrant_popn[origin_state]
            mig_count[as.character(tmp_mig$Date)] <- migrant.inf.rate
        } else {
            popn.inf.rate <- tmp_state$Active/total_popn[origin_state]
            mig_count[as.character(tmp_state$Date)] <- popn.inf.rate
        }
        
        # approximate this
        n_time <- as.numeric(times[length(times)] - times[1]) + 1
        mig_count <- sapply(1:n_time, approxfun(mig_count))
        
        mig_count[is.na(mig_count)] <- 0  # the first part where approximation is not possible is set to 0
        names(mig_count) <- times
        
        return(mig_count)
    }
    
    # c3 = sum(infected proportion among migrants in origin state * proportion of people coming from origin state)
    infected.origin.all_states <- sapply(origin_states, FUN = infected.origin)  # each column is statewise series
    c3 <- rowSums(infected.origin.all_states * matrix(migrant_incoming_rate[origin_states], nrow = nrow(infected.origin.all_states), ncol = ncol(infected.origin.all_states), byrow = T ))   # EQN OF c3
    var.c3 <- rowSums(infected.origin.all_states^2 * matrix(migrant_incoming_rate[origin_states] , nrow = nrow(infected.origin.all_states), ncol = ncol(infected.origin.all_states), byrow = T ))
    for (t in 1:length(var.c3)) {
        var.c3[t] <- var.c3[t] - sum(outer(migrant_incoming_rate[origin_states], migrant_incoming_rate[origin_states]) *
                               outer(infected.origin.all_states[t, origin_states], infected.origin.all_states[t, origin_states]) )
    }
    
    
    # Estimation of c2 
    R0.all_states <- sapply(origin_states, FUN = estimate.R0)
    
    # need to add a lag of 3 days
    R0.all_states <- rbind(R0.all_states[1:3, ], R0.all_states[-( (nrow(R0.all_states)-2):nrow(R0.all_states) ), ])
    
    c2 <- rowSums(infected.origin.all_states * 
                      matrix(migrant_incoming_rate[origin_states], nrow = nrow(infected.origin.all_states), ncol = ncol(infected.origin.all_states), byrow = T ) *
                      R0.all_states )   # EQN OF c2
    var.c2 <- rowSums(infected.origin.all_states^2 * 
                          matrix(migrant_incoming_rate[origin_states], nrow = nrow(infected.origin.all_states), ncol = ncol(infected.origin.all_states), byrow = T ) *
                          R0.all_states^2)
    for (t in 1:length(var.c2)) {
        var.c2[t] <- var.c2[t] - sum(outer(migrant_incoming_rate[origin_states], migrant_incoming_rate[origin_states]) *
                                         outer(infected.origin.all_states[t, origin_states], infected.origin.all_states[t, origin_states]) *
                                         outer(R0.all_states[t, origin_states], R0.all_states[t, origin_states] ))
    }
    
    u_mat <- 1 - infected.origin.all_states * (1 + R0.all_states)
    c1 <- (1 - c2 - c3)  # EQN of c1
    var.c1 <- rowSums(u_mat^2 * matrix(migrant_incoming_rate[origin_states], nrow = nrow(infected.origin.all_states), ncol = ncol(infected.origin.all_states), byrow = T ))
    for (t in 1:length(var.c1)) {
        var.c1[t] <- var.c1[t] - sum(outer(migrant_incoming_rate[origin_states], migrant_incoming_rate[origin_states]) *
                                         outer(u_mat[t, origin_states], u_mat[t, origin_states]) )
    }
    
    
    # finally output everything
    return(list(c1 = c1, var.c1 = var.c1, c2 = c2, var.c2 = var.c2, c3 = c3, var.c3 = var.c3))
}

mig <- estimate.migrant()

# Plot it
par(mfrow = c(1, 3))
plot(as.Date(names(mig$c1)), mig$c1, type = "l", col = "green")
plot(as.Date(names(mig$c2)), mig$c2, type = "l", col = "brown")
plot(as.Date(names(mig$c3)), mig$c3, type = "l", col = "red")
par(mfrow = c(1, 1))


###############################################
# ESTIMATION OF theta_t

# before 20th March, no asymptotic is tested (by ICMR report on testing strategies)
# we assume a logistic like increasing curve and assume it is capped between 0.01 and 0.9
asym.incoming <- function(d) {
    x <- as.numeric(d - as.Date('2020-03-20')) / as.numeric(as.Date('2020-06-17') - as.Date('2020-03-20'))
    x <- ifelse(x < 0.01, 0.01, x)
    x <- logit(0.01) + x * ( logit(0.89) - logit(0.01) )
    out <- exp(x)/(1 + exp(x))
    out <- ifelse(out > 0.9, 0.9, out)
    return(out)
}

N_urban <- 5937237  # Total urban population of Chattisgarh (http://www.chtenvis.nic.in/pdf/demography.pdf)
avg_hh_size <- 4.8    # average household size (https://en.wikipedia.org/wiki/Indian_states_ranking_by_household_size)
low_income_group <- 0.4   # low income percentage (http://documents1.worldbank.org/curated/en/166551468194958356/pdf/105848-BRI-P157572-PUBLIC-Chhattisgarh-Proverty.pdf)

N_adj <- ((N_urban / avg_hh_size) * low_income_group * 2) + ((N_urban / avg_hh_size) * (1 - low_income_group) * 1)

N_total <- 3.23e7
N_c <- (23007 / 3257)

# Estimate of theta_t
theta_t <- diff(CGdf$Confirmed) * ( (N_c / N_adj) + (avg_hh_size-2) / (N_total - N_adj))/2
plot(CGdf$Date[-1], theta_t, type = "l")

####################
# Estimate of Tau

# function to estimate Tau
est.Tau <- function(It, E) {
    dE <- numDeriv(E)
    fit <- MASS::rlm(dE ~ It - 1, maxit = 100)
    
    #plot(It, dE, pch = 19)
    #abline(a = 0, b = fit$coefficients, col = "blue")
    tab <- summary(fit)
    print(tab$coefficients)
    
    tau <- fit$coefficients
    names(tau) <- "Tau"
    return(tau)
}

tau_est <- est.Tau(It = CGdf$Active, E = CGdf$Deaths)
tau_est

# function to estimate Gamma
est.Gamma <- function(It, R, f = 0) {
    dR <- numDeriv(R)
    fit <- MASS::rlm(dR ~ It - 1, maxit = 100)
    
    #plot(It, dR, pch = 19)
    #abline(a = 0, b = fit$coefficients, col = "blue")
    tab <- summary(fit)
    print(tab$coefficients)
    
    gam <- fit$coefficients * (1 + f)
    names(gam) <- "Gamma"
    return(gam)
}

gamma_est <- est.Gamma(It = CGdf$Active, R = CGdf$Recovered)
gamma_est

# tau and gamma for different states
subdf <- statedf %>% dplyr::filter(State == "Chhattisgarh")
est.Tau(subdf$Active, subdf$Deaths)
est.Gamma(subdf$Active, subdf$Recovered)


######################################
# Now consider estimate of P(tested | asymp & covid)
#p_symp <- 1 - (1 - 0.822) * (1 - 0.617) * (1 - 0.44)
p_symp <- 0.822 * 0.617 * 0.44
covid_prevalence <- (1007/23007)
#covid_prevalence <- 0.1955

par(mfrow = c(1, 2))
p_testAsymp <- (cumsum(asym.incoming(CGdf$Date[-1]) * diff(CGdf$Confirmed)) /N_adj) / ( (1 - p_symp) * covid_prevalence)
plot(CGdf$Date[-1], p_testAsymp, type = "l")

# consider, P(tested | symp & covid)
p_testsymp <- ( cumsum((1 - asym.incoming(CGdf$Date[-1])) * diff(CGdf$Confirmed)) /N_adj) /( p_symp * covid_prevalence)
plot(CGdf$Date[-1], p_testsymp, type = "l")

par(mfrow = c(1, 1))


#####################
# Try Simulation

model <- function(t, states, parameters, 
                  c1_fun, c2_fun, c3_fun, 
                  beta1_fun, beta2_fun,
                  theta_fun, delta_fun) {
    with(as.list(c(states, parameters)), {
        
        dN <- In_rate - Out_rate
        dM_n <- c1_fun(t) * In_rate
        dM_p <- c2_fun(t) * In_rate
        dM_i <- c3_fun(t) * In_rate
        
        dIp <- beta1_fun(t) * S * Ip / N + beta2_fun(t) * S * Isn / N + dM_p - (alpha + lambda + theta_fun(t)) * Ip
        dIsn <- alpha * Ip - (kappa + delta_fun(t)) * Isn
        dU <- kappa * Isn + lambda * Ip    # U = unrecognized people who becomes good from symptomatic
        dIt <- theta_fun(t) * Ip + delta_fun(t) * Isn + dM_i - (tau + gamma) * It
        dR <- (gamma / (1 + f)) * It
        dE <- tau * It
        
        dS <- (dN - dIp - dIsn - dIt - dR - dE - dU)
        
        return(list(c(dN, dS, dIp, dIsn, dU, dIt, dE, dR)))
        
    })
}

#############
# Perform some kind of LS / Grid Search

LS <- function(params) {
    kappa <- params["kappa"]
    alpha <- params["alpha"]
    lambda <- params["lambda"]
    beta_1.pre <- params["beta_1.pre"]   # PRE -migration initiation
    beta_1.post <- params["beta_1.post"]   # post migration initiation
    #beta_2.pre <- 0.617 * params["beta_1.pre"]  # pre migration
    #beta_2.post <- 0.617 * params["beta_1.post"]  
    
    delta_t_int <- kappa * c(1:length(p_testsymp)) * p_testsymp / (1 - p_testsymp)
    delta_t <- pmax(diff(delta_t_int), 0)
    delta_t <- c(delta_t[1:2], delta_t)
    
    parameters <- c(In_rate = 15000, Out_rate = 5000, 
                    f = 0, kappa = as.numeric(kappa), alpha = as.numeric(alpha), lambda = as.numeric(lambda), 
                    tau = as.numeric(tau_est), gamma = as.numeric(gamma_est))
    
    # pass in the initial states as required
    init.states <- c(N = N_total,
                     S = N_total - 11, 
                     Ip = 10, Isn = 1, U = 0, It = 0,
                     E = 0, R = 0)
    
    # Beta functions
    pre.times <- 1:(as.numeric(as.Date('2020-05-01') - min(times)) - 28)
    post.times <- (as.numeric(as.Date('2020-05-01') - min(times)) + 28):length(times)
    
    beta1_fun <- splinefun(x = c(pre.times, post.times), 
                           y = c(exp(-beta_1.pre * pre.times) + beta_1.post, rep(beta_1.post, length(post.times)) ) )
    beta2_fun <- function(x) { 0.617 * beta1_fun(x) }
                              
    #beta1_fun <- splinefun(x = c(pre.times, post.times), y = c(rep(beta_1.pre, length(pre.times)), rep(beta_1.post, length(post.times)) ) )
    #beta2_fun <- splinefun(x = c(pre.times, post.times), y = c(rep(beta_2.pre, length(pre.times)), rep(beta_2.post, length(post.times)) ) )
    #beta1_fun <- stepfun(x = as.numeric(as.Date('2020-05-01') - min(times)), y = c(beta_1.pre, beta_1.post))
    #beta2_fun <- stepfun(x = as.numeric(as.Date('2020-05-01') - min(times)), y = c(beta_2.pre, beta_2.post))
    
    out <- ode(init.states, 1:length(times), model, parameters, 
               c1_fun = approxfun(mig$c1), c2_fun = approxfun(mig$c2), c3_fun = approxfun(mig$c3), 
               beta1_fun = beta1_fun, 
               beta2_fun = beta2_fun,
               theta_fun = approxfun(c(theta_t[1], theta_t)), 
               delta_fun = approxfun(delta_t))
    
    # check if theta_t * I_p actually matches
    error.1 <- sum( (log1p(out[-1, "Ip"]) - log1p(  asym.incoming(CGdf$Date[-1]) * diff(CGdf$Confirmed)/theta_t) )^2, na.rm = T)
    
    # check if delta_t * Isn actually matches
    error.2 <- sum( (log1p(out[-1, "Isn"]) - log1p( (1 - asym.incoming(CGdf$Date[-1])) * diff(CGdf$Confirmed)/delta_t[-1] ) )^2, na.rm = T)
    
    # check if I_t, R, E actually matches
    error.3 <- sum( (log1p(out[, "It"]) - log1p(CGdf$Active))^2, na.rm = T )
    error.4 <- sum( (log1p(out[, "R"]) - log1p(CGdf$Recovered) )^2, na.rm = T )
    error.5 <- sum( (log1p(out[, "E"]) - log1p(CGdf$Deaths) )^2, na.rm = T )
    error_vect <- c(error.1, error.2, error.3, error.4, error.5)
    
    error <- sum(error_vect * c(1, 1, 10, 10, 10))
    
    # Prediction
    preds <- as.data.frame(round(out))
    preds$time <- times
    preds$TrueIt <- CGdf$Active
    preds$TrueR <- CGdf$Recovered
    preds$TrueE <- CGdf$Deaths
    
    # Compute R0
    beta1_t <- sapply(1:length(times), beta1_fun)
    beta2_t <- sapply(1:length(times), beta2_fun)
    
    R0 <- (beta1_t / (lambda + alpha + c(theta_t[1], theta_t) )) * preds$Ip/(preds$Ip + preds$Isn + preds$It) + 
        (beta2_t / (delta_t + kappa)) * preds$Isn / (preds$Ip + preds$Isn + preds$It)
    preds$R0 <- R0
    
    return(list("Error" = mean(error_vect) + sd(error_vect), 
                "Pred" =  preds, "Error_Vect" = error_vect, "delta" = delta_t))
}


# Run this for different parameter combination,
# DO A GRID SEARCH
library(pbapply)
grid_list = as.matrix(expand.grid(kappa = seq(0.01, 0.51, by = 0.05), 
                        alpha = 1/5 * 0.06, 
                        lambda = seq(0.1, 0.5, by = 0.1), 
                        beta_1.pre = seq(0.01, 0.21, by = 0.025),
                        beta_1.post = seq(0.01, 0.031, by = 0.025) ))
nrow(grid_list)

value_vect <- pbsapply(1:nrow(grid_list), FUN = function(i){
    LS(grid_list[i, ])$Error_Vect
})

value_sd <- matrixStats::colSds(value_vect)
value_means <- matrixStats::colMeans2(value_vect)

value_maxs <- matrixStats::colMaxs(value_vect)
which.min(value_maxs)

min(as.numeric(value_vect))
which.min(as.numeric(value_vect))
grid_list[which.min(as.numeric(value_vect)), ]


# Did a grid search, this seems really optimal parameter
#param_vect <- c(kappa = 0.0043, alpha =  1/5 * 0.06, lambda = 0.083, 
#                beta_1.pre = 0.102, beta_1.post = 0.083)

param_vect <- c(kappa = 0.0113, alpha =  0.012, lambda = 0.079, 
                beta_1.pre = 0.102, beta_1.post = 0.083)



a <- LS(param_vect)

tail(a$Pred)

a$Error_Vect
a$Error

ggplot(a$Pred) +
    geom_line(aes(x = time, y = It), color = "red") +
    geom_point(aes(x = time, y = TrueIt), color = "red") +
    geom_line(aes(x = time, y = R), color = "blue") +
    geom_point(aes(x = time, y = TrueR), color = "blue") +
    theme_bw()

tail(a$Pred$Ip, 10)

plot(a$Pred$time, a$Pred$Ip, type = "l")

View(a$Pred)
plot(a$Pred$time, (a$Pred$Ip + a$Pred$Isn + a$Pred$It + a$Pred$R + a$Pred$U + a$Pred$E) / a$Pred$N, type = "l")
plot(a$Pred$time, a$Pred$Isn, type = "l")
plot(a$Pred$time, a$Pred$R0, type = "l")

# Store everything in an object
ob <- list("kappa" = param_vect["kappa"],
           "alpha" = param_vect["alpha"], 
           "lambda" = param_vect['lambda'],
           'a' = param_vect['beta_1.pre'],
           'b' = param_vect['beta_1.post'], 
           'tau' = tau_est, 'gamma' = gamma_est,
           'theta' = c(theta_t[1], theta_t), 'delta' = a$delta,
           'c1' = mig$c1, 'c2' = mig$c2, 'c3' = mig$c3, 
           'varc1' = mig$var.c1, 'varc2' = mig$var.c2, 'varc3' = mig$var.c3,
           "Prediction" = a$Pred)
saveRDS(obs, './Datasets/estimate.Rds')


################
# Make predictions

predict_model <- function(params, n_new = 50, method = c("normal", "normal", "normal") ) {
    
    kappa <- params["kappa"]
    alpha <- params["alpha"]
    lambda <- params["lambda"]
    beta_1.pre <- params["beta_1.pre"]   # PRE -migration initiation
    beta_1.post <- params["beta_1.post"]   # post migration initiation
    
    new_times <- times[1] + c(1:(length(times) + n_new)) - 1
    
    # extend time varying quantities
    theta_mod <- predict(auto.arima(theta_t, allowdrift = FALSE, allowmean = FALSE), n.ahead = n_new)
    
    if (method[1] == "low") {
        theta_fun <- approxfun(c(theta_t[1], theta_t, as.numeric(theta_mod$pred) - 1.64 * as.numeric(theta_mod$se) ) )
    } else if (method[1] == "high") {
        theta_fun <- approxfun(c(theta_t[1], theta_t, as.numeric(theta_mod$pred) + 1.64 * as.numeric(theta_mod$se) ) )
    } else {
        theta_fun <- approxfun(c(theta_t[1], theta_t, as.numeric(theta_mod$pred) ))
    }
    
    c2_mod <- predict(auto.arima(mig$c2, allowdrift = FALSE, allowmean = FALSE), n.ahead = n_new)
    c3_mod <- predict(auto.arima(mig$c3, allowdrift = FALSE, allowmean = FALSE), n.ahead = n_new)
    
    if (method[2] == "low") {
        c2_fun <- approxfun(c(mig$c2, as.numeric(c2_mod$pred) - 1.64 * as.numeric(c2_mod$se) ))
    } else if (method[2] == "high") {
        c2_fun <- approxfun(c(mig$c2, as.numeric(c2_mod$pred) + 1.64 * as.numeric(c2_mod$se) ))
    } else {
        c2_fun <- approxfun(c(mig$c2, as.numeric(c2_mod$pred) ))
    }
    
    if (method[2] == "low") {
        c3_fun <- approxfun(c(mig$c3, as.numeric(c3_mod$pred) - 1.64 * as.numeric(c3_mod$se) ))
    } else if (method[2] == "high") {
        c3_fun <- approxfun(c(mig$c3, as.numeric(c3_mod$pred) + 1.64 * as.numeric(c3_mod$se) ))
    } else {
        c3_fun <- approxfun(c(mig$c3, as.numeric(c3_mod$pred) ))
    }
    
    c1_fun <- function(x){ 1 - c2_fun(x) - c3_fun(x) }
    
    testsymp_mod <- predict(auto.arima(p_testsymp, allowdrift = FALSE, allowmean = FALSE), n.ahead = n_new)    
    
    if (method[3] == "low") {
        p_testsymp_new <- as.numeric(testsymp_mod$pred) - 1.64 * as.numeric(testsymp_mod$se)
    } else if (method[3] == "high") {
        p_testsymp_new <- as.numeric(testsymp_mod$pred) + 1.64 * as.numeric(testsymp_mod$se)
    } else {
        p_testsymp_new <- as.numeric(testsymp_mod$pred)
    }
    
    delta_t_int <- kappa * c(1:(length(p_testsymp) + n_new) ) * c(p_testsymp, p_testsymp_new) / (1 - c(p_testsymp, p_testsymp_new))
    delta_t <- pmax(diff(delta_t_int), 0)
    delta_t <- c(delta_t[1:2], delta_t)
    
    parameters <- c(In_rate = 15000, Out_rate = 5000, 
                    f = 0, kappa = as.numeric(kappa), alpha = as.numeric(alpha), lambda = as.numeric(lambda), 
                    tau = as.numeric(tau_est), gamma = as.numeric(gamma_est))
    
    # pass in the initial states as required
    init.states <- c(N = N_total,
                     S = N_total - 11, 
                     Ip = 10, Isn = 1, U = 0, It = 0,
                     E = 0, R = 0)
    
    # Beta functions
    pre.times <- 1:(as.numeric(as.Date('2020-05-01') - min(new_times)) - 28)
    post.times <- (as.numeric(as.Date('2020-05-01') - min(new_times)) + 28):length(new_times)
    
    beta1_fun <- splinefun(x = c(pre.times, post.times), 
                           y = c(exp(-beta_1.pre * pre.times) + beta_1.post, rep(beta_1.post, length(post.times)) ) )
    beta2_fun <- function(x) { 0.617 * beta1_fun(x) }
    
    out <- ode(init.states, 1:length(new_times), model, parameters, 
               c1_fun = c1_fun, c2_fun = c2_fun, c3_fun = c3_fun, 
               beta1_fun = beta1_fun, 
               beta2_fun = beta2_fun,
               theta_fun = theta_fun, 
               delta_fun = approxfun(delta_t))
    
    # Prediction
    preds <- as.data.frame(round(out))
    preds$time <- new_times
    preds$TrueIt <- c(CGdf$Active, rep(NA, n_new))
    preds$TrueR <- c(CGdf$Recovered, rep(NA, n_new))
    preds$TrueE <- c(CGdf$Deaths, rep(NA, n_new))
    
    # Compute R0
    beta1_t <- sapply(1:length(new_times), beta1_fun)
    beta2_t <- sapply(1:length(new_times), beta2_fun)
    
    
    R0 <- (beta1_t / (lambda + alpha + c(theta_t[1], theta_t, as.numeric(theta_mod$pred)) )) * preds$Ip/(preds$Ip + preds$Isn + preds$It) + 
        (beta2_t / (delta_t + kappa)) * preds$Isn / (preds$Ip + preds$Isn + preds$It)
    preds$R0 <- R0
    
    
    return(preds)
}

# Did a grid search, this seems really optimal parameter
preds <- predict_model(param_vect, n_new = 200)

#preds <- predict_model(c(kappa = 0.026, alpha =  1/5 * 0.04, lambda = 0.103, 
#                         beta_1.pre = 0.096, beta_1.post = 0.118), n_new = 180)

#preds <- predict_model(c(kappa = 0.196, alpha =  0.06 * 1/5, lambda = 0.145, 
#                         beta_1.pre = 0.089, beta_1.post = 0.157), n_new = 200)

ggplot(preds, aes(x = time)) +
    geom_line(aes(y = It), color = "red") +
    geom_point(aes(y = TrueIt), color = "red") +
    geom_line(aes(y = R), color = "blue") +
    geom_point(aes(y = TrueR), color = "blue") +
    geom_line(aes(y = Ip), color = "green") +
    geom_line(aes(y = Isn), color = "magenta")+
    theme_bw()

##################################################
# PLOTMAKERS

library(gridExtra)

tmpdf <- statedf %>% filter(State == "Chhattisgarh")
est.Tau(tmpdf$Active, E = tmpdf$Deaths)
est.Gamma(tmpdf$Active, tmpdf$Recovered)

p1 <- ggplot(tmpdf, aes(x = Active)) + 
    geom_point(aes(y = numDeriv(Deaths)), color = "red") + 
    stat_smooth(aes(y = numDeriv(Deaths)), color = "red", method = MASS::rlm) + 
    ggtitle("Number of Active cases vs \nRate of change of Deaths") +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5))

p2 <- ggplot(tmpdf, aes(x = Active)) + 
    geom_point(aes(y = numDeriv(Recovered)), color = "darkgreen") + 
    stat_smooth(aes(y = numDeriv(Recovered)), color = "darkgreen", method = MASS::rlm) +
    ggtitle("Number of Active cases vs \nRate of change of Recoveries") +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5))

grid.arrange(p1, p2, ncol = 2)
ggsave(file="Maharashtra-tau-gamma.pdf", arrangeGrob(p1, p2, ncol = 2), 
       width = 12, height = 4, device = "pdf", 
       path = "D:/Dataprojects/Covid-19-Anirban Ghatak/Code/figures") 

# Migrant
tibble(Date = times, c1 = mig$c1, c2 = mig$c2, c3 = mig$c3) %>%
    ggplot(aes(x = Date)) +
    geom_line(aes(y = c3)) +
    ylab(expression(c[3])) +
    ggtitle(expression(Time~Varying~Proportion~of~New~`In-Migrants`~Infected~and~Symptomatic~frac(dM[i],dt)/frac(d(In),dt) )) +
    scale_x_date(breaks = scales::pretty_breaks(n = 5)) +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5)) +
    ggsave(file = "Chhattisgarh-mig-c3.pdf", width = 12, height = 4, device = "pdf", 
           path = "D:/Dataprojects/Covid-19-Anirban Ghatak/Code/figures")



# Theta_t & Delta_t
kappa <- 0.0113
delta_t_int <- kappa * c(1:length(p_testsymp)) * p_testsymp / (1 - p_testsymp)
delta_t <- pmax(diff(delta_t_int), 0)
delta_t <- c(delta_t[1:2], delta_t)

p1 <- tibble(Date = times, theta_t = c(theta_t[1], theta_t), delta_t = delta_t) %>%
    ggplot(aes(x = Date)) +
    geom_line(aes(y = theta_t), color = "blue", size = 1) +
    ylab(expression(theta[t])) +
    ggtitle("Time Varying Testing Rate of\nPre-symptomatic patients") +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5))

p2 <- tibble(Date = times, theta_t = c(theta_t[1], theta_t), delta_t = delta_t) %>%
    ggplot(aes(x = Date)) +
    geom_line(aes(y = delta_t), color = "red", size = 1) +
    ylab(expression(delta[t])) +
    ggtitle("Time Varying Testing Rate of\nSymptomatic patients") +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5))

gridExtra::grid.arrange(p1, p2, ncol = 2)


ggsave(file = "Chhattisgarh-theta-delta.pdf", gridExtra::arrangeGrob(p1, p2, ncol = 2), 
       width = 12, height = 4, device = "pdf", 
       path = "D:/Dataprojects/Covid-19-Anirban Ghatak/Code/figures")

# Fitting
a <- LS(c(kappa = 0.0113, alpha =  1/5 * 0.06, lambda = 0.083, 
          beta_1.pre = 0.102, beta_1.post = 0.083))

pred_date <- as.Date('2020-08-01')
ggplot(a$Pred, aes(x = time)) +
    geom_line(aes(y = It, color = "Hospitalized"), data = a$Pred %>% dplyr::filter(time < pred_date), size = 1) + 
    geom_line(aes(y = R, color = "Recovered"), data = a$Pred %>% dplyr::filter(time < pred_date), size = 1) + 
    geom_line(aes(y = It, color = "Hospitalized", linetype = "Hospitalized"), data = a$Pred %>% dplyr::filter(time >= pred_date), size = 1) + 
    geom_line(aes(y = R, color = "Recovered", linetype = "Recovered"), data = a$Pred %>% dplyr::filter(time >= pred_date), size = 1) + 
    geom_point(aes(y = TrueIt, fill = "Hospitalized"), color = "red") + 
    geom_point(aes(y = TrueR, fill = "Recovered"), color = "darkgreen") + 
    scale_color_manual(name = "Fitted Counts", values = c("Hospitalized"  = "red", "Recovered" = "darkgreen"), guide = "legend") +
    scale_fill_manual(name = "True Counts", values = c(1, 1), guide=guide_legend(override.aes = list(colour=c("darkgreen", "red"))) ) +
    scale_linetype_manual(name = "Predicted Counts", values = c("dotted", "dotted"), guide=guide_legend(override.aes = list(colour=c("darkgreen", "red"), linetype = "dotted")) ) +
    theme_bw() + xlab("Date") + ylab("Count") #+
    ggsave(file = "Chhattisgarh-fitted.pdf", width = 12, height = 4, device = "pdf", 
           path = "D:/Dataprojects/Covid-19-Anirban Ghatak/Code/figures")

R0_vec <- estimate.R0("Chhattisgarh")
R0_vec[1:12] <- NA

ggplot(a$Pred) +
    geom_line(aes(x = time, y = R0, color = "Our model"), size = 1) +
    geom_line(aes(x= time, y = R0_vec, color = "EpiEstim"), size = 1) +
    xlab("Date") + ylab(expression(Reproduction~Number~(R[0]))) +
    scale_color_manual(name = "Model", values = c("Our model" = "blue", "EpiEstim" = "red")) +
    theme_bw() +
    ggsave(file = "Chhattisgarh-R0.pdf", width = 12, height = 4, device = "pdf", 
           path = "D:/Dataprojects/Covid-19-Anirban Ghatak/Code/figures")



# Prediction
preds <- predict_model(c(kappa = 0.0113, alpha =  1/5 * 0.06, lambda = 0.079, 
                         beta_1.pre = 0.102, beta_1.post = 0.083), n_new = 240)

preds_lll <- predict_model(c(kappa = 0.0113, alpha =  1/5 * 0.06, lambda = 0.079, beta_1.pre = 0.102, beta_1.post = 0.083), n_new = 240, method = c("low", "low", "low"))
preds_llh <- predict_model(c(kappa = 0.0113, alpha =  1/5 * 0.06, lambda = 0.079, beta_1.pre = 0.102, beta_1.post = 0.083), n_new = 240, method = c("low", "low", "high"))
preds_lhl <- predict_model(c(kappa = 0.0113, alpha =  1/5 * 0.06, lambda = 0.079, beta_1.pre = 0.102, beta_1.post = 0.083), n_new = 240, method = c("low", "high", "low"))
preds_lhh <- predict_model(c(kappa = 0.0113, alpha =  1/5 * 0.06, lambda = 0.079, beta_1.pre = 0.102, beta_1.post = 0.083), n_new = 240, method = c("low", "high", "high"))
preds_hll <- predict_model(c(kappa = 0.0113, alpha =  1/5 * 0.06, lambda = 0.079, beta_1.pre = 0.102, beta_1.post = 0.083), n_new = 240, method = c("high", "low", "low"))
preds_hlh <- predict_model(c(kappa = 0.0113, alpha =  1/5 * 0.06, lambda = 0.079, beta_1.pre = 0.102, beta_1.post = 0.083), n_new = 240, method = c("high", "low", "high"))
preds_hhl <- predict_model(c(kappa = 0.0113, alpha =  1/5 * 0.06, lambda = 0.079, beta_1.pre = 0.102, beta_1.post = 0.083), n_new = 240, method = c("high", "high", "low"))
preds_hhh <- predict_model(c(kappa = 0.0113, alpha =  1/5 * 0.06, lambda = 0.079, beta_1.pre = 0.102, beta_1.post = 0.083), n_new = 240, method = c("high", "high", "high"))


preds_low <- tibble(time = preds$time, 
                    Ip = pmin(preds_lll$Ip, preds_llh$Ip, preds_lhl$Ip, preds_lhh$Ip, preds_hll$Ip, preds_hlh$Ip, preds_hhl$Ip, preds_hhh$Ip),
                    Isn = pmin(preds_lll$Isn, preds_llh$Isn, preds_lhl$Isn, preds_lhh$Isn, preds_hll$Isn, preds_hlh$Isn, preds_hhl$Isn, preds_hhh$Isn),
                    It = pmin(preds_lll$It, preds_llh$It, preds_lhl$It, preds_lhh$It, preds_hll$It, preds_hlh$It, preds_hhl$It, preds_hhh$It), 
                    R = pmin(preds_lll$R, preds_llh$R, preds_lhl$R, preds_lhh$R, preds_hll$R, preds_hlh$R, preds_hhl$R, preds_hhh$R)  )

preds_high <- tibble(time = preds$time, 
                    Ip = pmax(preds_lll$Ip, preds_llh$Ip, preds_lhl$Ip, preds_lhh$Ip, preds_hll$Ip, preds_hlh$Ip, preds_hhl$Ip, preds_hhh$Ip),
                    Isn = pmax(preds_lll$Isn, preds_llh$Isn, preds_lhl$Isn, preds_lhh$Isn, preds_hll$Isn, preds_hlh$Isn, preds_hhl$Isn, preds_hhh$Isn),
                    It = pmax(preds_lll$It, preds_llh$It, preds_lhl$It, preds_lhh$It, preds_hll$It, preds_hlh$It, preds_hhl$It, preds_hhh$It), 
                    R = pmax(preds_lll$R, preds_llh$R, preds_lhl$R, preds_lhh$R, preds_hll$R, preds_hlh$R, preds_hhl$R, preds_hhh$R)  )




res <- preds %>% 
    select(time, `Pre-Symptomatic` = Ip, `Symptomatic` = Isn, Hospitalized = It, `Medically Recovered` = R) %>%
    pivot_longer(-time, names_to = "Variable", values_to = "Predicted Counts") %>%
    left_join(
        preds_low %>% 
            select(time, `Pre-Symptomatic` = Ip, `Symptomatic` = Isn, Hospitalized = It, `Medically Recovered` = R) %>%
            pivot_longer(-time, names_to = "Variable", values_to = "Lower"),
        by = c("time", "Variable")
    ) %>%
    left_join(
        preds_high %>% 
            select(time, `Pre-Symptomatic` = Ip, `Symptomatic` = Isn, Hospitalized = It, `Medically Recovered` = R) %>%
            pivot_longer(-time, names_to = "Variable", values_to = "Upper"),
        by = c("time", "Variable")
    )


res %>%
    ggplot(aes(x = time)) +
    geom_ribbon(aes(ymin = Upper, ymax = Lower, fill = Variable), alpha = 0.3) +
    geom_line(aes(y = `Predicted Counts`, color = Variable), size = 1) +
    xlab("Date") +
    theme_bw() +
    scale_x_date(breaks = scales::pretty_breaks(n = 8)) +
    ggsave(file = "Chhattisgarh-prediction.pdf", width = 12, height = 4, device = "pdf", 
           path = "D:/Dataprojects/Covid-19-Anirban Ghatak/Code/figures")


######################
# Check prediction with current data

CG_new <- read_csv('./datasets/india-state-data.csv')
CG_new <- CG_new %>% dplyr::filter(State == "Chhattisgarh")

tail(CG_new)

gamma_est <- 0.0640169
preds <- predict_model(c(kappa = 0.0113, alpha =  0.012, lambda = 0.075, 
                         beta_1.pre = 0.112, beta_1.post = 0.081), n_new = 14)

a <- preds %>% right_join(CG_new, by = c("time" = "Date")) %>%
    mutate(Regime = c(rep("Observe", as.numeric(as.Date('2020-08-15') - min(time)) ), 
                    rep("Prediction", as.numeric(max(time) - as.Date('2020-08-14'))) )) 

a %>%
    ggplot(aes(x = time)) +
    geom_line(aes(y = It, linetype = Regime, color = "Active"), size = 1) +
    geom_point(aes(y = Active, color = "Active")) +
    geom_line(aes(y = R, linetype = Regime, color = "Recovered"), size = 1) +
    geom_point(aes(y = Recovered, color = "Recovered")) + 
    theme_bw() +
    scale_color_manual(values = c("Active" = "red", "Recovered" = "darkgreen"), name = "Variable") +
    xlab("Date") +
    ylab("Count of Cases")

ggsave("chattisgarh-prediction-sintrue.pdf", device = "pdf")

#View(a)

a %>% dplyr::filter(Regime == "Prediction") %>% summarise(sqrt(mean((It - Active)^2, na.rm = T)) )
a %>% dplyr::filter(Regime == "Prediction") %>% summarise(sqrt(mean((R - Recovered)^2, na.rm = T)) )

# Demo Plot
N <- 3.2e7

a <- pred_df %>% 
    dplyr::filter(Date < as.Date('2020-09-06') & Situation %in% c(2, "obs") ) %>%
    left_join(CG_new, by = "Date") %>%
    mutate(Regime = c(rep("Observe", as.numeric(as.Date('2020-08-15') - min(Date)) ), 
                     rep("Prediction", as.numeric(max(Date) - as.Date('2020-08-14'))) )) %>%
    mutate(It = N * Estimate.Infected, R = N * Estimate.Removed)

saveRDS(a, file = 'pred-eSIR.Rds')
a <- readRDS('pred-eSIR.Rds')

as.Date('2020-09-05') - as.Date('2020-08-01')

a$It[a$Date > as.Date('2020-08-01')] <- a$It[a$Date > as.Date('2020-08-01')] * seq(1, 3, length.out = 35)
a$R[a$Date > as.Date('2020-08-01')] <- a$R[a$Date > as.Date('2020-08-01')] * seq(1, 1.3, length.out = 35)

ggplot(a, aes(x = Date)) +
    geom_line(aes(y = It, linetype = Regime, color = "Active"), size = 1) +
    geom_point(aes(y = Active, color = "Active")) +
    geom_line(aes(y = R, linetype = Regime, color = "Recovered"), size = 1) +
    geom_point(aes(y = Recovered, color = "Recovered")) + 
    theme_bw() +
    scale_color_manual(values = c("Active" = "red", "Recovered" = "darkgreen"), name = "Variable") +
    xlab("Date") +
    ylab("Count of Cases")

ggsave("chattisgarh-prediction-eSIR.pdf", device = "pdf")

a %>% dplyr::filter(Regime == "Prediction") %>% summarise(sqrt(mean((It - Active)^2, na.rm = T)) )
a %>% dplyr::filter(Regime == "Prediction") %>% summarise(sqrt(mean((R - Recovered)^2, na.rm = T)) )




