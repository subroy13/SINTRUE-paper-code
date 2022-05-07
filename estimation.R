library(dplyr)
library(readr)
library(tidyr)
library(EpiEstim)
setwd("D:/Academic/Projects/2020_Covid 19 Prediction Model_Anirban Ghatak_IIMV")

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
statedf <- read_csv('./datasets/india-state-data.csv')
statedf <- statedf %>% dplyr::filter(Date < max(Date))   # today's data is NA

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
    migrant_incoming_rate<-migrant_incoming_rate[1:12]
    
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
    
    # Estimation of c2 
    R0.all_states <- sapply(origin_states, FUN = estimate.R0)
    
    # need to add a lag of 3 days
    R0.all_states <- rbind(R0.all_states[1:3, ], R0.all_states[-( (nrow(R0.all_states)-2):nrow(R0.all_states) ), ])
    
    c2 <- rowSums(infected.origin.all_states * 
                  matrix(migrant_incoming_rate[origin_states], nrow = nrow(infected.origin.all_states), ncol = ncol(infected.origin.all_states), byrow = T ) *
                  R0.all_states )   # EQN OF c2
    
    c1 <- (1 - c2 - c3)  # EQN of c1
    
    
    # finally output everything
    return(list(c1 = c1, c2 = c2, c3 = c3))
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
#N_c <- 140

# Estimate of theta_t
theta_t <- diff(CGdf$Confirmed) * ( (14 / N_adj) + (avg_hh_size-2) / (N_total - N_adj))/2
plot(CGdf$Date[-1], theta_t, type = "l")


####################################
# ESTIMATION of I_p
I_p <- asym.incoming(CGdf$Date[-1]) * diff(CGdf$Confirmed) / theta_t
plot(CGdf$Date[-1], I_p, type = "l")


###################################
# ESTIMATION of Lambda

# Now consider estimate of P(tested | asymp & covid)
#p_symp <- 1 - (1 - 0.822) * (1 - 0.617) * (1 - 0.44)
p_symp <- 0.822 * 0.617 * 0.44
covid_prevalence <- 0.25

p_testAsymp <- (cumsum(asym.incoming(CGdf$Date[-1]) * diff(CGdf$Confirmed)) /N_adj) / ( (1 - p_symp) * covid_prevalence)
plot(CGdf$Date[-1], p_testAsymp, type = "l")

alpha_est <- 1/10
Lambda <- (theta_t / p_testAsymp) - theta_t - alpha_est

hist(Lambda[Lambda > 0], breaks = 25)
summary(Lambda[Lambda > 0])    # take the median
Lambda.est <- median(Lambda[Lambda > 0])


# consider, P(tested | symp & covid)
p_testsymp <- ( cumsum((1 - asym.incoming(CGdf$Date[-1])) * diff(CGdf$Confirmed)) /N_adj) /( p_symp * covid_prevalence)
plot(CGdf$Date[-1], p_testsymp, type = "l")


#########
# note that, delta_t / (delta_t + kappa) = P(tested | symp & covid)
# Estimation of delta_t, kappa, and I_sn

delta <- function(kappa){
    del <- (kappa * p_testsymp)/(1- p_testsymp)
    return(del)
}

# I_sn as a function of kappa
I_sn <- function(t, kappa){
    del <- delta(kappa) 
    int.factor <- approxfun(x = 0:t, c(0, cumsum(del + kappa)[1:t]) )
    
    I_p_approx <- approxfun(x = c(0, which(!is.na(I_p)) ), y = c(50, I_p[which(!is.na(I_p))] ))
    int <- integrate(f = function(x) { exp(int.factor(x)) * I_p_approx(x) } , lower = 0, upper=t, subdivisions = 2000, stop.on.error=FALSE)$value
    return( exp(-int.factor(t)) * int * alpha_est)
}

# delta_t * I_sn as a function of kappa
str <- function(kappa){
    del <- delta(kappa)
    Isn <- sapply(1:length(del), I_sn, kappa = kappa)
    return(del * Isn)
}

LS.kappa <- function(kappa){
    pred <- str(kappa)
    y <- (1 - asym.incoming(CGdf$Date[-1])) * diff(CGdf$Confirmed)
    out <- sum((pred - y)^2)
    return(log(out))
}

init<- 0.5
kappa<-optim(init,fn=LS.kappa,lower=0, upper=1,control = list(trace = TRUE), method = "L-BFGS-B")$par
kappa 

plot((1 - asym.incoming(CGdf$Date[-1])) * diff(CGdf$Confirmed), type = "l")
points(str(kappa), type = "l", col = "red")

plot(CGdf$Date[-1], sapply(1:(nrow(CGdf)-1), FUN = I_sn, kappa = kappa), type = "l")
plot(CGdf$Date[-1], delta(kappa), type = "l")

Isn <- sapply(1:(nrow(CGdf)-1), FUN = I_sn, kappa = kappa)

####################
# Estimate of Tau

# function to estimate Tau
est.Tau <- function(It, E) {
    dE <- numDeriv(E)
    fit <- MASS::rlm(dE ~ It - 1, maxit = 100)
    print(summary(fit))
    tau <- fit$coefficients
    names(tau) <- "Tau"
    return(tau)
}

est.Tau(It = CGdf$Active, E = CGdf$Deaths)

########################
# Estimation of beta_1t, beta_2t
# We know, dI_p/dt + (alpha + lambda + theta_t) I_p = beta_1t * S * I_p + beta_2t * S * I_sn
# And, R0_t = beta_1t/(lambda + alpha + theta_t) * (I_p / I) + beta_2t/(kappa + delta_t) * (Isn / I)

est.Beta <- function(state_name) {
    S <- N_total - c(I_p, NA) - c(Isn, NA) - CGdf$Confirmed + 15000 * mig$c1
    diffI_p <- c(numDeriv(I_p), NA)
    R0_t <- estimate.R0(state_name)
    del <- delta(kappa)
    I.total <- I_p + Isn + CGdf$Active
    
    solve.beta <- function(t) {
        tryCatch({
            const <- c(diffI_p[t] + (alpha_est + Lambda.est + theta_t[t]) * I_p[t], R0_t[t])
            coeff <- matrix(c(S[t] * I_p[t], S[t] * Isn[t],  
                              I_p[t] / (I.total[t] * (Lambda.est + alpha_est + theta_t[t]) ),
                              Isn[t] / (I.total[t] * (kappa + del[t]) ) ), nrow = 2, ncol = 2, byrow = TRUE)
            betas <- solve(coeff, const)
            return(betas)
        },
        error = function(cond) {
            return(c(NA, NA))   
        })
    }
    
    betas <- sapply(1:nrow(CGdf), FUN = solve.beta)
    return(betas)
}

x <- est.Beta("Chhattisgarh")

########################
# Estimation of beta_1t, beta_2t   (TYPE 2)
# R0_t = beta_1t/(lambda + alpha + theta_t) * (I_p / I) + beta_2t/(kappa + delta_t) * (Isn / I)
# and do 0 = beta1 *(1-diff(active)/active) - beta2

est.Beta <- function(state_name) {
    S <- N_total - c(I_p, NA) - c(Isn, NA) - CGdf$Confirmed + 15000 * mig$c1
    R0_t <- estimate.R0(state_name)
    del <- delta(kappa)
    I.total <- I_p + Isn + CGdf$Active
    
    solve.beta <- function(t) {
        tryCatch({
            const <- c( R0_t[t], 0)
            coeff <- matrix(
                c( I_p[t] / (I.total[t] * (Lambda.est + alpha_est + theta_t[t])),
                   Isn[t] / (I.total[t] * (kappa + del[t]) ),
                   1-  (abs(CGdf$Active[t]-CGdf$Active[t-1]))/CGdf$Active[t], -1),
                nrow = 2, ncol = 2, byrow = TRUE)
            betas <- solve(coeff, const)
            return(betas)
        },
        error = function(cond) {
            return(c(NA, NA))   
        })
    }
    
    betas <- sapply(2:nrow(CGdf), FUN = solve.beta)
    return(betas)
}

x <- est.Beta("Chhattisgarh")

beta1<-x[1,]
beta2<- x[2,]
summary(beta1[beta1>0])

summary(beta2[beta2>0])

which(beta1<0)



##################################
# Estimation of Gamma and f

derivS <- numDeriv(S)

derivR <- numDeriv(CGdf$Recovered)
## in_rate<- 15000
## out_rate<-100 (to be estimated)

f <- (derivS+ beta1_t * S * Ip/N_total + beta2_t * S * Isn/N_total - Lambda.est*I_p - kappa * Isn 
      + 100 - 15000 * mig$c1)/(derivR)

summary(f)


# function to estimate Tau
est.Gamma <- function(It, R, f = 0) {
    dR <- numDeriv(R)
    fit <- MASS::rlm(dR ~ It - 1, maxit = 100)
    print(summary(fit))
    gam <- fit$coefficients * (1 + f)
    names(gam) <- "Gamma"
    return(gam)
}

est.Gamma(It = CGdf$Active, R = CGdf$Recovered)





