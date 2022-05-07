library(dplyr)
library(readr)
library(tidyr)
library(lubridate)

url1 <- 'https://api.covid19india.org/csv/latest/state_wise_daily.csv'
url2 <- 'https://api.covid19india.org/csv/latest/statewise_tested_numbers_data.csv'

df1 <- read_csv(url1)
df2 <- read_csv(url2)

stateList <- c(
    "an" = "Andaman and Nicobar Islands", "ap" = "Andhra Pradesh", "ar" = "Arunachal Pradesh", "as" = "Assam",
    "br" = "Bihar", "ch" = "Chandigarh", "ct" = "Chhattisgarh", "dl" = "Delhi", "dn" = "Dadra and Nagar Haveli and Daman and Diu",
    "ga" = "Goa", "gj" = "Gujarat", "hp" = "Himachal Pradesh", "hr" = "Haryana", "jh" = "Jharkhand", "jk" = "Jammu and Kashmir",
    "ka" = "Karnataka", "kl" = "Kerala", "la" = "Ladakh", "ld" = "Lakshadweep", "mh" = "Maharashtra", "ml" = "Meghalaya",
    "mn" = "Manipur", "mp" = "Madhya Pradesh", "mz" = "Mizoram", "nl" = "Nagaland", "or" = "Odisha", "pb" = "Punjab",
    "py" = "Puducherry", "rj" = "Rajasthan", "sk" = "Sikkim", "tg" = "Telangana", "tn" = "Tamil Nadu",
    "tr" = "Tripura", "up" = "Uttar Pradesh", "ut" = "Uttarakhand", "wb" = "West Bengal")

df1 <- df1 %>% 
    pivot_longer(-c(Date, Status), names_to = "State", values_to = "Value") %>%
    pivot_wider(names_from = "Status", values_from = "Value") %>%
    mutate(Date = dmy(Date), State = stateList[tolower(State)]) %>%
    select(Date, State, Confirmed, Recovered, Deceased) %>%
    drop_na() %>%
    group_by(State) %>%
    arrange(Date) %>%
    mutate(Confirmed = cumsum(Confirmed), Recovered = cumsum(Recovered), Deaths = cumsum(Deceased)) %>%
    mutate(Active = Confirmed - Recovered - Deaths) %>%
    select(Date, State, Confirmed, Active, Recovered, Deaths)

df2 <- df2 %>%
    select(Date = `Updated On`, State, Tested = `Total Tested`, Positive) %>%
    mutate(Date = dmy(Date), TPR = Positive / Tested)

statedf <- df1 %>% full_join(df2, by = c("Date", "State"))
statedf <- statedf %>% group_by(Date, State) %>% slice(1)  # Uttar Pradesh 4th May has a double entry, remove that

write_csv(statedf, path = './datasets/india-state-data.csv')


# Districtwise Data (get Migrant things)

url3 <- "https://api.covid19india.org/csv/latest/districts.csv"   # districtwise new data
distdf <- read_csv(url3)

distdf <- distdf %>% mutate(Active = Confirmed - Recovered - Deceased) %>%
    select(Date, State, District, Confirmed, Active, Recovered, Deaths = Deceased)    

write_csv(distdf, path = './datasets/india-district-level.csv')









