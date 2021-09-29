library(tidyverse)
library(lubridate)
load("~/Dropbox/DengueEntomologicalEffects/data/Iquitos_foi.RData")

foi_data = data.frame(date = dates.foi, 
                      foi = data.foi %>%
                        colSums()) %>%
  mutate(date = ymd(date)) %>%
  mutate(year = year(date)) %>%
  group_by(year) %>%
  summarize(yearly_foi = sum(foi))



# par(mfrow=c(1,2))
plot(colSums(data.foi)~dates.foi, xlab="Date", type="l",
     ylab="Force of Infection")



# infection data
load("~/Dropbox/DengueEntomologicalEffects/data/Iquitos_infections.RData")
data.inf
plot(colSums(data.inf)~dates.inf, xlab="Date", type="l",
     ylab="Infection Incidence")
