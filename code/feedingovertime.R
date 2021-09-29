library(tidyverse)

baseline = read.csv("bloodmeal_baseline.csv")
int = read.csv("bloodmeal_intervention.csv")


baseline= baseline %>%
  mutate(prop_full = total_full / total_evaluated,
         prop_half =  total_half/ total_evaluated,
         prop_empty = total_empty / total_evaluated)


int= int %>%
  mutate(prop_full = total_full / total_evaluated,
         prop_half =  total_half/ total_evaluated,
         prop_empty = total_empty / total_evaluated)

baseline %>%
  pivot_longer(cols=prop_full:prop_empty, names_to = "quantity") %>%
  select(-contains("total")) %>%
  ggplot(aes(x=as.factor(Cluster), y=value, fill=quantity)) +
  geom_col()

int %>%
  pivot_longer(cols=prop_full:prop_empty, names_to = "quantity") %>%
  select(-contains("total")) %>%
  ggplot(aes(x=as.factor(Cluster), y=value, fill=quantity)) +
  geom_col()

