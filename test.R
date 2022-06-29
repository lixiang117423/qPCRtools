library(tidyverse)

iris %>%
  rstatix::t_test(Sepal.Length ~ Species, ref.group = "setosa") %>%
  dplyr::select(group2, p) %>%
  dplyr::mutate(signif = dplyr::case_when(
    p < 0.001 ~ "***",
    p > 0.001 & p < 0.01 ~ "**",
    p > 0.01 & p < 0.05 ~ "*",
    TRUE ~ "NS"
  )) %>%
  dplyr::add_row(group2 = )

