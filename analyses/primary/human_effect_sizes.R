library(tidyverse)


es_dir <- "data/human/derivatives/v3/547/effect_sizes/resolution_3.0/absolute/"
file <- file.path(es_dir, "effect_sizes.csv")
df_es_abs <- as_tibble(data.table::fread(file, header = TRUE))

mat_es_abs <- df_es_abs %>%
  select(-file) %>%
  as.matrix()

j <- 1:ncol(mat_es_abs)

j <- 1

get_stats <- function(j, x, stat = "t") {
  model <- summary(lm(x[, j] ~ 1))
  if (stat == "t") {
    out <- model$coefficients[1, "t value"]
  } else if (stat == "beta") {
    out <- model$coefficients[1, "Estimate"]
  }
  return(out)
}

betas <- map_dbl(
  .x = 1:ncol(mat_es_abs),
  .f = get_stats,
  x = mat_es_abs,
  stat = "beta"
)
tvals <- map_dbl(
  .x = 1:ncol(mat_es_abs),
  .f = get_stats,
  x = mat_es_abs,
  stat = "t"
)

df_stats <- tibble(beta = betas, tval = tvals) %>%
  mutate(
    pvals = 2 * (1 - pt(q = abs(tvals), df = nrow(mat_es_abs))),
    qvals = p.adjust(pvals, method = "fdr"),
    logq = -log10(qvals)
  )


ggplot(df_stats, aes(x = beta)) +
  geom_histogram()

pvals <- 2 * (1 - pt(q = abs(tvals), df = nrow(mat_es_abs)))
qvals <- p.adjust(pvals, method = "fdr")
df_qvals <- tibble(logq = -log10(qvals))
ggplot(df_qvals, aes(x = logq)) +
  geom_histogram() +
  geom_vline(xintercept = -log10(c(0.05, 0.1)))

qt(p = 0.9, df = 709)
dt(x = 2, df = 709)

rnorm()

qnorm(p = )
pnorm()
dnorm(x = 0)

?quantile()
