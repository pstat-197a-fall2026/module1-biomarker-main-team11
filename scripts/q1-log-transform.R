## QUESTION 1: Why log-transform the protein levels?
## Hint from the assignment: look at the distribution of the raw values
## for a sample of proteins.
##
## Run from the repository root (same convention as preprocessing.R).
## Reads biomarker-raw.csv only; does not modify any data files.

library(tidyverse)

dir.create('results', showWarnings = FALSE)

# ---- 1. read raw data (same parsing as preprocessing.R, but NO transforms) ----

var_names <- read_csv('biomarker-raw.csv',
                      col_names = F,
                      n_max = 2,
                      col_select = -(1:2)) %>%
  t() %>%
  as_tibble() %>%
  rename(name = V1,
         abbreviation = V2) %>%
  na.omit()

biomarker_raw <- read_csv('biomarker-raw.csv',
                          skip = 2,
                          col_select = -2L,
                          col_names = c('group',
                                        'empty',
                                        pull(var_names, abbreviation),
                                        'ados'),
                          na = c('-', '')) %>%
  filter(!is.na(group)) %>%
  select(-any_of(c('empty', 'ados')))

# raw protein levels in long format
raw_long <- biomarker_raw %>%
  pivot_longer(-group, names_to = 'protein', values_to = 'level')

# sanity checks: log10 requires strictly positive values
cat('subjects:', nrow(biomarker_raw), '\n')
cat('proteins:', ncol(biomarker_raw) - 1, '\n')
cat('missing values:', sum(is.na(raw_long$level)), '\n')
cat('smallest raw value:', min(raw_long$level, na.rm = TRUE), '\n')
cat('non-positive values:', sum(raw_long$level <= 0, na.rm = TRUE), '\n\n')

# ---- 2. histograms for a random sample of proteins: raw vs log10 ----

set.seed(101422)
protein_sample <- sample(unique(raw_long$protein), size = 9)

sample_long <- raw_long %>%
  filter(protein %in% protein_sample) %>%
  mutate(log10_level = log10(level)) %>%
  pivot_longer(c(level, log10_level),
               names_to = 'scale',
               values_to = 'value') %>%
  mutate(scale = factor(scale,
                        levels = c('level', 'log10_level'),
                        labels = c('raw', 'log10')))

p_raw <- sample_long %>%
  filter(scale == 'raw') %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 30, fill = 'steelblue', color = 'white') +
  facet_wrap(~ protein, scales = 'free') +
  labs(title = 'Raw serum levels: right-skewed, a few very large values',
       x = 'raw level', y = 'number of subjects') +
  theme_minimal()

p_log <- sample_long %>%
  filter(scale == 'log10') %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 30, fill = 'darkorange', color = 'white') +
  facet_wrap(~ protein, scales = 'free') +
  labs(title = 'After log10: much closer to symmetric',
       x = 'log10(level)', y = 'number of subjects') +
  theme_minimal()

ggsave('results/q1-hist-raw.png', p_raw, width = 9, height = 7, dpi = 150)
ggsave('results/q1-hist-log10.png', p_log, width = 9, height = 7, dpi = 150)

print(p_raw)
print(p_log)

# ---- 3. quantify skewness across ALL proteins, raw vs log10 ----

skew <- function(x) {
  x <- x[!is.na(x)]
  mean((x - mean(x))^3) / sd(x)^3
}

skew_tbl <- raw_long %>%
  mutate(log10_level = log10(level)) %>%
  group_by(protein) %>%
  summarize(skew_raw = skew(level),
            skew_log10 = skew(log10_level),
            .groups = 'drop')

cat('Skewness across all proteins (0 = symmetric):\n')
skew_tbl %>%
  summarize(across(c(skew_raw, skew_log10),
                   list(median = median,
                        q25 = ~ quantile(.x, 0.25),
                        q75 = ~ quantile(.x, 0.75)))) %>%
  pivot_longer(everything(), names_to = c('scale', 'stat'),
               names_pattern = 'skew_(raw|log10)_(.*)') %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  print()

cat('\nproportion of proteins with skewness > 1:\n')
cat('  raw  :', round(mean(skew_tbl$skew_raw > 1), 3), '\n')
cat('  log10:', round(mean(skew_tbl$skew_log10 > 1), 3), '\n\n')

p_skew <- skew_tbl %>%
  pivot_longer(-protein, names_to = 'scale', values_to = 'skewness') %>%
  mutate(scale = recode(scale, skew_raw = 'raw', skew_log10 = 'log10')) %>%
  ggplot(aes(x = skewness, fill = scale)) +
  geom_histogram(bins = 50, position = 'identity', alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  labs(title = 'Skewness of each protein, before and after log10',
       x = 'sample skewness', y = 'number of proteins') +
  theme_minimal()

ggsave('results/q1-skewness.png', p_skew, width = 7, height = 4.5, dpi = 150)
print(p_skew)

# ---- 4. mean-variance relationship: does log10 stabilize the spread? ----
# If SD grows with the mean on the raw scale, proteins with high abundance
# would dominate any analysis that compares or pools them. Logging should
# remove that dependence.

mv_tbl <- raw_long %>%
  mutate(log10_level = log10(level)) %>%
  group_by(protein) %>%
  summarize(mean_raw = mean(level, na.rm = TRUE),
            sd_raw = sd(level, na.rm = TRUE),
            mean_log10 = mean(log10_level, na.rm = TRUE),
            sd_log10 = sd(log10_level, na.rm = TRUE),
            .groups = 'drop')

cat('Correlation between protein mean and protein SD:\n')
cat('  raw  :', round(cor(mv_tbl$mean_raw, mv_tbl$sd_raw), 3), '\n')
cat('  log10:', round(cor(mv_tbl$mean_log10, mv_tbl$sd_log10), 3), '\n\n')

p_mv <- bind_rows(
  mv_tbl %>% transmute(scale = 'raw', mean = mean_raw, sd = sd_raw),
  mv_tbl %>% transmute(scale = 'log10', mean = mean_log10, sd = sd_log10)
) %>%
  mutate(scale = factor(scale, levels = c('raw', 'log10'))) %>%
  ggplot(aes(x = mean, y = sd)) +
  geom_point(alpha = 0.4) +
  facet_wrap(~ scale, scales = 'free') +
  labs(title = 'Protein SD vs. protein mean',
       x = 'mean', y = 'standard deviation') +
  theme_minimal()

ggsave('results/q1-mean-sd.png', p_mv, width = 9, height = 4.5, dpi = 150)
print(p_mv)

# ---- 5. takeaway ----
# Raw serum levels are strictly positive, right-skewed, and have a spread that
# grows with the mean (typical of concentration/intensity measurements, where
# effects are multiplicative). A log transform pulls in the long right tail,
# makes each protein's distribution roughly symmetric, and makes the spread
# similar across proteins. That better matches the assumptions behind the
# t-test and keeps a handful of extreme raw values from dominating means and
# standard deviations.

# What the raw data show
# 
# All 1,317 proteins have strictly positive values (smallest is 12.9), so log10 is valid.
# In the sampled histograms, most proteins have a long right tail, with a few very large values (e.g., ER, HIF-1a, Marapsin).
# Median skewness across all proteins is 2.26 raw, versus 0.83 after log10.
# 76% of proteins have skewness above 1 raw, versus 45% after log10.
# 
# Spread vs. level
# 
# On the raw scale, a protein’s SD tracks its mean (correlation 0.78), so high-abundance proteins have much larger variance.
# After log10 that relationship nearly disappears (correlation 0.10).
# 
# Why it matters
# 
# The t-tests assume roughly symmetric, normal-like data, and a few extreme raw values can dominate means and SDs.
# Comparable spread across proteins keeps the later analyses from being driven by whichever proteins have the biggest raw numbers.
# Concentration measurements tend to behave multiplicatively, and logs turn that into additive effects.

# Caveat for your write-up: log10 improves things a lot but doesn’t fully symmetrize every protein, since 45% still have skewness above 1. That residual skew is part of why outliers show up in Q2.