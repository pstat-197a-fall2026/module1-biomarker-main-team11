#combined script for questions 1-4. 
#all used libraries loaded here
library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)
library(glmnet)


#----Question 1----


dir.create('results', showWarnings = FALSE)

# ---- 1. read raw data (same parsing as preprocessing.R, but NO transforms) ----

var_names <- read_csv('data/biomarker-raw.csv',
                      col_names = F,
                      n_max = 2,
                      col_select = -(1:2)) %>%
  t() %>%
  as_tibble() %>%
  rename(name = V1,
         abbreviation = V2) %>%
  na.omit()

biomarker_raw <- read_csv('data/biomarker-raw.csv',
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


#----Question 2----


dir.create('results', showWarnings = FALSE)

# ---- 1. preprocessing WITHOUT trimming ----

var_names <- read_csv('data/biomarker-raw.csv',
                      col_names = F,
                      n_max = 2,
                      col_select = -(1:2)) %>%
  t() %>%
  as_tibble() %>%
  rename(name = V1,
         abbreviation = V2) %>%
  na.omit()

# same as preprocessing.R (log10, center, scale) but with the trim() step removed
biomarker_untrimmed <- read_csv('data/biomarker-raw.csv',
                                skip = 2,
                                col_select = -2L,
                                col_names = c('group',
                                              'empty',
                                              pull(var_names, abbreviation),
                                              'ados'),
                                na = c('-', '')) %>%
  filter(!is.na(group)) %>%
  mutate(across(.cols = -c(group, ados),
                ~ scale(log10(.x))[, 1])) %>%
  # subject id = row position in the data file (after dropping the blank row)
  mutate(subject = row_number()) %>%
  select(subject, group, ados, everything())

# the cutoff used by the trimming step in preprocessing.R
cutoff <- 3

n_subjects <- nrow(biomarker_untrimmed)
n_proteins <- ncol(biomarker_untrimmed) - 3
group_sizes <- biomarker_untrimmed %>% count(group)

cat('subjects:', n_subjects, ' proteins:', n_proteins, '\n')
print(group_sizes)

# long format of the standardized (untrimmed) values
long <- biomarker_untrimmed %>%
  select(-ados) %>%
  pivot_longer(-c(subject, group),
               names_to = 'protein',
               values_to = 'z')

# ---- 2. how many values are outlying overall? ----

long <- long %>%
  mutate(outlier = abs(z) > cutoff,
         direction = case_when(z > cutoff ~ 'high',
                               z < -cutoff ~ 'low',
                               TRUE ~ 'none'))

total_outliers <- sum(long$outlier)
cat('\nTotal outlying values (|z| >', cutoff, '):', total_outliers,
    'out of', nrow(long),
    sprintf('(%.2f%%)\n', 100 * total_outliers / nrow(long)))

# a normal distribution would give this many per subject, on average
expected_per_subject <- n_proteins * 2 * pnorm(-cutoff)
cat('Expected outliers per subject if values were normal:',
    round(expected_per_subject, 2), '\n')

cat('\nOutliers by direction:\n')
long %>% filter(outlier) %>% count(direction) %>% print()

# ---- 3. outliers per SUBJECT ----

subject_outliers <- long %>%
  group_by(subject, group) %>%
  summarize(n_outliers = sum(outlier),
            n_high = sum(direction == 'high'),
            n_low = sum(direction == 'low'),
            max_abs_z = max(abs(z)),
            .groups = 'drop') %>%
  arrange(desc(n_outliers))

cat('\nDistribution of outlying values per subject:\n')
print(summary(subject_outliers$n_outliers))

cat('\nTop 15 subjects by number of outlying values:\n')
print(subject_outliers %>% slice_head(n = 15), n = 15)

# tabulation requested in the hint: number of subjects at each outlier count
cat('\nTabulation: number of subjects by count of outlying values (binned):\n')
subject_outliers %>%
  mutate(bin = cut(n_outliers,
                   breaks = c(-Inf, 0, 2, 5, 10, 20, 50, Inf),
                   labels = c('0', '1-2', '3-5', '6-10',
                              '11-20', '21-50', '>50'))) %>%
  count(bin, group) %>%
  pivot_wider(names_from = group, values_from = n, values_fill = 0) %>%
  print()

# flag "outlier subjects": counts above the Tukey upper fence for the counts
q <- quantile(subject_outliers$n_outliers, c(0.25, 0.75))
fence <- q[2] + 1.5 * (q[2] - q[1])
subject_outliers <- subject_outliers %>%
  mutate(outlier_subject = n_outliers > fence)

cat('\nTukey upper fence on per-subject counts:', fence, '\n')
cat('Subjects flagged as outlier subjects:',
    sum(subject_outliers$outlier_subject), '\n')
print(subject_outliers %>% filter(outlier_subject))

# ---- 4. is one group more outlier-prone? ----

cat('\nOutlying values per subject, by group:\n')
group_summary <- subject_outliers %>%
  group_by(group) %>%
  summarize(n_subjects = n(),
            total_outliers = sum(n_outliers),
            mean_per_subject = mean(n_outliers),
            median_per_subject = median(n_outliers),
            max_per_subject = max(n_outliers),
            n_high = sum(n_high),
            n_low = sum(n_low),
            subjects_flagged = sum(outlier_subject),
            .groups = 'drop')
print(group_summary)

# rank-based comparison of per-subject counts (counts are very skewed)
cat('\nWilcoxon rank-sum test of per-subject counts, ASD vs TD:\n')
print(wilcox.test(n_outliers ~ group, data = subject_outliers,
                  exact = FALSE))

# permutation test on the difference in mean counts
set.seed(101422)
obs_diff <- with(subject_outliers,
                 mean(n_outliers[group == 'ASD']) -
                   mean(n_outliers[group == 'TD']))
perm_diffs <- replicate(10000, {
  g <- sample(subject_outliers$group)
  mean(subject_outliers$n_outliers[g == 'ASD']) -
    mean(subject_outliers$n_outliers[g == 'TD'])
})
cat('\nPermutation test, difference in mean outlier count (ASD - TD):\n')
cat('  observed difference:', round(obs_diff, 2), '\n')
cat('  two-sided p-value  :', mean(abs(perm_diffs) >= abs(obs_diff)), '\n')

# are flagged subjects concentrated in one group?
cat('\nFlagged outlier subjects by group (Fisher exact test):\n')
flag_tab <- table(subject_outliers$group, subject_outliers$outlier_subject)
print(flag_tab)
if (ncol(flag_tab) == 2) print(fisher.test(flag_tab))

# do the groups differ in whether a subject has ANY outlying value?
cat('\nSubjects with at least one outlying value, by group:\n')
any_tab <- table(subject_outliers$group, subject_outliers$n_outliers > 0)
print(any_tab)

# ---- 5. do outliers cluster in particular proteins instead? ----

protein_outliers <- long %>%
  group_by(protein) %>%
  summarize(n_outliers = sum(outlier), .groups = 'drop') %>%
  arrange(desc(n_outliers))

cat('\nProteins with the most outlying values:\n')
print(protein_outliers %>% slice_head(n = 10))

cat('\nShare of all outlying values coming from the flagged subjects:',
    round(sum(subject_outliers$n_outliers[subject_outliers$outlier_subject]) /
            total_outliers, 3), '\n')

# ---- 6. plots ----

# (a) outlier count per subject, ordered, colored by group
p_subject <- subject_outliers %>%
  arrange(desc(n_outliers)) %>%
  mutate(rank = row_number()) %>%
  ggplot(aes(x = rank, y = n_outliers, fill = group)) +
  geom_col() +
  geom_hline(yintercept = fence, linetype = 'dashed') +
  labs(title = 'Outlying values per subject (|z| > 3, no trimming)',
       subtitle = 'Dashed line: Tukey upper fence on per-subject counts',
       x = 'subjects, ordered by outlier count',
       y = 'number of outlying values') +
  theme_minimal()

# (b) distribution of per-subject counts by group
p_group <- subject_outliers %>%
  ggplot(aes(x = group, y = n_outliers, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  labs(title = 'Outlying values per subject, by group',
       x = NULL, y = 'number of outlying values') +
  theme_minimal() +
  theme(legend.position = 'none')

# (c) heatmap-style view for the most extreme subjects: which proteins?
top_subjects <- subject_outliers %>% slice_head(n = 5) %>% pull(subject)

p_top <- long %>%
  filter(subject %in% top_subjects) %>%
  ggplot(aes(x = z, fill = group)) +
  geom_histogram(bins = 60) +
  geom_vline(xintercept = c(-cutoff, cutoff), linetype = 'dashed') +
  facet_wrap(~ paste('subject', subject, '-', group), scales = 'free_y') +
  labs(title = 'Standardized values for the 5 subjects with the most outliers',
       x = 'z-score (untrimmed)', y = 'number of proteins') +
  theme_minimal()

ggsave('results/q2-outliers-by-subject.png', p_subject,
       width = 9, height = 5, dpi = 150)
ggsave('results/q2-outliers-by-group.png', p_group,
       width = 5, height = 5, dpi = 150)
ggsave('results/q2-top-subjects.png', p_top,
       width = 9, height = 6, dpi = 150)

print(p_subject)
print(p_group)
print(p_top)

# ---- 7. save the summary tables for the report ----

write_csv(subject_outliers, 'results/q2-subject-outlier-counts.csv')
write_csv(group_summary, 'results/q2-group-summary.csv')


#----Question 3-a----


load('data/biomarker-clean.RData')
set.seed(197)

# Using 80% for training and 20% for testing (matching in-class ratio)
biomarker_split <- initial_split(biomarker_clean, prop = 0.8, strata = group)
train_data <- training(biomarker_split)
test_data  <- testing(biomarker_split)

cat("\n--- Data Partition Complete ---\n")
cat("Training Set Count:", nrow(train_data), "\n")
cat("Testing Set Count :", nrow(test_data), "\n")


# STEP 1: MULTIPLE TESTING (STRICTLY ON train_data)
test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = FALSE)
}

ttests_out <- train_data %>%
  select(-ados) %>%
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  nest(data = c(level, group)) %>% 
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  arrange(p_value) %>%
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m*hm*p_value/rank)

# Select top 10 significant proteins from training set
proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 10) %>%
  pull(protein)



# STEP 2: RANDOM FOREST (STRICTLY ON train_data)
rf_predictors <- train_data %>%
  select(-c(group, ados))

rf_response <- train_data %>% pull(group) %>% factor()

set.seed(101422)
rf_out <- randomForest(x = rf_predictors, 
                       y = rf_response, 
                       ntree = 1000, 
                       importance = TRUE)

# Compute top 10 important proteins from training set
proteins_s2 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)



# STEP 3: FEATURE INTERSECTION & LOGISTIC REGRESSION FITTING
# Combine selected features from training set only
proteins_sstar <- intersect(proteins_s1, proteins_s2)

cat("\n--- Selected Proteins (Training Set Intersection) ---\n")
print(proteins_sstar)

train_sstar <- train_data %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(is_asd = (group == 'ASD')) %>%
  select(-group)

# Fit logistic regression model ONLY on training data
fit_train <- glm(is_asd ~ ., 
                 data = train_sstar, 
                 family = 'binomial')



# STEP 4: EVALUATE ACCURACY ON HELD-OUT TEST DATA
# Prepare test subset using the selected proteins
test_sstar <- test_data %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(is_asd = (group == 'ASD')) %>%
  select(-group)


# Calculate predictions on the held-out test set
test_results <- test_sstar %>%
  add_predictions(fit_train, type = "response") %>%
  mutate(
    truth_fac = factor(is_asd, levels = c(FALSE, TRUE)),
    pred_class = factor(pred > 0.5, levels = c(FALSE, TRUE))
  )

# Classification metrics
classification_results <- yardstick::metric_set(
  yardstick::sensitivity,
  yardstick::specificity,
  yardstick::accuracy
)(
  data = test_results,
  truth = truth_fac,
  estimate = pred_class,
  event_level = "second"
)

# ROC AUC
auc_result <- yardstick::roc_auc(
  data = test_results,
  truth = truth_fac,
  !!!rlang::syms("pred"),
  event_level = "second"
)

cat("\n--- Classification Metrics on Test Set ---\n")
print(classification_results)

cat("\n--- ROC AUC on Test Set ---\n")
print(auc_result)

cat("\n--- Number of Test Predictions ---\n")
print(nrow(test_results))


#----Question 3-b----


load("data/biomarker-clean.RData")


#splitting data before feature selection

set.seed(101422)

biomarker_split <- initial_split(
  biomarker_clean,
  prop = 0.8,
  strata = group
)

train_df <- training(biomarker_split)
test_df  <- testing(biomarker_split)

#function to select n proteins

select_proteins <- function(train_df, n) {
  
  #t test selection
  
  test_fn <- function(.df) {
    t_test(
      .df,
      formula = level ~ group,
      order = c("ASD", "TD"),
      alternative = "two-sided",
      var.equal = FALSE
    )
  }
  
  ttests_out <- train_df %>%
    select(-ados) %>%
    pivot_longer(
      -group,
      names_to = "protein",
      values_to = "level"
    ) %>%
    nest(data = c(level, group)) %>%
    mutate(ttest = map(data, test_fn)) %>%
    unnest(ttest) %>%
    arrange(p_value) %>%
    mutate(
      m = n(),
      hm = log(m) + 1 / (2 * m) - digamma(1),
      rank = row_number(),
      p.adj = m * hm * p_value / rank
    )
  
  proteins_s1 <- ttests_out %>%
    slice_min(p.adj, n = n, with_ties = FALSE) %>%
    pull(protein)
  
  #random forest select
  predictors <- train_df %>%
    select(-c(group, ados))
  
  response <- factor(train_df$group)
  
  set.seed(101422)
  
  rf_out <- randomForest(
    x = predictors,
    y = response,
    ntree = 1000,
    importance = TRUE
  )
  
  proteins_s2 <- rf_out$importance %>%
    as_tibble() %>%
    mutate(protein = rownames(rf_out$importance)) %>%
    slice_max(
      MeanDecreaseGini,
      n = n,
      with_ties = FALSE
    ) %>%
    pull(protein)
  #intersect both methods
  
  selected_proteins <- intersect(
    proteins_s1,
    proteins_s2
  )
  
  # Return the selected proteins and both individual lists
  list(
    selected = selected_proteins,
    ttest = proteins_s1,
    random_forest = proteins_s2
  )
}

#fitting/eval classifier
evaluate_model <- function(train_df, test_df, proteins) {
  
  # Build model datasets using only selected proteins
  train_model <- train_df %>%
    mutate(
      class = factor(
        group == "ASD",
        levels = c(FALSE, TRUE),
        labels = c("TD", "ASD")
      )
    ) %>%
    select(class, all_of(proteins))
  
  test_model <- test_df %>%
    mutate(
      class = factor(
        group == "ASD",
        levels = c(FALSE, TRUE),
        labels = c("TD", "ASD")
      )
    ) %>%
    select(class, all_of(proteins))
  
  # Fit logistic regression
  fit <- glm(
    class ~ .,
    data = train_model,
    family = binomial()
  )
  
  # Predict probabilities on the held-out test set
  predictions <- test_model %>%
    mutate(
      pred = predict(
        fit,
        newdata = test_model,
        type = "response"
      ),
      pred_class = factor(
        if_else(pred > 0.5, "ASD", "TD"),
        levels = c("TD", "ASD")
      )
    )
  
  # Calculate classification metrics
  class_metrics <- metric_set(
    sensitivity,
    specificity,
    accuracy,
    roc_auc
  )
  
  metrics <- predictions %>%
    class_metrics(
      truth = class,
      estimate = pred_class,
      pred,
      event_level = "second"
    )
  
  list(
    model = fit,
    predictions = predictions,
    metrics = metrics
  )
}

#top 20 vs 10

# Select proteins using training data only
selection_10 <- select_proteins(train_df, 10)
selection_20 <- select_proteins(train_df, 20)

# Fit and evaluate each model on the SAME test set
results_10 <- evaluate_model(
  train_df,
  test_df,
  selection_10$selected
)

results_20 <- evaluate_model(
  train_df,
  test_df,
  selection_20$selected
)

#print tests

cat("\nTOP-10 ANALYSIS\n")
cat("Number selected:", length(selection_10$selected), "\n")
print(selection_10$selected)

cat("\nTOP-20 ANALYSIS\n")
cat("Number selected:", length(selection_20$selected), "\n")
print(selection_20$selected)

#comparing

comparison <- bind_rows(
  results_10$metrics %>%
    mutate(approach = "Top 10"),
  
  results_20$metrics %>%
    mutate(approach = "Top 20")
) %>%
  select(approach, .metric, .estimate) %>%
  pivot_wider(
    names_from = .metric,
    values_from = .estimate
  ) %>%
  select(
    approach,
    n_selected = everything(),
    sensitivity,
    specificity,
    accuracy,
    roc_auc
  )

# Add the number of proteins selected
comparison <- bind_rows(
  results_10$metrics %>% mutate(approach = "Top 10"),
  results_20$metrics %>% mutate(approach = "Top 20")
) %>%
  select(approach, .metric, .estimate) %>%
  pivot_wider(names_from = .metric, values_from = .estimate) %>%
  mutate(
    n_selected = c(
      length(selection_10$selected),
      length(selection_20$selected)
    ),
    .before = sensitivity
  ) %>%
  select(approach, n_selected, sensitivity, specificity, accuracy, roc_auc)

print(comparison)


#----Question 3-c----


dir.create('results', showWarnings = FALSE)

select_on_train_only <- FALSE

top_n_values <- c(10, 20, 30)
seeds <- c(101422, 1:9)

# ---- data -------------------------------------------------------------------

load("data/biomarker-clean.RData")

stopifnot(exists('biomarker_clean'))

biomarker_clean <- biomarker_clean %>%
  select(-any_of('ados'))


ttest_ranks <- function(df) {
  df %>%
    pivot_longer(-group, names_to = 'protein', values_to = 'level') %>%
    group_by(protein) %>%
    summarize(p_value = t.test(level ~ group, var.equal = FALSE)$p.value,
              .groups = 'drop') %>%
    arrange(p_value) %>%
    mutate(m = n(),
           hm = log(m) + 1 / (2 * m) - digamma(1),
           step = row_number(),
           p_adj = m * hm * p_value / step) %>%
    arrange(p_adj) %>%
    mutate(sel_rank = row_number()) %>%
    select(protein, sel_rank)
}

rf_ranks <- function(df) {
  set.seed(101422)
  rf <- randomForest(x = select(df, -group),
                     y = factor(df$group),
                     ntree = 1000,
                     importance = TRUE)
  tibble(protein = rownames(rf$importance),
         gini = rf$importance[, 'MeanDecreaseGini']) %>%
    arrange(desc(gini)) %>%
    mutate(sel_rank = row_number()) %>%
    select(protein, sel_rank)
}

lasso_ranks <- function(df) {
  x <- as.matrix(select(df, -group))
  y <- factor(df$group)
  fit <- glmnet(x, y, family = 'binomial', alpha = 1)
  beta <- as.matrix(fit$beta)
  entry <- apply(beta != 0, 1, function(r) if (any(r)) which(r)[1] else Inf)
  size <- abs(beta[, ncol(beta)])
  tibble(protein = rownames(beta), entry = entry, size = size) %>%
    arrange(entry, desc(size)) %>%
    mutate(sel_rank = row_number()) %>%
    select(protein, sel_rank)
}

get_ranks <- function(df) {
  list(ttest = ttest_ranks(df),
       rf = rf_ranks(df),
       lasso = lasso_ranks(df))
}

top_set <- function(ranks, n) {
  ranks %>% filter(sel_rank <= n) %>% pull(protein)
}

fuzzy_intersect <- function(sets, k) {
  counts <- table(unlist(sets))
  if (length(counts) == 0) return(character(0))
  names(counts)[as.vector(counts) >= k]
}

class_metrics <- metric_set(sensitivity, specificity, accuracy, roc_auc)

evaluate <- function(proteins, train, test) {
  if (length(proteins) == 0) {
    return(tibble(n_proteins = 0L, sensitivity = NA_real_,
                  specificity = NA_real_, accuracy = NA_real_,
                  roc_auc = NA_real_, separated = NA))
  }
  
  prep <- function(d) {
    d %>%
      select(group, any_of(proteins)) %>%
      mutate(class = group == 'ASD') %>%
      select(-group)
  }
  tr <- prep(train)
  te <- prep(test)
  
  fit <- suppressWarnings(glm(class ~ ., data = tr, family = 'binomial'))
  
  te_pred <- te %>%
    mutate(pred = suppressWarnings(predict(fit, newdata = te, type = 'response')),
           est = factor(pred > 0.5, levels = c(FALSE, TRUE)),
           truth = factor(class, levels = c(FALSE, TRUE)))
  
  class_metrics(te_pred, truth = truth, pred, estimate = est,
                event_level = 'second') %>%
    select(.metric, .estimate) %>%
    pivot_wider(names_from = .metric, values_from = .estimate) %>%
    mutate(n_proteins = length(proteins),
           separated = any(fit$fitted.values < 1e-6 | fit$fitted.values > 1 - 1e-6))
}

ranks_all <- if (!select_on_train_only) get_ranks(biomarker_clean) else NULL

run_once <- function(seed) {
  set.seed(seed)
  split <- initial_split(biomarker_clean, prop = 0.8)
  train <- training(split)
  test <- testing(split)
  
  ranks <- if (select_on_train_only) get_ranks(train) else ranks_all
  
  baseline <- tibble(
    rule = 'baseline: t-test AND RF',
    top_n = 10,
    proteins = list(intersect(top_set(ranks$ttest, 10), top_set(ranks$rf, 10)))
  )
  
  three_method <- expand_grid(top_n = top_n_values, k = c(3, 2)) %>%
    mutate(rule = if_else(k == 3, 'hard: all 3 methods', 'fuzzy: at least 2 of 3'),
           proteins = map2(top_n, k, function(n, k) {
             fuzzy_intersect(list(top_set(ranks$ttest, n),
                                  top_set(ranks$rf, n),
                                  top_set(ranks$lasso, n)), k)
           })) %>%
    select(rule, top_n, proteins)
  
  bind_rows(baseline, three_method) %>%
    mutate(perf = map(proteins, evaluate, train = train, test = test),
           seed = seed) %>%
    unnest(perf)
}

results <- map_dfr(seeds, run_once)

# ---- results ----------------------------------------------------------------

cat('\n=== Baseline split (seed 101422): which proteins were kept ===\n')
results %>%
  filter(seed == 101422) %>%
  mutate(proteins = map_chr(proteins, ~ paste(.x, collapse = ', '))) %>%
  select(rule, top_n, n_proteins, accuracy, roc_auc, proteins) %>%
  print(n = Inf, width = Inf)

cat('\n=== Averaged over', length(seeds), 'train/test splits ===\n')
summary_tbl <- results %>%
  select(-proteins) %>%
  group_by(rule, top_n) %>%
  summarize(n_empty = sum(n_proteins == 0),
            n_proteins = mean(n_proteins),
            across(c(accuracy, sensitivity, specificity, roc_auc),
                   ~ mean(.x, na.rm = TRUE)),
            .groups = 'drop')
print(summary_tbl, n = Inf, width = Inf)

write_csv(summary_tbl, 'results/q3c-summary.csv')

p <- summary_tbl %>%
  ggplot(aes(x = factor(top_n), y = accuracy, fill = rule)) +
  geom_col(position = 'dodge') +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = 'Test accuracy: hard vs fuzzy intersection',
       x = 'size of each method\'s top-protein set', y = 'mean test accuracy') +
  theme_minimal()

ggsave('results/q3c-accuracy.png', p, width = 8, height = 4.5, dpi = 150)
print(p)


#----Question 4----


load('data/biomarker-clean.RData')


# Prepare data 
biomarker_all <- biomarker_clean %>%
  select(-ados) %>%
  mutate(class = factor(group == 'ASD', levels = c(FALSE, TRUE), labels = c("TD", "ASD"))) %>%
  select(-group)

# . Same seed/split as the test script. Selection only uses the training set
set.seed(101422)
biomarker_split <- initial_split(biomarker_all, prop = 0.8)
train_data <- training(biomarker_split)

# 3. Rank proteins
rank_proteins <- function(df) {
  x <- as.matrix(select(df, -class))
  fit <- glmnet(x, df$class, family = "binomial", alpha = 1, nlambda = 200)
  beta <- as.matrix(fit$beta)                                   # proteins x lambdas
  entry <- apply(beta != 0, 1, function(v) if (any(v)) which(v)[1] else Inf)
  ord <- order(entry, -abs(beta[, ncol(beta)]))                 # ties: bigger final coef first
  names(entry)[ord][is.finite(entry[ord])]
}

# 4. Final panel: top 5 proteins from the lasso ranking on the training set
panel_size <- 5
proteins_lasso <- head(rank_proteins(train_data), panel_size)

cat("\nTop", panel_size, "proteins. Copy into my_proteins:\n\n")
cat("my_proteins <- c(", paste0('"', proteins_lasso, '"', collapse = ", "), ")\n", sep = "")

# ---------------- PLOT: LASSO PATH AND SELECTED PANEL ----------------
x_train <- as.matrix(select(train_data, -class))
lasso_fit <- glmnet(x_train, train_data$class, family = "binomial", alpha = 1, nlambda = 200)

# Standardize coefficients (coef x protein SD) so paths are comparable across proteins
beta_std <- as.matrix(lasso_fit$beta) * apply(x_train, 2, sd)

path_df <- as.data.frame(beta_std) %>%
  rownames_to_column("protein") %>%
  pivot_longer(-protein, names_to = "step", values_to = "coef") %>%
  mutate(step = as.integer(str_remove(step, "^s")),
         log_lambda = log(lasso_fit$lambda[step + 1]),
         in_panel = protein %in% proteins_lasso)

# Penalty at which the last of the 5 panel proteins enters the model
entry_step <- apply(as.matrix(lasso_fit$beta)[proteins_lasso, , drop = FALSE] != 0, 1,
                    function(v) which(v)[1])
panel_log_lambda <- log(lasso_fit$lambda[max(entry_step)])

lasso_path_plot <- ggplot() +
  # all other proteins in the background
  geom_line(data = filter(path_df, !in_panel),
            aes(log_lambda, coef, group = protein),
            colour = "grey80", linewidth = 0.3) +
  # the 5 selected proteins on top
  geom_line(data = filter(path_df, in_panel),
            aes(log_lambda, coef, colour = protein), linewidth = 1) +
  geom_vline(xintercept = panel_log_lambda, linetype = "dashed") +
  annotate("text", x = panel_log_lambda, y = Inf, vjust = 1.5, hjust = -0.05,
           label = paste0("all ", panel_size, " panel proteins in model"), size = 3.5) +
  scale_x_reverse() +   # strong penalty (few proteins) on the left, weak on the right
  labs(title = "Lasso coefficient paths on the training set",
       subtitle = paste0("Top ", panel_size, " proteins by entry order highlighted; grey = all other proteins"),
       x = "log(penalty)",
       y = "Standardized coefficient",
       colour = "Panel protein") +
  theme_minimal()

print(lasso_path_plot)
ggsave("results/report_files/lasso_path.png", lasso_path_plot, width = 8, height = 5, dpi = 300)

#testing the panels
load('data/biomarker-clean.RData')

# ---- edit this to change proteins between 5 and 3
my_proteins <- c("DERM", "IgD", "14-3-3 protein zeta/delta", "M2-PK", "Coagulation Factor IX")

missing <- setdiff(my_proteins, names(biomarker_clean))
if (length(missing) > 0) stop("Not found in data: ", paste(missing, collapse = ", "))

biomarker_panel <- biomarker_clean %>%
  select(group, all_of(my_proteins)) %>%
  mutate(class = factor(group == 'ASD', levels = c(FALSE, TRUE), labels = c("TD", "ASD"))) %>%
  select(-group)

# Same seed and proportion as before
set.seed(101422)
biomarker_split <- initial_split(biomarker_panel, prop = 0.8)
train_data <- training(biomarker_split)
test_data  <- testing(biomarker_split)

# Logistic regression 
fit <- glm(reformulate(sprintf("`%s`", my_proteins), "class"),
           data = train_data,
           family = 'binomial')

class_metrics <- metric_set(sensitivity, specificity, accuracy, roc_auc)

test_data %>%
  mutate(
    pred = predict(fit, newdata = test_data, type = 'response'),
    pred_class = factor(ifelse(pred > 0.5, "ASD", "TD"), levels = c("TD", "ASD"))
  ) %>%
  class_metrics(
    truth = class,
    estimate = pred_class,
    pred,
    event_level = 'second'
  )

#number of test_data
nrow(test_data)

# For reference, the in-class baseline:
#   sensitivity 0.812 | specificity 0.733 | accuracy 0.774 | roc_auc 0.883
