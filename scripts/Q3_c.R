## QUESTION 3(c): fuzzy intersection instead of a hard intersection

library(glmnet)
library(randomForest)
library(tidyverse)
library(tidymodels)

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
