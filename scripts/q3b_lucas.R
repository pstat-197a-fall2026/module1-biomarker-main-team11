library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)

load("../data/biomarker-clean.RData")


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