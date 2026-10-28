library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)

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



# INTERPRETATION & DISCUSSION (For Report / Submission)

# 1. Methodological Justification & Data Leakage Prevention:
#
# In the initial in-class baseline analysis, feature selection (t-tests and
# Random Forest importance) was conducted on the full dataset prior to
# partitioning. This could introduce data leakage because information from
# the held-out test cases would influence predictor selection and potentially
# lead to overly optimistic performance estimates.
#
# To reduce data leakage, we implemented the following procedure:
# - Dataset Allocation: The dataset contained 154 samples, with 122 assigned
#   to the training set and 32 assigned to the test set.
# - Leakage-Free Feature Selection: T-tests with multiple-testing correction
#   and Random Forest feature importance ranking were performed using only
#   the training data.
# - Feature Selection: The intersection of the two sets of top 10 proteins
#   identified three predictors: DERM, IgD, and TrkC.
# - Model Fitting and Evaluation: Logistic regression was fitted using the
#   selected proteins in the training set and evaluated on the held-out
#   test set.

# 2. Performance & Interpretation:
#
# On the held-out test set, the model achieved the following results:
# - ROC AUC: 0.812
# - Accuracy: 0.750
# - Sensitivity: 0.812
# - Specificity: 0.688
#
# The test set contained 16 ASD samples and 16 TD samples. The model correctly
# classified 24 of the 32 samples, corresponding to an accuracy of 0.750.
# Sensitivity of 0.812 indicates that 13 of the 16 ASD samples were correctly
# identified, while specificity of 0.688 indicates that 11 of the 16 TD
# samples were correctly identified.
#
# Overall, the results suggest that the model shows some ability to
# distinguish ASD from TD samples in this held-out test set. However, the
# test set contains only 32 samples, so the performance estimates may be
# uncertain. These results provide preliminary evidence of predictive
# performance but do not establish that the model will generalize reliably
# to other populations or serve as a clinically useful early-detection tool.
