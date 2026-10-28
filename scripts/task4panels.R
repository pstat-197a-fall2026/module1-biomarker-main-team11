library(tidyverse)
library(tidymodels)
library(yardstick)
load('data/biomarker-clean.RData')


#baseline
#The in class analysis uses 5 proteins: "DERM"       "RELT"       "MRC2"       "IgD"        "Cadherin-5"
#sensitivity 0.812 | specificity 0.733 | accuracy 0.774 | roc_auc 0.883

# alternative panel:
# 5 proteins: c("DERM", "IgD", "14-3-3 protein zeta/delta", "M2-PK", "Coagulation Factor IX")
# sensitivity: 0.688, specificity: 0.9333, accuracy: 0.806, roc_auc: 0.933
#analysis: this panel of proteins gives more weight to specificity. trades lower sensitivity for higher specificity and overall accuracy. 

#attempted smaller panel: 
# 3 proteins:  c("DERM", "IgD", "14-3-3 protein zeta/delta")
# sensitivity: 0.562, specificity: 0.8, accuracy: 0.677, roc_auc: 0.8
# 4 proteins: c("DERM", "IgD", "14-3-3 protein zeta/delta", "M2-PK")
# sensitivity: 0.5, specificity: 0.933, accuracy: 0.710, roc_auc: 0.825
# the results for the smaller panel are noticeably worse.

#summary: 
# used lasso to find an alternative panel with the same number of proteins. 
#2 shared proteins with the base panel: DERM AND IgD
#results tested on the same logistic regression.
#the alternative panel with 5 proteins was able to achieve a better result 
#trades lower sensitivity for higher specificity and overall roc_auc Unfortunately, the smaller panels from this method were not able to achieve comparable classification accuracy, dropping 6% ROC_AUC compared to baseline when going down to 4 proteins. Loses even more ROC_AUC when going down to 3, 2, or 1 proteins. 
#caveats: the difference in accuracy from this split is 1 correct classification (small test set)





# ---- edit this to change proteins between 5 and 3
my_proteins <- c("DERM")

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
