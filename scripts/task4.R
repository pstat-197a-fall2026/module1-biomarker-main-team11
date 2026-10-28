
# This script selects the protein panel using Lasso.
#task4panels.R has the panels being tested and a brief analysis to help with write up


library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)
load('data/biomarker-clean.RData')


library(glmnet)

library(glmnet)

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