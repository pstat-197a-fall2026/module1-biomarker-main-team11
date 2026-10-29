
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
ggsave("lasso_path.png", lasso_path_plot, width = 8, height = 5, dpi = 300)
#for combined report, move the path to results/report_files

