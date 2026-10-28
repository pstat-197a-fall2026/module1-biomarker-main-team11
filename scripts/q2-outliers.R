## QUESTION 2: Outliers
## Temporarily remove the outlier trimming from preprocessing and explore the
## outlying values. Are there specific SUBJECTS that seem to be outliers? If
## so, are outliers more frequent in one group or the other?
##
## Run from the repository root (same convention as preprocessing.R).
## Reads data/biomarker-raw.csv only. The untrimmed data are kept in memory and
## are NOT saved, so data/biomarker-clean.RData is left untouched.

library(tidyverse)

dir.create('results', showWarnings = FALSE)

# ---- 1. preprocessing WITHOUT trimming ----

var_names <- read_csv('biomarker-raw.csv',
                      col_names = F,
                      n_max = 2,
                      col_select = -(1:2)) %>%
  t() %>%
  as_tibble() %>%
  rename(name = V1,
         abbreviation = V2) %>%
  na.omit()

# same as preprocessing.R (log10, center, scale) but with the trim() step removed
biomarker_untrimmed <- read_csv('biomarker-raw.csv',
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


# Are specific subjects outliers? Yes.
# 
# Under a normal distribution you’d expect about 3.6 outlying values (|z| > 3) per subject, and the typical subject has a median of 8.5. A small set of subjects have far more.
# 
# 13 subjects exceed the upper fence (about 29 outlying values).
# Together they account for 46% of all outlying values.
# The most extreme are subject 154 (TD, 157 outliers), 108 (TD, 127), 9 (ASD, 126), 121 (TD, 122), 52 (ASD, 121), 77 (TD, 114), and 147 (TD, 77).
# The outliers are subject-driven, not protein-driven. No single protein has more than 7 outlying values.
# Subject IDs are row positions after dropping the blank first row.
# 
# Are outliers more frequent in one group? Not significantly.
# 
# TD has more by raw counts: 1,372 outlying values vs. 1,007 for ASD, and a mean of 17.6 vs. 13.2 per subject.
# 9 of the 13 flagged subjects are TD, and 4 are ASD.
# The medians are identical (8.5), so the TD excess comes from a few extreme subjects.
# None of the tests find a difference: Wilcoxon p = 0.82, permutation test on mean counts p = 0.29, and Fisher exact test on flagged subjects p = 0.25.
# 
# Conclusion for your write-up: a handful of subjects, mostly TD, have unusually many extreme values. There is no statistical evidence that one group is more outlier-prone.
# 
# The cutoff of 3 matches the trimming rule in preprocessing.R; a different cutoff would change the counts.
