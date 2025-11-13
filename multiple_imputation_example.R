library(Hmisc)

set.seed(123)

n <- 200

# Generate synthetic dataset
df <- data.frame(
  age = sample(20:70, n, replace = TRUE),
  income = sample(30000:100000, n, replace = TRUE),
  gender = factor(sample(c("M", "F"), n, replace = TRUE)),
  smoker = factor(sample(c("Yes", "No"), n, replace = TRUE)),
  education = factor(sample(c("HighSchool", "Bachelors", "Masters", "PhD"), n, replace = TRUE))
)

# Introduce some missing values (~10% missing at random)
set.seed(456)
for (col in names(df)) {
  missing_idx <- sample(1:n, size = round(0.1 * n))
  df[missing_idx, col] <- NA
}

# Check the dataset
summary(df)


# Perform multiple imputation
imp <- aregImpute(
  ~ age + income + gender + smoker + education,
  data = df,
  n.impute = 5,  # generate 5 imputed datasets
  nk = 0         # use default spline knots for numeric variables
)


# Extract the first completed dataset as a list
imputed_list <- impute.transcan(
  x = imp,
  imputation = 1,  # first imputed dataset
  data = df,
  list.out = TRUE,  # important: list of vectors for each variable
  pr = FALSE
)

# Create a copy of the original dataframe
completed_df <- df

# Replace missing values in original dataframe with imputed values
for (var in names(imputed_list)) {
  completed_df[[var]][is.na(completed_df[[var]])] <- imputed_list[[var]]
}

# Check that there are no missing values
summary(completed_df)
