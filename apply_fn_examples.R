# Author: Jingyang (Judy) Zhang
# Date: Dec. 9th, 2025
# Purpose: an example of the "apply" family functions.


# 1. apply(): 
## Sum of each row.
A <- matrix(1:9, 3)
A
apply(A, 1, sum)
apply(A, 2, sum)

# 2. lapply():
B <- list(a = 1:3, b = 4:6)
B
lapply(B, mean)

# 3. sapply()
sapply(B, mean)


# 4. vapply()
vapply(B, mean, numeric(1))

# 5. mapply()
C1 <- c(1, 2, 3)
C2 <- c(10, 20, 30)
mapply(sum, C1, C2)


# 6. tapply()
data(iris)
tapply(iris$Sepal.Length, iris$Species, mean)

# 7. by()
by(iris[,1:4], iris$Species, colMeans)
library(tidyverse)
iris %>% group_by(Species) %>% summarise(across(where(is.numeric), mean))

aggregate(.~Species, data = iris, FUN = mean)

library(data.table)
dt <- as.data.table(iris)
dt[, lapply(.SD, mean), by = Species]