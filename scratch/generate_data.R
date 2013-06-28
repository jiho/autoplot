
# PCA data
# contains several quantitative variables (some to be used as supplementary variables), one qualitative variable (to be used as supplementary)
# contains enough rows so that some can be used as supplementary objects

set.seed(123)
d = data.frame(
  a=runif(5), b=rnorm(5), c=rnorm(5)*3,
  d=c("bob","bob","joe","joe","joe")
)
row.names(d) = paste("foo", 1:5)

write.table(d, file="pca-data.csv", row.names=F, sep=",")



# CA data
# contains three quantitative variables which sum to 100 on each row
# contains enough rows so that some can be used as supplementary objects

d = data.frame(a=c(50, 20, 30, 30, 20), b=c(40, 40, 20, 10, 40), c=c(10, 40, 50, 60, 40))
row.names(d) = paste("foo", 1:5)

library("plyr")
aaply(d, 1, sum, .expand=F)

write.table(d, file="ca-data.csv", row.names=F, sep=",")
