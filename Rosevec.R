# Rosevec.R -- try to vectorize Rosenbrock function
# 20231201 -- Thanks to Ivan for sorting this out for me, JN
onerose<-function(i,x, scale){
   if ( ! (i %in% 1:2)) stop("bad index i")
   if (i == 1) {
       res <- scale*(x[2] - x[1]^2)
   } else {
       res <- (1 - x[1])
   }
   res
}

onerose_v <- Vectorize(onerose, 'i')
allres <- function(x, scale) onerose_v(1:2, x, scale)
x <- c(-1.2, 1)
library(adagio)
print(x)
cat("adagio::fnRosenbrock(x):",fnRosenbrock(x),"\n")
aa <- allres(x, scale=10)
print(aa)
print(sum(aa^2))

