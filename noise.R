library(parallel)
library(pbapply)
loadNamespace('depcache')
stopifnot(`Please run this script with a cluster object registered` = !is.null(cl <- getDefaultCluster()))

source("definitions.R")
kw <- 1e-14
HEPTOT <- 0.05
# For maximum compatibility, choose the oldest version of R's RNG for which
# RNGversion doesn't warn
RNGversion('3.6.0')
set.seed(1365)

params <- expand.grid(
	replication = seq_len(10),
	pK_true = replicate(10, sort(runif(3, -3, 9)), FALSE),
	noise_ratio = 10^seq(-4, -.5, len = 8),
	noise = list(
		normal = rnorm,
		Cauchy = rcauchy,
		uniform = function(n) runif(n, -1, 1)
	),
	n_points = rev(c(10, 20, 50, 100, 200, 500, 1000, 2000))
)

getrows <- function(d) lapply(seq_len(nrow(d)), function(i) d[i,])
vnorm <- function(x) sqrt(sum(x^2))

# Set the noise in stone before we start using parallel operations.
# Far from the most efficient way to do this, but at least it's reproducible.
params$vnoise <- pblapply(getrows(params), function(p) p$noise[[1]](p$n_points))
# Send over all variables except the giant 'params' and the cluster itself.
clusterExport(cl, setdiff(ls(), c('cl','params')))
#pboptions(use_lb = TRUE)

params$m <- pblapply(getrows(params), function(p) {
	SID <- seq(-5e-2, 4e-2, length.out = p$n_points)
	pH <- pHm(SID, HEPTOT, p$pK_true[[1]][1], p$pK_true[[1]][2], p$pK_true[[1]][3])
	pH <- pH + p$vnoise[[1]] * vnorm(pH) / vnorm(p$vnoise[[1]]) * p$noise_ratio
	m <- depcache::cache(nlsr::nlxb(
		pH ~ pHm(SID, HEPTOT, pK1, pK2, pK3),
		data.frame(pH = pH, HEPTOT=HEPTOT, SID = SID),
		start = c(pK1 = -1, pK2 =  3, pK3 =  7.55),
		lower = c(pK1 = -3, pK2 = -3, pK3 = -3   ),
		upper = c(pK1 =  9, pK2 =  9, pK3 =  9   ),
		control = list(japprox = 'jacentral')
	))
	m$jacobian <- NULL # too large to store
	m$weights <- NULL # we'll just remember they are all ones
	m
}, cl = cl)
params$pK_est <- pblapply(params$m, coef)
params$noise_kind <- names(params$noise)
saveRDS(
	params[c("pK_true", "noise_ratio", "n_points", "pK_est", "noise_kind")],
	'noise.rds', version = 2, compress = 'xz'
)
