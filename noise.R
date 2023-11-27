source("definitions.R")
kw <- 1e-14
HEPTOT <- 0.05
# For maximum compatibility, choose the oldest version of R's RNG for which
# RNGversion doesn't warn
RNGversion('3.6.0')
set.seed(1365)

params <- expand.grid(
	pK_true = replicate(16, sort(runif(3, -3, 9)), FALSE),
	replication = 1:16,
	noise_ratio = 10^seq(-4, -.5, len = 8),
	n_points = c(10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000),
	noise = list(
		normal = rnorm,
		Cauchy = rcauchy,
		uniform = function(n) runif(n, -1, 1)
	)
)

vnorm <- function(x) sqrt(sum(x^2))

params$m <- lapply(seq_len(nrow(params)), function(i) {
	p <- params[i,]
	SID <- seq(-5e-2, 4e-2, length.out = p$n_points)
	pH <- pHm(SID, HEPTOT, p$pK_true[[1]][1], p$pK_true[[1]][2], p$pK_true[[1]][3])
	noise <- p$noise[[1]](length(pH))
	pH <- pH + noise * vnorm(pH) / vnorm(noise) * p$noise_ratio
	nlsr::nlxb(
		pH ~ pHm(SID, HEPTOT, pK1, pK2, pK3),
		data.frame(pH = pH, HEPTOT=HEPTOT, SID = SID),
		start = c(pK1 = -1, pK2 =  3, pK3 =  7.55),
		lower = c(pK1 = -3, pK2 = -3, pK3 = -3   ),
		upper = c(pK1 =  9, pK2 =  9, pK3 =  9   ),
		control = list(japprox = 'jacentral')
	)
})

saveRDS(params, 'noise.rda', version = 2, compress = 'xz')
