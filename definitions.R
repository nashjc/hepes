# The charge balance function for a weak acid with three given pKs.  The
# root of this function gives the equilibrium value of [H] at the given
# conditions.
HEPESFUNC <- function(H,SID,HEPTOT,pK1,pK2,pK3) {
	XX <- (H^3/(10^-pK3*10^-pK2*10^-pK1)+H^2/(10^-pK3*10^-pK2)+H/(10^-pK3))
	IV <- HEPTOT/(1+XX)
	I <- IV*H^3/(10^-pK3*10^-pK2*10^-pK1)
	II <- IV*H^2/(10^-pK3*10^-pK2)
	III <- IV*H/10^-pK3
	H - kw/H + SID + I*2 + II - IV
}

# Return a vector of equilibrium pHs for the given SIDs and other
# constant parameters.
pHm <- Vectorize(function(SID, HEPTOT, pK1, pK2, pK3)
  -log10(uniroot(HEPESFUNC,c(1e-20,1),tol=1e-20,maxiter=1e4,
                 HEPTOT=HEPTOT,SID = SID, pK1=pK1,pK2=pK2,pK3=pK3)$root))
