decomp.contour.cod <- function(c1, c2, ages=c(0, 15, 40, 65, 80, 111), FUN = ex.per, ...) {

  # decomposition by CoD country1 and country 2

  # input arguments: 
  # c1, c2 - data for calculations, data frames. First column ages, nex columns - rates by CoD 
  # FUN - function to calculate statistics of interest
  #
  # (c) Dmitri A. Jdanov
  # jdanov@demogr.mpg.de
  # last revised 07.05.2024
  
  
  require(gtools)
  Age = c1[, 1]
  agroups=cut(Age,ages,labels=FALSE, right = FALSE)
  na = length(Age)

  
  # vectors of rates
  n = dim(c1)[2]
  n1 = ((n - 1) / 2) + 1
  a = c1[, 2 : n1]
  b = c2[, 2 : n1]
  A = c1[, (n1 + 1) : n]
  B = c2[, (n1 + 1) : n]
  nC = n1 - 1
  
  # vectors of age-specific components of differences  
  #
  codperm = permutations(n = nC, r = nC, v = 1:nC)
  
  # A->B  
  dAB = decomp.contour.cod1(A, B, a, b, Age, na, nC, codperm, dir = 1, FUN = ex.per, ...)
  dBA = decomp.contour.cod1(B, A, b, a, Age, na, nC, codperm, dir = -1, FUN = ex.per, ...)
  D = (dBA - dAB) / 2
  write.table(D, file = "res.csv", quote = F, sep = ",", row.names = F, col.names = F)  
  return(D)
  
}

decomp.contour.cod1 <- function(A, B, a, b, Age, na, nC, codperm, dir = 1, FUN = ex.per, ...) {
  dfAB = matrix(0,na,1)
  dfABc  = matrix(0, na, 3 * (nC + 1))
  np = dim(codperm)[1]
  r = matrix(0, np, 4 * nC)
  for (x in 1:na) {
    fBA_1 = FUN(A, Age, ...)
    Aa = decomp.cod1(A, a, x, nC, Age, codperm, FUN = ex.per, ...)
    A[x, ] = a[x, ]
    ab = decomp.cod1(A, b, x, nC, Age, codperm, FUN = ex.per, ...)
    A[x, ] = b[x, ]
    bB = decomp.cod1(A, B, x, nC, Age, codperm, FUN = ex.per, ...)
    if (dir < 0) {
      dfABc[x, ] = c(bB, ab, Aa) 
    }
    else {
      dfABc[x, ] = c(Aa, ab, bB) 
    }
    A[x, ] = B[x, ]
    fBA = FUN(A, Age, ...)
    dfAB[x] = fBA - fBA_1 #component of conventional decomposition (for control)
  }
  return(cbind(dfAB, -1 * dfABc))
}

decomp.cod1 <- function(A, a, x, nC, Age, codperm, FUN = ex.per, ...) {
  fBA_1 = FUN(A, Age, ...)
  np = dim(codperm)[1]
  r = matrix(0, np, nC)
  for (p in 1 : np) {
    cp = codperm[p, ]
    A1 = A 
    fBA1_1 = fBA_1 #
    for (i in 1 : nC) {
      A1[x, cp[i]] = a[x, cp[i]]
      fAa = FUN(A1, Age, ...)
      r[p, cp[i]] = fBA1_1 - fAa
      fBA1_1 = fAa
    }
  }  
  A1[x, ] = a[x, ]
  fAa = FUN(A1, Age, ...)
  return(c(fBA_1 - fAa, colMeans(r)))
}





decomp.contour.cod1l <- function(A, B, a, b, Age, na, nC, codperm, FUN = ex.per, ...) {
  dfAB = matrix(0,na,1)
  dfABc  = matrix(0, na, 4 * nC)
  np = dim(codperm)[1]
  r = matrix(0, np, 4 * nC)
  for (x in 1:na) {
    fBA_1 = FUN(A, Age, ...)
    for (p in 1 : np) {
      cp = codperm[p, ]
      A1 = A 
      fBA1_1 = fBA_1 #
      for (i in 1 : nC) {
        A1[x, cp[i]] = a[x, cp[i]]
        fBa = FUN(A1, Age, ...)
        r[p, 4 * (cp[i] - 1) + 1] = fBA1_1 - fBa
        
        A1[x, cp[i]] = b[x, cp[i]]
        fBb = FUN(A1, Age, ...)
        r[p, 4 * (cp[i] - 1) + 2] = fBa - fBb
        
        A1[x, cp[i]] = B[x, cp[i]]
        fBB = FUN(A1, Age, ...)
        r[p, 4 * (cp[i] - 1) + 3] = fBb - fBB
        
        r[p, 4 * (cp[i] - 1) + 4] = fBA1_1 - fBB
        fBA1_1 = fBB
      }
      dfABc[x, ] = colMeans(r)
    }  
    A[x, ] = B[x, ]
    fBA = FUN(A, Age, ...)
    dfAB[x] = fBA - fBA_1 #component of conventional decomposition (for control)
  }
  return(cbind(dfAB, -1 * dfABc))
}