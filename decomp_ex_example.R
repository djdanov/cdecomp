# an example of contour decomposition difference 
# in life expectancy at birth by age and cause of death 
# (c) Dmitri Jdanov 
# jdanov@demogr.mpg.de
# last revised 7.05.2024

cdecomp.cod.ex<-function(fname1 = "USA", fname2 = "ENW", years = c(2019, 2020), sex = 2) {

  # prepare the data
  c1 = mx.select(fname1, years, sex)
  c2 = mx.select(fname2, years, sex)

  # decomposition

  r = decomp.contour.cod(c1, c2, ages=c(seq(0, 100, by = 5), 111), FUN = ex.per, sex = sex)
  return(r)
}

ex.per <- function(mx, age, sex = '') {
  # calculates life table from mx (period data, HMD method, w/o smoothing) 
  # usage: ex.per(mx, sex = 'm')
  # input arguments: mx - vector of mortality rates (one year age group), sex - "m" for males and "f" for females
  # output: life expectancy
  # 
  # Last revised: 07.05.2024
  mx = rowSums(mx)
  mx [mx < 0] = 0 #NaN
  last = length(mx)
  ax = 0 * mx + 0.5
  nx = c(diff(age), 1)
  if (sex == '') {
    cat('warning: sex==m')
    sex = 'm'
  }
  # a(0) & a(1)
  if (mx[1] >= 0.107) {
    if (sex == 1) {
      ax[1] = 0.33
      ax[2] = 1.352 / 4
    }
    else {
      ax[1] = 0.35
      ax[2] = 1.361 / 4
    }
  }
  else {
    if (sex == 1) {
      ax[1] = 0.045 + 2.684 * mx[1]
      ax[2] = (1.653 - 3.013 * mx[1]) / 4
    }
    else {
      ax[1] = 0.053 + 2.8 * mx[1]
      ax[2] = (1.524 - 1.627 * mx[1]) / 4
    }
  }
  # last age
  ax[last] = 1 / mx[last]
  # q(x)
  qx = nx * mx / (1 + nx * (1 - ax) * mx)
  qx[is.na(qx)] = 1
  n = which.max(qx>1)
  if (n > 1) {
    qx = c(qx[1:(n - 1)], 1)
    ax = c(ax[1:(n - 1)], 0.5)
    nx = nx[1:n]
  }
  last=length(qx)
  qx[last] = 1 
  ax[is.na(ax)] = 0.5
  qx[is.na(qx)] = 1
  ax[last] = 1 / mx[last]
  px = 1 - qx
  lx = c(1, cumprod(px))
  lx = head(lx, -1)
  dx = lx * qx
  dx[length(dx)] = lx[length(lx)]
#  Lx = lx - (1 - ax) * dx
  Lx = nx * lx - nx * (1 - ax) * dx
  Lx[last] = lx[last] * ax[last]
  # return e(0)
  return(sum(Lx))
}

mx.select<-function(fname, country, years, sex) {
    # selects data from the HCD file
    # input arguments: fname, years (vector), causes (matrix)
    # returns data frame with columns Age, mx1, mx2, 
    # where mx1 and mx2 are age specific mortality rates for years specified in vector years
    ages = c(0, 1, seq(5, 85, 5))
    
    a = length(ages)
    dta = read.table(fname, header = T, sep = ",",  strip.white=TRUE, na.string=".")
    dta1 = cbind(ages, dta[dta$Country == country & dta$Year == years[1] & dta$Sex == sex, 6 : 9] / 100000)
    dta2 = dta[dta$Country == country & dta$Year == years[2] & dta$Sex == sex, 6 : 9] / 100000
    r = data.frame(cbind(dta1, dta2))
    return(r)
}

