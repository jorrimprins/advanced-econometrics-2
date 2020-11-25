use "/Users/jorrimprins/Google Drive/MSc Econometrics/Advanced Econometrics 2/Week 4/Assignment/AE2prodfunc.dta", clear
*ssc install xtabond2
*ssc install estout
*mata: mata set matafavor speed, perm

xtset id year
* "create dummies for all years and drop 1"
tab year, gen(YEAR)

* OLS

* normal regression of y on l and k + lagged variables
reg y l l.l k l.k l.y
* delete entity fixed effects by within estimation
xtreg y l l.l k l.k l.y, fe
* "work with lambda_t, year fixed effects"
xtreg y l l.l k l.k l.y YEAR*, fe
* "same model but with clustered errors"
xtreg y l l.l k l.k l.y YEAR*, fe cluster(id)


* IV
*we always use robust standard errors (Windmeijer 2005) but have shown the difference it makes in first 2 regressions
 
* iv regression where y is seen as predetermined and l and k as endogeneous, no robust errors to show difference
ivreg2 d.y (d.l ld.l d.k ld.k ld.y = l2.l l3.l l2.k l3.k l2.y) YEAR*, gmm2s noconstant
* iv regression where y is seen as predetermined and l and k as endogeneous
ivreg2 d.y (d.l ld.l d.k ld.k ld.y = l2.l l3.l l2.k l3.k l2.y) YEAR*, gmm2s noconstant robust cluster(id)
* same iv regression but use extra instruments, we would want even more though so we use recycling
ivreg2 d.y (d.l ld.l d.k ld.k ld.y = l2.l l3.l l4.l l2.k l3.k l4.k l2.y l3.y) YEAR*, gmm2s noconstant robust cluster(id)


* GMM
* note that we don't have to use cluster(id) as this would give exactly the same results,
* start by doing same regression as iv, but using recycled instruments
* we do not have to use specific instruments for the lags of l and k as those instruments are already there

* t-2 as recycled instrument for l k and y
xtabond2 y l l.l k l.k l.y YEAR*, gmm(l k y, lag(2 2)) iv(YEAR*) noleveleq twostep robust

* collapse lags for y, well identified
xtabond2 y l l.l k l.k l.y YEAR*, gmm(l, lag(2 2)) gmm(k, lag(2 2)) gmm(y, lag(2 2) c) iv(YEAR*) noleveleq twostep robust

* t-2 and t-3 as recycled instruments for l
xtabond2 y l l.l k l.k l.y YEAR*, gmm(l, lag(2 3)) gmm(k, lag(2 2)) gmm(y, lag(2 2) c) iv(YEAR*) noleveleq twostep robust


quietly reg y l l.l k l.k l.y
est store spec1
qui xtreg y l l.l k l.k l.y, fe 
est store spec2
qui xtreg y l l.l k l.k l.y YEAR*, fe
est store spec3
qui xtreg y l l.l k l.k l.y, fe cluster(id)
est store spec4

** Next, print key estimation results (b and se) for the fta regressor only

esttab spec1 spec2 spec3 spec4, b se keep(l L.l k L.k L.y) sfmt(%9.3f) star(* 0.10 ** 0.05) mtitles nonumbers compress title(Table 1: Just an example)

quietly xtabond2 y l l.l k l.k l.y YEAR*, gmm(l k y, lag(2 2)) iv(YEAR*) noleveleq twostep robust
testnl (_b[l] + _b[l.l] + _b[k] + _b[l.k])/(1 - _b[l.y]) = 1

* collapse lags for y, well identified
quietly xtabond2 y l l.l k l.k l.y YEAR*, gmm(l, lag(2 2)) gmm(k, lag(2 2)) gmm(y, lag(2 2) c) iv(YEAR*) noleveleq twostep robust
testnl (_b[l] + _b[l.l] + _b[k] + _b[l.k])/(1 - _b[l.y]) = 1

* t-2 and t-3 as recycled instruments for l
quietly xtabond2 y l l.l k l.k l.y YEAR*, gmm(l, lag(2 3)) gmm(k, lag(2 2)) gmm(y, lag(2 2) c) iv(YEAR*) noleveleq twostep robust
testnl (_b[l] + _b[l.l] + _b[k] + _b[l.k])/(1 - _b[l.y]) = 1










