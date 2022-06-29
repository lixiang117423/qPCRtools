# generate NAMESPACE and build and check
cp DESCRIPTION CRAN/qPCRtools/DESCRIPTION
rm CRAN/qPCRtools/NAMESPACE
Rscript roxygen2.R
R CMD build CRAN/qPCRtools
R CMD check *.tar.gz --as-cran

# add function CalExp2ddCt
cp deve/R/CalRTable.R  CRAN/qPCRtools/R/

git add *
git commit -m "add function CalRTable"
git push origin


# add function CalCurve
cp deve/R/CalCurve.R  CRAN/qPCRtools/R/

git add *
git commit -m "add function CalCurve"
git push origin

# add function CalExpCurve
cp deve/R/CalExpCurve.R  CRAN/qPCRtools/R/

git add *
git commit -m "add function CalExpCurve"
git push origin

# add function CalExp2ddCt
cp deve/R/CalExp2ddCt.R  CRAN/qPCRtools/R/

git add *
git commit -m "add function CalExp2ddCt"
git push origin

# edit README
git add *
git commit -m "Edit README"
git push origin


# add function CalExpRqPCR
cp deve/R/CalExpRqPCR.R  CRAN/qPCRtools/R/

git add *
git commit -m "add function CalExp2ddCt"
git push origin
