# generate NAMESPACE and build and check
rm CRAN/qPCRtools/NAMESPACE
Rscript roxygen2.R
R CMD build CRAN/qPCRtools
R CMD check *.tar.gz --as-cran


# add function CalRTable
cp deve/R/CalRTable.R  CRAN/qPCRtools/R/

git add *
git commit -m "add function CalRTable"
git push origin


# add function CalCurve
cp deve/R/CalCurve.R  CRAN/qPCRtools/R/

git add *
git commit -m "add function CalCurve"
git push origin


# edit README
git add *
git commit -m "Edit README"
git push origin
