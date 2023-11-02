# create description
golem::fill_desc(
  pkg_name = "qPCRtools", # The Name of the package containing the App
  pkg_title = "Tools for qPCR", # The Title of the package containing the App
  pkg_description = "PKG_DESC.", # The Description of the package containing the App
  author_first_name = "Xiang", # Your First Name
  author_last_name = "LI", # Your Last Name
  author_email = "lixiang117423@gmail.com", # Your Email
  repo_url = "https://github.com/lixiang117423/qPCRtools" # The URL of the GitHub Repo (optional)
)

# use license
usethis::use_mit_license( "Xiang LI" )

# insert packages
usethis::use_package("dplyr")
usethis::use_package("magrittr")
usethis::use_package("ggplot2")
usethis::use_package("ggpmisc")
usethis::use_package("stringr")
usethis::use_package("data.table")
usethis::use_package("multcomp")
usethis::use_package("broom")
usethis::use_package("readxl")
usethis::use_package("reshape2")
usethis::use_package("xlsx")
usethis::use_package("sjmisc")

# other
usethis::use_lifecycle_badge( "Experimental" )
# 版本号，DESCRIPTION文件发生变化
usethis::use_version("major")
usethis::use_version("minor")
usethis::use_version("patch")


# create functions
# golem::add_module("test")

# copy files
file.copy("./inst/", "./CRAN/qPCRtools/", recursive = TRUE)
file.copy("./man/", "./CRAN/qPCRtools/", recursive = TRUE)
file.copy("./DESCRIPTION", "./CRAN/qPCRtools/", recursive = TRUE)
file.copy("./R/", "./CRAN/qPCRtools/", recursive = TRUE)
file.copy("./LICENSE", "./CRAN/qPCRtools/", recursive = TRUE)

# check and build
library(roxygen2)

roxygen2::roxygenize('../qPCRtools')

system('R CMD build ../qPCRtools --binary')

system('R CMD check ggmotif_0.1.3.tar.gz --as-cran')

devtools::install_local('ggmotif_0.1.3.tar.gz')



