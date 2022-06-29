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
usethis::use_package("rstatix")
usethis::use_package("ggthemes")
usethis::use_package("tibble")
usethis::use_package("tidyr")

# other
usethis::use_lifecycle_badge( "Experimental" )
usethis::use_version("patch")
usethis::use_version("minor")
#usethis::use_version("minor")


# create functions
golem::add_module("test")

# check
library(roxygen2)

roxygen2::roxygenize('CRAN/qPCRtools')

system('R CMD build CRAN/qPCRtools --binary')

system('R CMD check qPCRtools_0.1.1.tar.gz --as-cran')

devtools::install_local('qPCRtools_0.1.1.tar.gz')
