loadOrInstall <- function (packageName, type = "CRAN") {
  if (type == "CRAN") {
    isPackageInstalled <- packageName %in% rownames(installed.packages())
    if (!isPackageInstalled) {
      install.packages(packageName)
    } 
    library(packageName, character.only = TRUE)
  }
  
  if (type == "git") {
    isPackageInstalled <- strsplit(packageName,split = "/",fixed = TRUE)[[1]][2] %in% rownames(installed.packages())
    if (!isPackageInstalled) {
      install_github(packageName)
    } 
    library(strsplit(packageName,split = "/",fixed = TRUE)[[1]][2], 
            character.only = TRUE)
  }
}

onlyInstall <- function (packageName, type = "CRAN") {
  isPackageInstalled <- packageName %in% rownames(installed.packages())
  if (!isPackageInstalled) {
    if (type == "CRAN") {
      install.packages(packageName)
    } 
  }
}

cranPackages <- c(
  'shiny',
  'rhandsontable',
  'shinyjs',
  'openxlsx',
  'shinythemes',
  "tibble",
  "purrr",
  "dplyr",
  "cowplot",
  "ggplot2",
  "ggiraph",
  "ggraph",
  "igraph",
  "devtools",
  "colorspace"
)

cranPackages2 <- c(
  "colorspace"
)

gitPackages <- c(
  "sandmanns/clevRvis"
)

for (package in cranPackages) {
  loadOrInstall(package)
}
#for (package in cranPackages2) {
#  onlyInstall(package)
#}
for (package in gitPackages) {
  loadOrInstall(package,"git")
}


source("util/Phylogeny.R")#code for creating a phylogeny
source("util/CheckTable.R")#code for checking phylogeny table
source("util/Mutations.R")#code for mutation creation
source("util/CNVchecks.R")#code for checking cnv table


#for colored phylogeny table
color_renderer = "
    function(instance, td, row, col, prop, value, cellProperties) {
        Handsontable.renderers.TextRenderer.apply(this, arguments);
        if (instance.params) {
            clr = instance.params.ColorCode
            clr = clr instanceof Array ? clr : [clr]
            td.style.background =  clr[row]
            td.style.fontWeight = 'bold'
            clr2 = instance.params.ColorCode2
            clr2 = clr2 instanceof Array ? clr2 : [clr2]
            td.style.color =  clr2[row]
        }
    }"













