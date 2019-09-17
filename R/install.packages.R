if (!require("plyr")){
  install.packages("plyr", repos='http://cran.us.r-project.org', dependencies = TRUE)
  library(plyr)
}
if (!require("dplyr")){
  install.packages("dplyr", repos='http://cran.us.r-project.org', dependencies = TRUE)
  library(dplyr)
}
if (!require("ggplot2")){
  install.packages("ggplot2", repos='http://cran.us.r-project.org', dependencies = TRUE)
  library(ggplot2)
}

if (!require("EnvStats")){
  install.packages("EnvStats", repos='http://cran.us.r-project.org', dependencies = TRUE)
  library(EnvStats)
}
if (!require("psych")){
  install.packages("psych", repos='http://cran.us.r-project.org', dependencies = TRUE)
  library(psych)
}
if (!require("ggbeeswarm")){
  install.packages("ggbeeswarm", repos='http://cran.us.r-project.org', dependencies = TRUE)
  library(gplots)
}
if (!require("ggpubr")){
  install.packages("ggpubr", repos='http://cran.us.r-project.org', dependencies = TRUE)
  library(ggpubr)
}
