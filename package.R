library(devtools)
#create("/home/malab8/zjj/FTgenePrediction/RAP")
setwd("/home/malab14/research/DeepRNAMethy/script/PEA")
dir()
load_all()
document(roclets = "namespace")
build()
check()
install.packages("/home/malab14/research/DeepRNAMethy/script/PEA_1.0.tar.gz", repos = NULL, type="source")


pack <- "PEAm5C"
path <- find.package(pack)
system(paste(shQuote(file.path(R.home("bin"), "R")),
             "CMD", "Rd2pdf", shQuote(path)))
