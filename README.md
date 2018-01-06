# **PEAm5C**: An integrated R toolkit for plant m5C analysis. </br>
![](https://halobi.com/wp-content/uploads/2016/08/r_logo.png "R logo")
![](https://encrypted-tbn2.gstatic.com/images?q=tbn:ANd9GcSvCvZWbl922EJkjahQ5gmTpcvsYr3ujQBpMdyX-YG99vGWfTAmfw "linux logo")
![](https://tctechcrunch2011.files.wordpress.com/2014/06/apple_topic.png?w=220)
<br>
We developed PEA-m5C, an accurate transcriptome-wide m5C modification predictor under machine learning framework with random forest algorithm. PEA-m5C was trained with features from the flanking sequences of m5C modifications. In addition, we also deposited all the candidate m5C modification sites in the Ara-m5C database (http://bioinfo.nwafu.edu.cn/software/Ara-m5C.html) for follow-up functional mechanism researches. Finally, in order to maximize the usage of PEA-m5C, we implement it into a cross-platform, user-friendly and interactive interface and an R package named “PEA-m5C” based R statistical language and JAVA programming language, which may advance functional researches of m5C.
<br>
## Version and download <br>
* [Version 0.11--R](https://github.com/cma2015/PEA-m5C/blob/master/PEAm5C_0.1.1.tar.gz) -First version released on January, 6th, 2018<br>
*  [Version 0.1--JAVA](https://github.com/cma2015/PEA-m5C/blob/master/PEA-m5C-java.zip) -First version released on January, 6th, 2018<br>
## Depends
#### R environment <br>
* [R](https://www.r-project.org/) (>= 3.3.1) <br>
* [randomForst](https://cran.r-project.org/web/packages/randomForest/index.html) (>= 0.6) <br>
* [seqinr](https://cran.rstudio.com/web/packages/seqinr/index.html) (>= 3.4-5) <br>
* [stringr](https://cran.r-project.org/web/packages/stringr/index.html) (>= 1.2.0) <br>
* [pROC](https://cran.rstudio.com/web/packages/pROC/index.html) (>= 1.10.0) <br>
* [ggplot2](https://bioconductor.org/packages/release/bioc/html/ggplot2.html) (>= 2.2.1) <br>
* [FSelector](https://cran.r-project.org/web/packages/FSelector/) (>= 0.21) <br>
#### Global software environment <br>
* [JAVA1.8](http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html) Environmentally dependent <br>
<br>

#### Dependency installation <br>
```bash
## Install rJAVA
sudo apt-get update
sudo apt-get install r-cran-rjava r-cran-rweka
```
```R
## Install R Dependency
dependency.packages <- c("randomForest", "seqinr", "stringr", "FSelector", "bigmemory", "ggplot2", "PRROC", "pROC")
install.packages(dependency.packages)
```
<br>

## Installation <br>
```R
install.packages("Download path/PEAm5C_0.11.tar.gz",repos = NULL, type = "source")
```
## Contents <br>
#### Predicting m5C sites <br>
* Read FASTA file and motif scanning <br>
* Feature encoding  of sequences <br>
* m5C prediction using Random Forest models <br>
<br>
#### user-defined model<br>
* Provide positive and negative sample information<br>
* Automatic verification of the training process<br>
* Prediction using user-defined models <br>
 <br>
## Quick start <br>
The basic data set can be finded in [data](https://github.com/cma2015/PEA-m5C/tree/master/data). <br>
More details can be seen from [user manual](https://github.com/cma2015/PEA-m5C/blob/master/PEAm5c.pdf).
<br>

#### 1.Predicting m5C sites <br>
* 1.1 Read FASTA file and motif scanning <br>
```R
aaa <- extra_motif_seq(input_seq_dir = paste0(system.file(package = "PEAm5c"),"/data/cdna.fa"),up = 5)
aaa <- lapply(aaa, c2s)
```
* 1.2 Feature encoding  of sequences <br>
```R
bbb <- FeatureExtract(aaa)
```
* 1.3 m5C prediction using Random Forest models  <br>
```R
ccc <- predict_m5c(bbb)
```
#### 2.User-defined model <br>
* 2.1 Provide positive and negative sample information <br>
```R
load(paste0(system.file(package = "PEAm5c"),"/data/samples.Rds"))
### The positive and negative sequence can be read and identified by extra_motif_seq and  feature encoding by FeatureExtract 
```
* 2.2 Automatic verification of the training process <br>
```R
aaa <- PEA_ml(pos_sample = pos_sample,neg_sample = neg_sample)
ddd <- extra_model(res = aaa)
ddd

```
* 2.3 Prediction using user-defined models <br>
```R
eee <- predict_self_model(models = ddd,sequence_dir = paste0(system.file(package = "PEAm5c"),"/data/cdna.fa"))
table(eee[,4])
```


## Ask questions
Please use [PEAm5C/issues](https://github.com/cma2015/PEAm5C/issues) for how to use PEAm5C and reporting bugs.
