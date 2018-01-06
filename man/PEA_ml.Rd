\name{PEA_ml}
\alias{PEA_ml}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transcriptome-wide m5C predictor training under machine learning framework.
}
\description{
This function used for transcriptome-wide m5C predictor construction. First, the fixed number (parameter "independent_num") of independent samples (positive and negative samples) are randomly sampled from training samples. The the k-fold cross-validation would be performed based on the training samples but excluding the independent samples. Finally, the m5C predictor, performance evaluation on independent test datasets and cross-validation results will be returned.
}
\usage{
PEA_ml(pos_sample,neg_sample,independent_num=100,ig="ALL",
                   ratio = 1,modeltype = "RFC",cvnum = 5,repeatTimes = 1, ntree=200,over_sampling = F)
                   }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pos_sample}{A numeric matrix recording the features for positive sample.}
  \item{neg_sample}{A numeric matrix recording the features for nagative sample.}
  \item{independent_num}{A numeric value, the number of independent sample }
  \item{feature_num}{A numeric value, the number of selected features based the top of information gain rank, the "ALL" means all features}
  \item{modeltype}{A character string, which specifies machine learing method.}
  \item{cvnum}{An integer value, the number of fold for cross validation.}
  \item{repeatTimes}{An integer value,If the negative sample is larger than the limit of the positive sample, the number of the negative samples and the number of samples of the positive sample is repeated
}
  \item{over_sampling}{Logical value, where TRUE represents balance the positive and negative samples according to the ratio based smote simulation
}
  \item{ratio}{A numeric value, where 1 represents balance the positive and negative sample.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of result.

The first level is used feature num group.

The second level is cross validation group.

The third level is the detail information including  positives.test.score.id, negatives.test.score.id, positives.test.score,negatives.test.score, positives.test, negatives.test, auc_test, auc_test_id
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
   Jie Song, Jingjing Zhai, Chuang Ma
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
load(paste0(system.file(package = "PEAm5c"),"/data/samples.Rds"))
aaa <- PEA_ml(pos_sample = pos_sample,neg_sample = neg_sample)
aaa
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
