\name{PEA_ml}
\alias{PEA_ml}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Machine learning-based Randomforest algorithm and information gain rank for sequences
}
\description{
First of all, the function separate the positive and negative sample to the independent test set automatically according to setting and training set, concentrated training shall be carried out in accordance with the cross validation in training and assessment, in the results returned to characteristics of the evaluation, prediction score and evaluation as well as the sample group. In this case, the function can define the proportion of positive and negative samples (when negative sample size can be selected to maintain multiple models), and provide two methods of random forest and support vector machine.
}
\usage{
PEA_ml(pos_sample,neg_sample,independent_num=100,ig="ALL",
                   times = 1,modeltype = "RFC",cvnum = 5,repeatTimes = 1, ntree=200,over_sampling = F)
                   }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pos_sample}{A numeric matrix recording the features for positive sample.}
  \item{neg_sample}{A numeric matrix recording the features for nagative sample.}
  \item{independent_num}{A numeric value, how many independent samples what you want}
  \item{ig}{A numeric value,vector or "ALL", how many features what you want base the top of information gain rank}
  \item{modeltype}{A character string, machine learing method.}
  \item{cvnum}{An integer value, cross-times cross validation.}
  \item{repeatTimes}{An integer value,If the negative sample is larger than the limit of the positive sample, the number of the negative samples and the number of samples of the positive sample is repeated
}
  \item{over_sampling}{F or T, where TRUE represents balance the positive and negative samples according to the ratio based smote simulation
}
  \item{times}{A numeric value, where 1 represents balance the positive and negative sample.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of result.

The first level is used feature num group.

The second level is cross validation group.

The third level is the detail information including

    positives.test.score.id

    negatives.test.score.id

    positives.test.score

    negatives.test.score

    positives.test

    negatives.test

    auc_test

    auc_test_id
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
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
