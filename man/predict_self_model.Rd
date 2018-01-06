\name{predict_self_model}
\alias{predict_self_model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predicting CMRs sites by user-specific model
}
\description{
Predicting the m5C through user-specific squences.
}
\usage{
predict_self_model(models,sequence_dir,end = 5,up = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{models}{A dataframe self models and selected feature based extra_model.}
   \item{sequence_dir}{A path representing the filename of the sequence in FASTA format.}
   \item{up}{A integer number, the length of the upstream sequence required.}
   \item{end}{A integer number, the length of the downstream sequence required.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
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
ddd <- extra_model(res = aaa)
ddd
#
eee <- predict_self_model(models = ddd,sequence_dir = paste0(system.file(package = "PEAm5c"),"/data/cdna.fa"))
table(eee[,4])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

