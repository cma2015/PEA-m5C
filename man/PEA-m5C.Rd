\name{predict_m5c}
\alias{predict_m5c}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predicting m5C sites by PEA-m5C
}
\description{
Predicting m5C information of the sequences by the proposed m5C models based Random Forest.
}
\usage{
predict_m5c(sample_feature)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{sample_feature}{A dataframe or list for target sequence.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A matrix with transcript , its sites, predicted score and m5C level.
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
aaa <- extra_motif_seq(input_seq_dir = paste0(system.file(package = "PEAm5c"),"/data/cdna.fa"),up = 5)
aaa <- lapply(aaa, c2s)
bbb <- FeatureExtract(aaa)
ccc <- predict_m5c(bbb)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ PEA-m5C }% use one of  RShowDoc("KEYWORDS")
