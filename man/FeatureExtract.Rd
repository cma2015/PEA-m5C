\name{FeatureExtract}
\alias{FeatureExtract}
\title{
Feature encoding
}
\description{
  This function contains three feature encoding scheme, binary, k-mer and PseDNC. For binary encoding scheme, a vector of 404 (4*101) features is generated through assigning 'A', 'C', 'G', 'U' and 'N' with (1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,1) and (0,0,0,0), respectively. Here 'N' is a gap used to ensure the fixed features of each sample, if an m6A/non- m6A site occurs near the initiation or termination of the transcript. For K-mer encoding, the composition of short sequence with different lengths was considered to encoding samples. For PseDNC (pseudo dinucleotide composition) encoding, the local and global sequence-order information along the RNA sequence was used for scoring the each sample.
}
\usage{
  FeatureExtract(RNAseq, lambda = 6, w = 0.9)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RNAseq}{A list containing the FASTA format sequences.}
     \item{lambda}{The lambda parameter for the PseDNC-related features, default is 6.}
  \item{w}{The weighting parameter for PseDNC-related features, default is 0.9.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  A matrix with features.
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
bbb[1:10,]

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ feature encoding }% use one of  RShowDoc("KEYWORDS")
\keyword{ PseDNC }% __ONLY ONE__ keyword per line
