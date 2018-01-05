\name{extra_motif_seq}
\alias{extra_motif_seq}
\title{
Scanning the motifs in the transcripts.
}
\description{
For a given motif, all transcript sequences will be scanned and then according to the size of "upstream" "downstream",
the output fits the sequence of length.
}
\usage{
extra_motif_seq(input_seq_dir,text='c',up=5,end=5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{input_seq_dir}{A path character, representing the file name of the sequence in FASTA format.}
   \item{text}{A string character, which specifies the motif to be searched.}
   \item{up}{A integer number, the length of the upstream sequence required.}
   \item{end}{A integer number, the length of the downstream sequence required.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of sequences around motif.
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
aaa
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ extra motif seq }% use one of  RShowDoc("KEYWORDS")
