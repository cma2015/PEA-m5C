\name{extra_motif_seq}
\alias{extra_motif_seq}
\title{
Scanning specific motifs in the transcripts.
}
\description{
For a given motif, all transcript sequences will be scanned according to the user-specific flanking sequence length (upstream/downstream), then the fixed length of sequences centered on motif will be returned.
}
\usage{
extra_motif_seq(input_seq_dir,text='c',up=5,end=5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{input_seq_dir}{A string character, representing the directory of the sequence in FASTA format.}
   \item{text}{A string character, which specifies the motif to be searched.}
   \item{up}{A integer number, the length of upstream sequence to be extracted.}
   \item{end}{A integer number, the length of downstream sequence to be extracted.}
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
