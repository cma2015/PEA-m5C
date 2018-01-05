\name{extra_model}
\alias{extra_model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
extracting model from PEA_ml
}
\description{
Proposing the desired model from the established result information
}
\usage{
extra_model(res,ignum)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{a list , result obtained from PEA_ml}
  \item{ignum}{a integer number,which indicates the number of features}}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
models and selected feature
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
ddd <- extra_model(res = aaa,ignum=150)
ddd
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
