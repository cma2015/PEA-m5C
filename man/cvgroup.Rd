\name{cvgroup}
\alias{cvgroup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cross validation grouping
}
\description{
Grouping datasets for cross validation.
}

\usage{
cvgroup(data,cvnum,seed=1234)
}


\arguments{
  \item{data}{A string vector representing the names of training samples or a string vector recording the names of training samples.}
  \item{cvnum}{which specifies the number of fold in cross-validation.}
  \item{seed}{A numeric value used for generating random seed.}
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
Jingjing Zhai, Chuang Ma, Jie Song.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
cvgroup(1:100,5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
