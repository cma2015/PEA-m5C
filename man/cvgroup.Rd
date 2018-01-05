\name{cvgroup}
\alias{cvgroup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cross validation grouping
}
\description{
Grouping data sets by cross validation
}



\usage{
cvgroup(data,cvnum,seed=1234)
}


\arguments{
  \item{data}{A numeric vector, the number of samples (such as 1:1000)}
  \item{cvnum}{A numeric value, the number of Cross validation group (k-fold)}
  \item{seed}{Arguments passed to set.seed()}
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
