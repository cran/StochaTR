\name{solveSLP}
\alias{solveSLP}
\title{
Solving stochastic linear programs with a single risk constraint
}
\description{
Solving stochastic linear programs with a single risk constraint. Coherent distortion risk measures are used.    
}
\usage{
solveSLP(fname = "SLPparams.dat", fdir = getwd())
}
\arguments{
  \item{fname}{
The name of the data file.
}
  \item{fdir}{
A path to the directory where the input and output files should be located. The default value is the \R working directory.
}
}
\details{
The points must be in the general position.
}

\value{
\item{\file{slp-WMTR.dat}}{The solution.} 
}
\references{
Bazovkin, P. and Mosler, K. Stochastic linear programming with distortion risk constraints. 201x. Working paper.
}
\author{
Pavel Bazovkin and Karl Mosler
}

\examples{
ls()
}

\keyword{ programming }
\keyword{ IO }
\keyword{ multivariate }
\keyword{ nonparametric }

