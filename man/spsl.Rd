\name{spsl}
\alias{spsl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  Function to exploit Independent SPSL model to estimate the dynamic gene coexpression among all the gene pairs systematically.
Fitting SPSL model
}
\description{
spsl is used to estimate the dynamic gene coexpression among all gene pairs with independent spike-and-slab prior in a Bayesian variable selection framework.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
spsl(x, z, n.iter, n.chains=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The matrix of gene expression data, where each row represents the genes of a patient and each column represents a gene.
%%     ~~Describe \code{x} here~~
}
\item{z}{
A vector of independent variable, where each element represents the cancer recurrence-free survival probability of a patient.
}
\item{n.iter}{number of iterations in MCMC}
\item{n.chains}{number of the MCMC chains}
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
Wenda Zhang
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
# set tau0<-(0,0,0,0,0,0,0,0,0,0)
# set tau1<-(0,0,0,0,0,0,0.1,0.2,0.3,0.5)
tau0<-rep(0,q)
wh<-rep(7:10,(q%/%10))+rep(0:((q%/%10)-1),rep(4,(q%/%10)))*10
tau1<-rep(0,q)
tau1[wh]<-c(0.1,0.2,0.3,0.5)
data <- get_data(n=100, nb=5, tau0=tau0, tau1=tau1)
x <- data[[1]]
z <- data[[2]]
# Implement SPSL model
result <- spsl(x, z, n.iter=100, n.chains=1)
summary(result[[1]])
cat("The computation time for MCMC is", result[[2]], "secs", "\n")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
