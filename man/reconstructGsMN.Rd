% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reconsructGsMN.R
\name{reconstructGsMN}
\alias{reconstructGsMN}
\title{Reconstuction of the specific-organism genome-scale metabolic network}
\usage{
reconstructGsMN(metabolic.data, RefData = RefDbcache, threshold = 10,
  is.gaint = TRUE)
}
\arguments{
\item{metabolic.data, }{df or a character vector. More details see function 
\code{getOrgMetabolicData} and \code{details}}

\item{RefData}{The reference metabolic data. It does not need reference data 
While organism metabolic data was collected from KEGG database, and RefData 
is set to NULL. Otherwise, RefDbCache, an internal dataset in this package,
was taken as the Reference metabolic data for Genome scale metabolic reconstruction.}

\item{threshold}{numeric, Nodes belonging to components with fewer than the 
value of threshold nodes will be ignored. This is a good option for networks
that contain many small and trivial components. Default is 10.}

\item{is.gaint}{logical, Ignore all nodes except those in the giant component:
selecting the only main largest component (connected set of nodes) of the
network. All smaller components will be ignored. This is a good option for
networks with a dominant component. Default is TRUE.}
}
\value{
igraph object
}
\description{
Reconstruction of genome-scale metabolic network (GsMN) whose 
 nodes represents compounds and whose edges represents reactions.
}
\details{
The input of this function can be of two forms. If organims is 
 collected in KEGG database, it can be obtained with 
 \code{getOrgMetabolicData} which is a data frame. Otherwise, 
 \code{metabolic.data} could be a  character vecotr which contains the KEGG 
 Orthology annotated information on this organism, e.g. we can download this 
 KO annotation profile in the \url{https://img.jgi.doe.gov} website for 
 species detected in a human microbime which not contained in KEGG organism 
 database. Several functions, such as \code{link{read.table}} and 
 \code{\link{read.delim}} could help us to read KO annotation profile.
}
\examples{
## not run (organism in KEGG)
## metabolic.data <- getOrgMetabolicData("buc")
## g <- reconstructGsMN(metabolic.data)

## species detected in a human microbiome
annodir <- system.file("extdata","koanno.tab",package = "RevEcoR")
metabolic.data <- read.delim2(file=annodir,stringsAsFactors=FALSE)
##load the reference metabolic data
data(RefDbcache)
g2 <- reconstructGsMN(metabolic.data, RefData = RefDbcache)

}
\seealso{
\code{\link{getOrgMetabolicData}}
}

