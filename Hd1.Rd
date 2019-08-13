\name{Hd1}
\alias{Hd1}

\title{
association between variable and run in cohort
}
\description{

\usage{
Hd1(m)
}

\arguments{
  \item{m}{
}
}
\details{
}
\value{

}
\references{

}
\author{

}
\note{

}



\seealso{
}
\examples{

batch<- rep(gl(10,50,lables=c(1:10)),1)
m<-data.frame(batch,matrix(rnorm(60000),nrow=500))

Hd1 <- function(m)
{
  a <- sweep(m[-1],2, colMeans(m[-1]), "-");
  b<- svd(a);
  x <- b$u
  y <- b$v
  c<-svd(m[-1])
  z<-c$u
  corr <- sapply(1:ncol(x),function(i){
    fit <- lm(x[,i]~as.factor(m[,1]))
    return( summary(fit)$adj.r.squared  )
  })
  corr.nscale <- sapply(1:ncol(z),function(i){
    fit <- lm(z[,i]~as.factor(m[,1]))
    return( summary(fit)$adj.r.squared  )
  })
  pdf("corPlot_f2_Raw.pdf")
  par (mfrow=c(2,1))
  p1<-plot(seq(along=corr), corr, xlab="PC",main="scaled")
  p2<- plot(seq(along=corr.nscale), corr.nscale, xlab="PC",main="Non-scaled")
  names(corr)<- names(m[-1])
  names(corr.nscale)<- names(m[-1])
  write.csv(corr,file="correlation_data_f2.csv")
  write.csv(corr.nscale,file="correlation_data_f2_nscaled.csv")
  print(p1)
  print(p2)
  dev.off()
}
 }

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
