#!/broad/tools/apps/R-2.6.0/bin/Rscript

args <- commandArgs(TRUE)
verbose <- TRUE
file.name <- args[1]
out.file.name <- args[2]
header <- TRUE
sep <- ""
extra.arg.delim <- "="
cex <- 1
height.scale <- 1
width.scale <- 1
value.delim <- ','
K <- 0.08
se.correct <- FALSE

is.qt <- FALSE
only.map <- FALSE

gene.col <- 'GENE'
p.col <- 'P'
maf.col <- 'MAF'
car.freq.col <- NULL
beta.col <- 'BETA'
se.col <- 'SE'

if (length(args) > 2)
{
  for (extra.arg in args[3:length(args)])
  {
    cur.arg <- as.character(unlist(strsplit(extra.arg, extra.arg.delim)))
    cur.arg[1] <- as.character(cur.arg[1])
    if (length(cur.arg) == 2)
    {
      if (cur.arg[1] == "sep")
      {
        sep <- cur.arg[2]
      } else if (cur.arg[1] == "header")
      {
        header <- as.logical(cur.arg[2])
      } else if (cur.arg[1] == "gene.col")
      {
        gene.col <- cur.arg[2]
      } else if (cur.arg[1] == "p.col")
      {
        p.col <- cur.arg[2]
      } else if (cur.arg[1] == "maf.col")
      {
        maf.col <- cur.arg[2]
      } else if (cur.arg[1] == "car.freq.col")
      {
        car.freq.col <- cur.arg[2]
      } else if (cur.arg[1] == "beta.col")
      {
        beta.col <- cur.arg[2]
      } else if (cur.arg[1] == "se.col")
      {
        se.col <- cur.arg[2]
      } else if (cur.arg[1] == "se.correct")
      {
        se.correct <- as.logical(cur.arg[2])
      } else if (cur.arg[1] == "prevalence")
      {
        K <- as.numeric(cur.arg[2])
      } else if (cur.arg[1] == "only.map")
      {
        only.map <- as.logical(cur.arg[2])
      } else if (cur.arg[1] == "is.qt")
      {
        #this is not implemented yet -- we need to define analog of func.Vg but for QTs
        #is.qt <- as.logical(cur.arg[2])
      }
    }
  }
}

if (sep == '\\t')
{
  sep = "\t"
}

to.read <- file.name
if (file.name == "/dev/stdin")
{
  to.read <- file(description="stdin")
}

to.write <- out.file.name
if (out.file.name == "/dev/stdout")
{
  to.write <- file(description="stdout")
}


table <- read.table(to.read, header=header, sep=sep)

get.cols <- function(table, initial.col)
{
  complement <- FALSE
  if (substr(initial.col, 1, 1) == "-")
  {
    complement <- TRUE
    initial.col <- substr(initial.col, 2, nchar(initial.col))
  }

  initial.cols <- as.character(unlist(strsplit(initial.col, value.delim)))
  for (i in 1:length(initial.cols))
  {
    match <- colnames(table) == initial.cols[i] | colnames(table) == make.names(initial.cols[i])

    if (sum(match) > 0)
    {
      initial.cols[i] <- (1:length(colnames(table)))[match][1]
    } else if (is.na(as.integer(initial.cols[i])))
    {
      stop(paste("Bad column",initial.cols[i]))
    }
  }
  initial.cols <- unique(as.integer(initial.cols))

  if (complement)
  {
    initial.cols <- (1:ncol(table))[-initial.cols]
  }
  return(initial.cols)
}

gene.col <- get.cols(table, gene.col)
p.col <- get.cols(table, p.col)
if (!is.null(car.freq.col)) {
  car.freq.col <- get.cols(table, car.freq.col)
} else
{
  maf.col <- get.cols(table, maf.col)
}

beta.col <- get.cols(table, beta.col)
se.col <- get.cols(table, se.col)

func.Vg <- function (PA,RR1,RR2,K) {
  Paa = (1-PA)^2
  PAa = 2*PA*(1-PA)
  PAA = PA^2
  muaa=0
  faa= K/(Paa + PAa*RR1 + PAA*RR2)
  fAa= RR1*faa
  fAA= RR2*faa
  T = qnorm(1-faa)
  muAa = T-qnorm(1-fAa)
  muAA = T-qnorm(1-fAA)
  mean.all= PAa*muAa+ PAA*muAA
  Vg= Paa*(muaa-mean.all)^2 + PAa*(muAa-mean.all)^2+ PAA*(muAA-mean.all)^2
  actual.Vg =  Vg/(1+Vg)
  VR = 1-actual.Vg
  actual.T = Paa*sqrt(VR)*qnorm(1-faa) + PAa*sqrt(VR)*qnorm(1-fAa) + PAA*sqrt(VR)*qnorm(1-fAA)
  actual.muaa = actual.T - sqrt(VR) * qnorm(1-faa)
  actual.muAa = actual.T - sqrt(VR) * qnorm(1-fAa)
  actual.muAA = actual.T - sqrt(VR) * qnorm(1-fAA)

  res <- list(Vg=actual.Vg,muaa=actual.muaa, muAa = actual.muAa, muAA=actual.muAA)
  res

}

##example
#func.Vg(PA=0.3,RR1=1.3,RR2=1.3^2,K=0.01)

convert.to.rr <- function(maf, effect) {
  f1 <- 2 * maf * (1-maf)
  f2 <- maf ** 2

  OR1 <- exp(effect)
  OR2 <- OR1 ** 2

  RR1 <- 1.1
  RR2 <- 1.1**2
  tol <- 1e-6
  max.it <- 100
  prev.RR1 <- NULL
  prev.RR2 <- NULL
  for (j in seq(max.it)) {
    f.aa <- K / ((1-f1-f2) + f1 * RR1 + f2 * RR2)
    RR1 <- OR1 / (1 + f.aa * (OR1 - 1))
    RR2 <- OR2 / (1 + f.aa * (OR2 - 1))
    if (!is.null(prev.RR1) && abs(prev.RR1 - RR1) < tol && abs(prev.RR2 - RR2) < tol) {
      break
    }
  }
  return(c(RR1,RR2))
}

result <- data.frame(Gene=c(), P=c(), MAF=c(), Effect=c(), SE=c(), LVE=c(), LVE_MAP=c())

for (i in seq(nrow(table))) {

  effect <- table[i,beta.col]
  se <- table[i,se.col]
  if (is.null(car.freq.col)) {
    maf <- table[i,maf.col]
  } else
  {
    car.freq <- table[i,car.freq.col]
    a <- -1
    b <- 2
    c <- -car.freq

    maf <- (-b + sqrt(b**2 - 4 * a * c)) / 2 * a
  }

  if (maf > .5) {
    maf <- 1 - maf
    effect <- -effect
  }

 integrand <- function(x) {
    
    #return(dnorm(x,table[i,beta.col],table[i,se.col]))

    vg <- sapply(x, function(y) { RR <- convert.to.rr(maf, y); func.Vg(maf,RR1=RR[1],RR2=RR[2],K=K)$Vg*dnorm(y,table[i,beta.col],table[i,se.col])})
    vg[is.na(vg)] <- 0

    return(vg)
  }

  RR <- convert.to.rr(maf, effect)
  RR.se<- convert.to.rr(maf, se)
  lve.map <- func.Vg(maf,RR1=RR[1],RR2=RR[2],K=K)$Vg
  lve.se.map <- func.Vg(maf,RR1=RR.se[1],RR2=RR.se[2],K=K)$Vg  
  if (is.qt) {
    lve <- func(maf, effect, table[i,se.col])
    lve.se <- func(maf, se, table[i,se.col])
  } else {
    if (!only.map) {
      lve <- integrate(integrand, -Inf, Inf)$value
      lve.se <- lve.se.map
    } else
    {
      lve <- NA
    }
  }
  gene.name <- table[i,gene.col]

  se.correct.term <- 0
  se.correct.map.term <- 0
  if (se.correct) {
    se.correct.term <- lve.se
    se.correct.map.term <- lve.se.map
  }    

  cur.results <- c(gene.name, table[i,p.col], maf, effect, table[i,se.col], lve - se.correct.term, lve.map - se.correct.map.term)
  if (se.correct) {
    cur.results <- c(cur.results, lve, lve.map)
  }
  result <- rbind(result, cur.results)

}
cols <- c("Gene", "P", "MAF", "Effect", "SE", "LVE", "LVE_MAP")
if (se.correct)
{
  cols <- c(cols, "LVE_uncorrected", "LVE_MAP_uncorrected")
}
colnames(result) <- cols

result$Gene <- table[,gene.col]

write.table(result, file=out.file.name, quote=FALSE, sep="\t", row.names=FALSE)
