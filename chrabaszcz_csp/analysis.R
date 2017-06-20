library(dplyr)
library(tidyr)
library(ggplot2)
library(foreign)
library(rstan)
library(data.table)
library(foreach)
library(reshape2)
options(stringsAsFactors = FALSE)
theme_set(theme_classic())

dat <- read.csv("data.csv")

### Munge
dat <- data.table(dat)
cntl <- dat[,grep("24|25|18",names(dat)),with=FALSE] # names are reversed
expmt <- dat[,grep("62|63|26",names(dat)),with=FALSE]
tmp <- c("p1","p2","p3","est.low","est.hi")
setnames(cntl, tmp)
setnames(expmt, tmp)
dat.joe <- data.table(dat[,1:2,with=FALSE],
                      rbind(cntl,expmt),
                      condition=rep(c("Experimental", "Control"),each=nrow(dat)))

# Convert judgments to probabilities
dat.joe[,c("p1","p2","p3"):=data.table(cbind(p1,p2,p3)/100)]

# Define quantiles
dat.joe[, q1:=est.low+(est.hi-est.low)/6]
dat.joe[, q2:=(est.hi+est.low)/2]
dat.joe[, q3:=est.hi-(est.hi-est.low)/6]

# Objective Function
optimFn <- function(par, probs, quants) {
  sum((probs - pnorm(quants,par[1],par[2]))^2)  
}

# Fitting Function
fitFn <- function(subj, optimFn) {
  probs <- unlist(subj[,list(p1,p2,p3)])
  quants <- unlist(subj[,list(q1,q2,q3)])
  
  if(!is.na(sum(probs)) & !is.na(sum(quants))) {
    mu <- seq(from=mean(quants)*.5, to=mean(quants)*1.5, length.out=10)
    sigma <- seq(from=var(quants)*.5, to=var(quants)*1.5, length.out=10)
    start.values <- data.frame(expand.grid(mu,sigma))
    names(start.values) <- c("mu","sigma")
    
    fit.current <- optim(par=c(start.values$mu[1],start.values$sigma[1]),
                         fn=optimFn,
                         probs=probs, 
                         quants=quants)  
    
    foreach(j=2:nrow(start.values)) %do% {
      fit <-  optim(par=c(start.values$mu[j],start.values$sigma[j]),
                    fn=optimFn,
                    probs=probs, 
                    quants=quants
      )  
      if (fit$value < fit.current$value) { fit.current <- fit }
    }
    
    out <- c(fit.current$par,fit.current$value)
  }  else {
    out <- c(rep(NA,times=3))
  }    
  return(out)
}

# Fit estimates
out <- foreach(i=1:nrow(dat.joe)) %dopar% {
    fitFn(dat.joe[i,],optimFn) 
}
out <- do.call(rbind,out)
dat.joe[,c("mu","sigma","sse"):=data.table(out)]

# Error
dat.joe[,SS.tot := 2*apply(dat.joe[,.SD,.SDcols=c("q1","q2","q3")],
                         1,var,na.rm=TRUE)]
dat.joe[, SS.fit := SS.tot - sse]
dat.joe[, R.sq := 1-(sse/SS.tot)]

### Plots

getX <- function(par.1, par.2, steps) {
  seq(par.1-4*par.2, par.1+4*par.2, length=steps)
}

# Get by subject distributions 
steps <- 100
dat.joe <- dat.joe[R.sq>.5,]
dat.joe[, idx:=1:.N]
dat.dists <- dat.joe[, list(x = x <- getX(mu, sigma, steps),
                            F.x = pnorm(x, mu, sigma),
                            exp,
                            condition,
                            mu,
                            sigma),
                     by=idx]

# Get consensus distributions
dat.con <- dat.joe[,list(mu=median(mu),sigma=median(sigma)), 
                   by=c("exp","condition")]
dat.con[,idx:=1:.N]
dat.con <- dat.con[,list(x = x <- getX(mu,sigma,steps),
                            F.x = pnorm(x, mu, sigma),
                            exp,
                            condition,
                            mu,
                            sigma),
                     by=idx]

endpoints <- dat.con[,list(low=quantile(x, .02), 
                            high=quantile(x, .98)),by=exp]
setkey(dat.dists, exp)
setkey(endpoints, exp)
dat.CDF <- foreach(e = endpoints$exp, .combine='rbind') %do% {
  dat.dists[e, ][x>endpoints[e,low] & x <endpoints[e,high],]  
}
dat.CDF[,model:='raw']

p.CDF <- ggplot(aes(x=x, y=F.x, group=idx), data=dat.CDF) + 
         geom_line(alpha=.5, color="grey") + 
         labs(title="Consensus Distributions", x="x", y="F(x)") +
         facet_wrap(~exp, scales="free", nrow=2)

p.con <- p.CDF + geom_line(aes(x=x, y=F.x, 
                      group=condition, color=condition), 
                  data=dat.con) + theme_classic(base_size = 36)

ggsave(file = "figure/cdfs.pdf", width = 40, height = 10, limitsize = FALSE)

# PDFs

# Get by subject distributions 
dat.d  <- dat.joe[, list(x = x <- getX(mu, sigma, steps),
                            f.x = dnorm(x, mu, sigma),
                            exp,
                            condition,
                            mu,
                            sigma),
                     by=idx]

# Get consensus distributions
dat.con.d <- dat.joe[,list(mu=median(mu),sigma=median(sigma)), 
                   by=c("exp","condition")]
dat.con.d[,idx:=1:.N]
dat.con.d <- dat.con.d[,list(x = x <- getX(mu,sigma,steps),
                            f.x = dnorm(x, mu, sigma),
                            exp,
                            condition,
                            mu,
                            sigma),
                     by=idx]

endpoints <- dat.con.d[,list(low=quantile(x, .02), 
                            high=quantile(x, .98)),by=exp]
setkey(dat.d, exp)
setkey(endpoints, exp)
dat.PDF <- foreach(e = endpoints$exp, .combine='rbind') %do% {
  dat.d[e, ][x>endpoints[e,low] & x <endpoints[e,high],]  
}

p.PDF <- ggplot(aes(x=x, y=f.x, group=idx), data=dat.PDF) + 
         geom_line(alpha=.5, color="grey") + 
         labs(title="Consensus Distributions", x="x", y="F(x)") +
         facet_wrap(~exp, scales="free", ncol=2)

p.con.d <- p.PDF + geom_line(aes(x=x, y=f.x, 
                      group=condition, color=condition), 
                  data=dat.con.d) 

### median effect size aggregation

obs <- read.csv("study.csv")
obs %<>% group_by(exp) %>%
  mutate(es = (exp.mean - control.mean) / sqrt(control.sd ^ 2 + exp.sd ^ 2))

temp <- dat.joe %>%
  select(V1, exp, condition, mu, sigma) %>%
  gather(parameter, value, -V1, -exp, -condition) %>%
  dcast(V1 + exp ~ condition + parameter) %>%
  group_by(V1, exp) %>%
  summarise(prior.es = (Experimental_mu - Control_mu) / sqrt(Control_sigma ^ 2 + Experimental_sigma ^ 2)) %>%
  group_by(exp) %>%
  summarise(prior.es = median(prior.es, na.rm = TRUE))

obs <- inner_join(obs, temp)

# effect on t
t.effect <- cbind(obs$t,
  abs((obs$es * obs$n + obs$prior.es) * sqrt(obs$n - 2) / (obs$n + 1))
)

save(dat.joe, file = "jeff.RData")

### stan code
# each question is a separate
# expect intermittent initialization failures, some models may need to be rerun (20 June 2017)

rm(list = ls())
load("jeff.RData")
chains <- 5
iters <- 4000
mc.cores <- 3

for (pick.study in unique(dat.joe$exp)) {
  jtmp <- filter(dat.joe, exp == pick.study, est.low < est.hi, p1 < p2, p2 < p3, p1 < p3) %>%
    filter(V1 %in% V1[duplicated(V1)]) %>% # select subjects who were monotonic for both conditions
    dplyr::select(p1:p3, condition, q1:q3, V1) %>%
    melt(id.vars = c("condition", "V1")) %>%
    mutate(variable = factor(grepl("p", variable), labels = c("prompt", "prob"))) %>%
    group_by(variable) %>%
    mutate(idx = 1:n()) %>%
    spread(variable, value) 
  mu <- mean(jtmp$prompt)
  sigma <- sd(jtmp$prompt)
  jtmp$prompt <- (jtmp$prompt - mu) / sigma

  jtmp1 <- list(prompt = jtmp$prompt,
             prob = jtmp$prob,
             condition = unclass(factor(jtmp$condition)),
             idx = unclass(factor(jtmp$V1)),
             N = nrow(jtmp),
             J = length(unique(jtmp$condition)),
             K = length(unique(jtmp$V1)),
             a = c(1, 1, 10))

  sfit <- stan(data = jtmp1, file = "beta2_hier.stan", chains = chains*3, iter = iters, cores = mc.cores)
  save(sfit, mu, sigma, file = paste0("q", pick.study, ".RData"))
  if (pick.study == 4) {
    rm(sfit)
    sfit <- stan(data = jtmp1, file = "beta2_hier_extra.stan", chains = chains, iter = iters, cores = mc.cores)
    save(sfit, mu, sigma, file = paste0("q", pick.study, "_pooling.RData"))
  }
  rm(sfit)
}
