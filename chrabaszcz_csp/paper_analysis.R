library(dplyr)
library(tidyr)
library(ggplot2)
library(foreign)
library(rstan)
library(data.table)
library(foreach)
library(reshape2)
library(gdata)
options(stringsAsFactors = FALSE)
theme_set(theme_classic())

load("jeff.RData")
load("abstract.RData")
dat.joe$exp <- as.character(dat.joe$exp)
qs <- c("1", "4", "5", "6", "7", "8", "9", "11")

### Abstract elicitation stuff

abstract$prob <- abstract$Q39_1 / 100
abstract$exp <- factor(abstract$exp, labels = 1:8)

abs.mean <- data.frame(exp = levels(abstract$exp), prob = round(tapply(abstract$prob, abstract$exp, mean), 2))

# histograms of abstract responses
abs.plot <- ggplot(abstract, aes(prob)) +
  geom_density() +
  geom_text(data = abs.mean, x = .5, y = .1, aes(label = prob)) +
  # geom_histogram(binwidth = .05) +
  facet_wrap(~exp, nrow = 2) +
  theme_classic(base_size = 11) +
  xlab("Estimated statement that effect is true") +
  ylab("Density") +
  scale_x_continuous(breaks = c(.25, .75), expand = c(0, 0)) +           
  scale_y_continuous(expand = c(0, 0))

### Comparing csprior to original studies

# read in original study data, calculate effect sizes
obs <- read.csv("studies.csv")
obs <- obs %>%
  group_by(exp) %>%
  mutate(emp.diff = exp.mean - control.mean,
         emp.sd = sqrt((control.sd ^ 2 + exp.sd ^ 2)))
obs$exp <- as.character(obs$exp)
obs$t <- ifelse(obs$exp.mean > obs$control.mean, obs$t, -obs$t)

# function to calculate mean diff and pooled sd from bayesian models
esFun <- function(x) {
  load(paste0("q", x, ".RData", sep = ""))
  out <- summary(sfit, pars = c("mu", "sigma"))$summary[,1]
  out <- out * sigma
  out[1:2] <- out[1:2] + mu
  data.frame(exp = x,
             bayes.diff = out[2] - out[1],
             bayes.sd = sqrt((sum(out[3] ^ 2) + sum(out[4] ^ 2))))
}

# mean diff and pooled sd from ML process
tmp <- dat.joe %>%
  group_by(exp, condition) %>%
  summarise(mu = median(mu), sigma = median(sigma)) %>%
  gather(param, value, -exp, -condition) %>%
  dcast(exp ~ condition + param) %>%
  group_by(exp) %>%
  summarise(ml_diff = Experimental_mu - Control_mu,
            ml_sd = sqrt((Control_sigma ^ 2 + Experimental_sigma ^ 2)))

# estiamte Bayesian model effect sizes
bayes.pool <- do.call("rbind", lapply(qs, esFun))
names(bayes.pool)[2] <- "bayes.diff"
obs <- inner_join(obs, bayes.pool) %>% inner_join(tmp)

# ES table for paper
effest <- cbind(obs[c("emp.diff", "ml_diff", "bayes.diff")] / obs[c("emp.sd", "ml_sd", "bayes.sd")])
effest <- cbind(obs$emp.diff, obs$emp.sd, effest[,1], obs$ml_diff, obs$ml_sd, effest[,2], obs$bayes.diff, obs$bayes.sd, effest[,3])

# effect on t
simEffEst <- function(x, obs) {
  load(paste0("q", x, ".RData", sep = ""))
  obs <- obs[obs$exp == x,]
  tmp <- melt(as.array(sfit)) %>%
    filter(parameters %in% c("mu[1]", "mu[2]", "sigma[1]", "sigma[2]")) %>%
    spread(parameters, value)
  tmp <- tmp[,3:6] * sigma
  tmp[,1] <- (obs$n * obs$control.mean + tmp[,1] + mu) / obs$n + 1
  tmp[,2] <- (obs$n * obs$exp.mean + tmp[,2] + mu) / obs$n + 1
  tmp[,3] <- sqrt((obs$n * obs$control.sd ^ 2 + tmp[,3] ^ 2) / obs$n + 1)
  tmp[,4] <- sqrt((obs$n * obs$exp.sd ^ 2 + tmp[,4] ^ 2) / obs$n + 1)


  median((tmp[,2] - tmp[,1]) / sqrt((tmp[,3] ^ 2 + tmp[,4] ^ 2) / (obs$n - 2)))
}

obs <- dat.joe %>%
  group_by(exp, condition) %>%
  summarise(mu = median(mu), sigma = median(sigma)) %>%
  gather(param, value, -exp, -condition) %>%
  dcast(exp ~ condition + param) %>%
  group_by(exp) %>%
  summarise(ml_c_mean = Control_mu,
            ml_e_mean = Experimental_mu,
            ml_c_sd   = Control_sigma,
            ml_e_sd   = Experimental_sigma) %>%
  left_join(dplyr::select(obs, control.mean, control.sd, exp.mean, exp.sd, n), .) %>%
  transmute(ml.t = (((n * exp.mean + ml_e_mean) - (n * control.mean + ml_c_mean)) / (n + 1)) / sqrt((sqrt((n * control.sd ^ 2 + ml_c_sd ^ 2) / n + 1) ^ 2 + sqrt((n * exp.sd ^ 2 + ml_e_sd ^ 2) / n + 1) ^ 2) / (n - 2))) %>%
  left_join(obs, .)

obs$bayes.t <- sapply(qs, simEffEst, obs)

ttab <- obs[c("exp", "t", "ml.t", "bayes.t")]
ttab$exp <- 1:8
names(ttab) <- c("Q", "Empirical", "MLE", "HBE")

### type m/s calculations

retrodesign <- function(A, s, alpha=.05, df=Inf, n.sims=10000){
  z <- qt(1-alpha/2, df)
  p.hi <- 1 - pt(z-A/s, df)
  p.lo <- pt(-z-A/s, df)
  power <- p.hi + p.lo
  typeS <- p.lo/power
  estimate <- A + s*rt(n.sims,df)
  significant <- abs(estimate) > s*z
  exaggeration <- mean(abs(estimate)[significant])/A
  return(data.frame("Pwr" = power, "Sn"=typeS, "Mag"=exaggeration))
}

tmp <- abs(obs[c("emp.diff", "ml_diff", "bayes.diff")])
tmp2 <- obs[c("emp.sd", "ml_sd", "bayes.sd")]

smerror <- lapply(1:ncol(tmp), function(x) do.call("rbind",
  Map(retrodesign, unlist(tmp[,x]), 
                   unlist(tmp2[,x]) / sqrt(obs$n),
                   df = obs$n - 2)))

### sample size planning

nFun <- function(dif, sig) {
  round((2.8 * sig / dif) ^ 2)
}

tmp <- do.call("cbind", Map(nFun, tmp, tmp2))

samp.size <- cbind(1:8, obs$n, tmp)
colnames(samp.size) <- c("Q", "N", "Empirical", "MLE", "HBM")

### what % convinced?

# function to aggregate individual ES estimates from HBMs
indESFun <- function(x) {
  load(paste0("q", x, ".RData", sep = ""))
  tmp <- summary(sfit, pars = c("mu_idx", "sigma_idx"))$summary[,1]
  tmp <- data.frame(val = tmp, pars = names(tmp))
  tmp <- separate(tmp, pars, into = c("par", "idx", "condition"), sep = "[\\[|,]") %>%
    dcast(idx ~ par + condition, value.var = "val")
  data.frame(es = (tmp$"mu_idx_2]" - tmp$"mu_idx_1]") / sqrt((tmp$"sigma_idx_1]" ^ 2 + tmp$"sigma_idx_2]" ^ 2)), exp = x, method = "HBM")
}

bayes.idx <- do.call("rbind", lapply(qs, indESFun))

# individual ML ESs
tmp <- dat.joe %>%
  dplyr::select(V1, exp, condition, mu, sigma) %>%
  gather(par, val, -V1, -exp, -condition) %>%
  dcast(exp + V1 ~ par + condition, value.var = "val") %>%
  mutate(es = (mu_Experimental - mu_Control) / sqrt((sigma_Control ^ 2 + sigma_Experimental ^ 2)),
         method = "ML") %>%
  dplyr::select(es, exp, method)

tmp <- rbind(bayes.idx, tmp)

### Big aggregation plots

# raw estimates

tmp <- dat.joe
tmp$exp <- factor(tmp$exp, labels = as.character(1:8))
# best.guess <- ggplot(tmp, aes(q2, fill = condition)) +
best.guess <- ggplot(tmp, aes(q2, linetype = condition)) +
  scale_x_log10(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_density(alpha = .5, adjust = 3) +
  facet_wrap(~exp, nrow = 2, scales = "free_y") +
  theme_classic(base_size = 11) +
  xlab("Response") +
  ylab("Density") +
  theme(legend.position = "top", legend.title = element_blank()) +
  scale_fill_brewer(type = "qual", palette = 6)

# big stuff

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
dat.dists$Method <- "Median"

# get Bayesian stuff
bayesDists <- function(x) {
  load(paste0("q", x, ".RData", sep = ""))
  tmp <- melt(as.array(sfit))
  tmp <- drop.levels(tmp[grepl("mu_", tmp$parameters) | grepl("sigma_idx", tmp$parameters),]) %>%
    group_by(parameters) %>%
    summarise(est = mean(value)) %>%
    separate(parameters, into = c("par", "idx", "condition"), sep = "_|,", )
  tmp <- dcast(tmp, idx + condition ~ par, value.var = "est")
  tmp$mu <- tmp$mu * sigma + mu
  tmp$sigma <- tmp$sigma * sigma
  tmp$exp <- x
  tmp
}

tmp <- as.data.table(do.call("rbind", lapply(qs, bayesDists)))
tmp[, idx:=1:.N]
tmp <- tmp[, list(x = x <- getX(mu, sigma, steps),
                            F.x = pnorm(x, mu, sigma),
                            exp,
                            condition,
                            mu,
                            sigma),
                     by=idx]
tmp$Method <- "Bayes"
dat.dists <- rbind(dat.dists, tmp)
# end bayes

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
dat.con$Method <- "Median"

# more bayes
bayesCon <- function(x) {
  load(paste0("q", x, ".RData", sep = ""))
  tmp <- summary(sfit)$summary[4:7,1]
  tmp <- data.frame(par = names(tmp), value = tmp)
  tmp <- separate(tmp, par, into = c("par", "condition"), sep = "\\[")
  tmp <- dcast(tmp, condition ~ par, value.var = "value")
  tmp$exp <- x
  tmp$condition <- as.character(factor(tmp$condition, labels = c("Control", "Experimental")))
  tmp$mu <- tmp$mu * sigma + mu
  tmp$sigma <- tmp$sigma * sigma
  tmp
}
tmp <- as.data.table(do.call("rbind", lapply(qs, bayesCon)))
tmp[,idx:=1:.N]
tmp <- tmp[,list(x = x <- getX(mu,sigma,steps),
                            F.x = pnorm(x, mu, sigma),
                            exp,
                            condition,
                            mu,
                            sigma),
                     by=idx]
tmp$Method <- "Bayes"
dat.con <- rbind(dat.con, tmp)
# end bayes                     

dat.con$exp <- as.character(factor(dat.con$exp, labels = 1:8))
dat.dists$exp <- as.character(factor(dat.dists$exp, labels = 1:8))

endpoints <- dat.con[,list(low=quantile(x, .02), 
                            high=quantile(x, .98)),by=exp]
setkey(dat.dists, exp)
setkey(endpoints, exp)
dat.CDF <- foreach(e = endpoints$exp, .combine='rbind') %do% {
  dat.dists[e, ][x>endpoints[e,low] & x <endpoints[e,high],]  
}
dat.CDF[,model:='raw']

p.CDF <- ggplot(aes(x=x, y=F.x, group=idx), data=dat.CDF[Method == "Bayes",]) + 
         # geom_line(alpha=.5, color="grey") + 
         geom_line(color="grey") + 
         labs(x="x", y="F(x)") +
         facet_wrap(~exp, scales="free", nrow = 2)

p.con <- p.CDF +
         geom_line(aes(x=x, y=F.x, 
                      group=condition, color=condition), 
                      data=dat.con[Method == "Bayes",], size = .5) +
         scale_color_brewer(type = "qual", palette = 6) +
         theme_classic(base_size = 11) +
         theme(panel.spacing = unit(.5, "lines"),
               panel.spacing.x = NULL, 
               panel.spacing.y = NULL,
               legend.position = "top", legend.title = element_blank()) +
         scale_y_continuous(expand = c(0, 0)) +
         scale_x_continuous(expand = c(0, 0))

p.CDF.ml <- ggplot(aes(x=x, y=F.x, group=idx), data=dat.CDF[Method == "Median",]) + 
         # geom_line(alpha=.5, color="grey") + 
         geom_line(color="grey") + 
         labs(x="x", y="F(x)") +
         facet_wrap(~exp, scales="free", nrow = 2)

p.con.ml <- p.CDF.ml +
         geom_line(aes(x=x, y=F.x, 
                      group=condition, color=condition), 
                      data=dat.con[Method == "Median",], size = .5) +
         scale_color_brewer(type = "qual", palette = 6) +
         theme_classic(base_size = 11) +
         scale_y_continuous(expand = c(0, 0)) +
         scale_x_continuous(expand = c(0, 0)) +
         theme(panel.spacing = unit(.5, "lines"),
               panel.spacing.x = NULL, 
               panel.spacing.y = NULL,
               legend.position = "top", legend.title = element_blank())

p.con + labs(title = "Hierarchical Bayesian")
ggsave(file = "figure/hb_cdf.pdf", width = 14.5, height = 3)
p.con.ml + labs(title = "Maximum Likelihood")
ggsave(file = "figure/ml_cdf.pdf", width = 14.5, height = 3)

### pooling coefficients

load("q4_pooling.RData")

tmp <- melt(as.array(sfit))
subj_sig <- tmp %>%
  filter(grepl("sigma_mu\\[1\\]", parameters)) %>%
  rename(sigma = value) %>%
  dplyr::select(iterations, chains, sigma)
mu <- tmp %>%
  filter(grepl("e_mu", parameters)) %>%
  left_join(subj_sig) %>%
  group_by(iterations, chains) %>%
  summarise(lambda = pmin((sd(value)/unique(sigma)) ^ 2, 1))

subj_sig <- tmp %>%
  filter(grepl("sigma_sigma\\[1\\]", parameters)) %>%
  rename(sigma = value) %>%
  dplyr::select(iterations, chains, sigma)
sigma <- tmp %>%
  filter(grepl("e_sigma", parameters)) %>%
  left_join(subj_sig) %>%
  group_by(iterations, chains) %>%
  summarise(lambda = pmin((sd(value)/unique(sigma)) ^ 2, 1))

mu$coef <- "mu"
sigma$coef <- "sigma"
tmp <- rbind(mu, sigma)

facet_names <- c("mu" = expression(mu),
                 "sigma" = expression(sigma))

pool.plot <- ggplot(tmp, aes(lambda)) +
  geom_histogram() +
  facet_wrap(~coef, labeller = label_parsed) +
  theme_classic(base_size = 11) +
  xlab("Pooling Coefficient") +
  ylab("Frequency") +
  scale_y_continuous(expand = c(0, 0))

### Supplement 1

# ML fits
dat.joe$exp <- factor(dat.joe$exp, labels = 1:8)

ml.plot <- ggplot(dat.joe, aes(sse)) +
  geom_histogram() +
  theme_classic() +
  facet_wrap(~exp, nrow = 2) +
  labs(x = "Sum of squared error", y = "Count")

ml.rsq <- mean(dat.joe$R.sq > .99)

# HBM fits

bayesRHat <- function(x) {
  load(paste0("q", x, ".RData"))
  tmp <- summary(sfit)$summary[,10]
  data.frame(rhat = tmp, exp = x)
}

tmp <- as.data.table(do.call("rbind", lapply(qs, bayesRHat)))
tmp$exp <- factor(tmp$exp, labels = 1:8)

hbm.rhat <- ggplot(tmp, aes(rhat)) +
  geom_histogram() +
  facet_wrap(~exp, nrow = 2) +
  theme_classic() +
  labs(x = "Rhat", y = "Count")

### Supplement sim

load("qSupp2.RData")

tmp <- melt(as.array(sfit))
tmp <- drop.levels(tmp[!grepl("idx|y|alpha|beta|lp", tmp$parameters),])
levels(tmp$parameters)[6:9] <- c("tau[1]", "tau[2]", "gamma[1]", "gamma[2]")

pars <- data.frame(parameters = levels(tmp$parameters),
                   value      = c(mu, phi, tau, gamma, sigma))

tmp <- tmp %>%
  group_by(parameters) %>%
  summarise(l95 = quantile(value, probs = c(.025)),
            l80 = quantile(value, probs = c(.1)),
            h80 = quantile(value, probs = c(.9)),
            h95 = quantile(value, probs = c(.975)),
            value = quantile(value, probs = c(.5)))

sim.plot <- ggplot(tmp, aes(parameters, value)) +
  geom_point() +
  geom_errorbar(width = 0, aes(ymin = l95, ymax = h95)) +
  geom_errorbar(width = 0, size = 1.5, aes(ymin = l80, ymax = h80)) +
  geom_point(data = pars, color = "red") +
  theme_classic()

# source and save, loaded in paper .Rnw
save.image(file = "paper.RData")