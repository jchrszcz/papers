# optional: get gemmR
library("devtools")
install_github("jchrszcz/gemmR", subdir = "gemmR")

# load required packages
library("MASS")
library("plyr")
library("reshape2")
library("ggplot2")
library("bestglm")
library("Hmisc")
library("gemmR")
library("parallel")

# defaults

reps <- 500
chains <- 10
gens <- 10
mccount <- 12
MC_CORES <- 12

# data

dat1 <- read.csv("Henry 2009.csv")
dat1$s.gini <- scale(dat1$gini_usethis)
dat1$s.pasture <- scale(dat1$percent_pasture)
dat1$s.gnp <- scale(dat1$GNPpercapita)
dat1$log.mr <- log(dat1$murderrate) + abs(min(log(dat1$murderrate)))
dat1$sqrt.mr <- sqrt(dat1$murderrate)
dat1$sqrt.gini <- sqrt(dat1$gini_usethis)
dat1$log.gini <- log(dat1$gini_usethis)
dat1$sqrt.gnp <- sqrt(dat1$GNPpercapita)
dat1$log.gnp <- log(dat1$GNPpercapita)
dat1$sqrt.pasture <- sqrt(dat1$percent_pastures)
dat1$log.pasture <- log(dat1$percent_pastures)

dat2 <- read.csv("IATDataSiegel.csv")
dat2 <- as.data.frame(scale(dat2))
dat2$f.atbmean <- factor(dat2$atbmean)

# outlier methods

outlier.mod <- lm(atbmean ~ raceampdiff + pollampdiff + RaceIAT + poliatDmean + emsmean + imsmean + pamean + Rstroopeffect + ssrt, data = dat2)

dat2$full <- 1
dat2$univariate <- ifelse(abs(dat2$atbmean) > 3, 0, 1)
dat2$cooks <- ifelse(cooks.distance(outlier.mod) > 4/nrow(dat2), 0, 1)
dat2$dffits <- ifelse(dffits(outlier.mod) > 2 * sqrt((length(coef(outlier.mod)) - 1)/nrow(dat2)), 0, 1)

race.select <- list(full = dat2$full,
  univariate = dat2$univariate,
  cooks = dat2$cooks,
  dffits = dat2$dffits)

coh.select <- list(full = rep(1, times = nrow(dat1)))

# model definitions

race.mods <- list(mod1 = atbmean ~ raceampdiff + pollampdiff + RaceIAT + poliatDmean + emsmean + imsmean + pamean + Rstroopeffect + ssrt)

coh.mods <- list(full = murderrate ~ percent_pastures + gini_usethis + GNPpercapita,
  log.crit = log.mr ~ percent_pastures + gini_usethis + GNPpercapita,
  log.all = log.mr ~ log.pasture + log.gini + log.gnp,
  sqrt.crit = sqrt.mr ~ percent_pastures + gini_usethis + GNPpercapita,
  sqrt.all = sqrt.mr ~ sqrt.pasture + sqrt.gini + sqrt.gnp)

# saturated models

sat.list <- list(false = FALSE)

# helper functions

mikeBic <- function(r, k, n) {
  return(n * log(1 - r^2) + k * log(n))
}

bicTau <- function(tau, p, k, n) {
  rp <- sin(pi/2 * tau * (n - p - 1)/n)
  return(n * log(1 - rp^2) + k * log(n))
}

# data generation functions

selectData <- function(dat, cv, ...) {
  if (cv) {
    picker <- sample(1:nrow(dat), nrow(dat)/2, replace = FALSE)
    return(list(est = dat[picker,], cross = dat[-picker,]))
  } else {
    return(list(est = dat))
  }
}

genModRace <- function(mod, dat, sat, ...) {
  ols <- lm(mod, data = dat)
  out <- list(ols = ols,
    ord.log = polr(update(mod, factor(.) ~ .), dat, Hess = TRUE),
    huber = rlm(mod, dat, psi = psi.huber),
    bisquare = rlm(mod, dat, psi = psi.bisquare),
    hampel = rlm(mod, dat, psi = psi.hampel)
  )
  if (sat) {
    out$gemmr <- gemm(mod, data = dat, fit.metric = "tau", n.chains = chains, n.gens = gens)
  } else {
    out$gemmr <- gemm(mod, data = dat, n.chains = chains, n.gens = gens)
    out$ols.bic <- bestglm(cbind(ols$model[,-1], ols$model[,1,drop=FALSE]))
    out$ols.bic$BestModel$ols <- ols
  }
  return(out)
}

genModCoh <- function(mod, dat, sat, ...) {
  ols <- lm(mod, data = dat)
  out <- list(ols = ols,
    ols.bic = bestglm(cbind(ols$model[,-1], ols$model[,1,drop=FALSE])),
    gemmr = gemm(mod, data = dat, n.chains = chains, n.gens = gens),
    poiss = glm(mod, data = dat, family = "poisson"),
    qpoiss = glm(mod, data = dat, family = "quasipoisson"),
    nb = glm.nb(mod, data = dat)
  )
  out$ols.bic$BestModel$ols <- ols 
  return(out)
}

# fitting functions

getFits <- function(mod, sat, cross, ...) {
  if (class(mod)[1] == "polr") {
    p <- nrow(coef(summary(mod)))
    k <- sum(ifelse(abs(coefficients(summary(mod))[,3]) > 1.96, 1, 0), na.rm = TRUE)
    if (!sat)
      mod$coefficients <- ifelse(abs(coefficients(summary(mod))[1:length(coef(mod)),3]) > 1.96, mod$coefficients, 0)
    if (k == 0) {
      if (is.null(cross) & missing(cross)) {
        return(rep(NA, 5))
      } else {
        return(rep(NA, 7))
      }
    }
    r <- cor(as.numeric(as.character(mod$model[,1])), as.numeric(as.character(predict(mod, type = "class"))))
    tau <- cor(as.numeric(as.character(mod$model[,1])), as.numeric(as.character(predict(mod, type = "class"))), method = "kendall")
  } else {
    if (class(mod)[1] == "bestglm") {
      p <- length(mod$BestModel$ols$coefficients) - 1
      k <- length(coef(mod$BestModel)) - 1
      mod <- mod$BestModel
    } else if (class(mod)[1] == "gemm") {
      p <- ncol(mod$model) - 1
      k <- sum(coef(mod)[1,-1] != 0, na.rm = TRUE)
    } else {
      p <- ncol(mod$model) - 1
      k <- sum(abs(coefficients(summary(mod)))[-1,3] > 1.96, na.rm = TRUE)
      if (!sat)
        mod$coefficients[-1] <- ifelse(abs(coefficients(summary(mod)))[-1,3] > 1.96, mod$coefficients[-1], 0)
    }
    if (k == 0) {
      if (is.null(cross) | missing(cross)) {
        return(rep(NA, 5))
      } else {
        return(rep(NA, 7))
      }
    }
    r <- cor(mod$model[,1], predict(mod, type = "response"))
    tau <- cor(mod$model[,1], predict(mod, type = "response"), method = "kendall")
  }
  if (sat) {
    k <- p
  }
  bict <- bicTau(tau, p, k, length(predict(mod)))
  bic <- mikeBic(r, k, length(predict(mod)))
  if(!is.null(cross) & !missing(cross)) {
    if (class(mod)[1] == "polr") {
      cross.r <- cor(cross[,names(cross) == unlist(strsplit(names(mod$model)[1], "[[:punct:]]"))[2]], as.numeric(as.character(predict(mod, cross, type = "class"))))
      cross.tau <- cor(cross[,names(cross) == unlist(strsplit(names(mod$model)[1], "[[:punct:]]"))[2]], as.numeric(as.character(predict(mod, cross, type = "class"))), method = "kendall")
    } else if (names(mod$model)[1] == "y") {
      cross.r <- cor(cross[,names(cross) == names(mod$ols$model)[1]], predict(mod, cross, type = "response"))
      cross.tau <- cor(cross[,names(cross) == names(mod$ols$model)[1]], predict(mod, cross, type = "response"), method = "kendall")
    } else {
      cross.r <- cor(cross[,names(cross) == names(mod$model)[1]], predict(mod, cross, type = "response"))
      cross.tau <- cor(cross[,names(cross) == names(mod$model)[1]], predict(mod, cross, type = "response"), method = "kendall")
    }
    return(c(bict= bict, bic = bic, tau = tau, r = r, k = k, cross.tau = cross.tau, cross.r = cross.r))
  } else {
    return(c(bict= bict, bic = bic, tau = tau, r = r, k = k))
  }
}

pRec <- function(mod, ...) {
  if (any(class(mod) == "gemm")) {
    rec <- ifelse(coef(mod)[1,-1] != 0, 1, 0)
  }  else if (any(class(mod) == "polr")) {
    rec <- ifelse(abs(coefficients(summary(mod))[1:length(coef(mod)),3]) > 1.96, 1, 0)
  }  else if (any(class(mod) == "bestglm")) {
    rec <- unlist(mod$BestModels[mod$BestModels$Criterion == min(mod$BestModels$Criterion), -ncol(mod$BestModels)])
  } else {
    rec <- ifelse(abs(coefficients(summary(mod))[-1,3]) > 1.96, 1, 0)
  }
  return(rec)
}

getCoef <- function(mod, ...) {
    if (any(class(mod) == "gemm")) {
    co <- coef(mod)[1,-1]
  }  else if (any(class(mod) == "polr")) {
    co <- coef(mod)
  }  else if (any(class(mod) == "bestglm")) {
    co <- unlist(mod$BestModels[mod$BestModels$Criterion == min(mod$BestModels$Criterion), -ncol(mod$BestModels)])
    co[co] <- coef(mod$BestModel)[-1]
  } else {
    co <- coef(mod)[-1]
  }
  return(co)
}

allSums <- function(mod, sat, cross, ...) {
  out <- list(fits = getFits(mod, sat, cross, ...))  
  out$rec <- pRec(mod, ...)
  out$coefs <- getCoef(mod, ...)
  return(out)
}

# controller functions

genModHelp <- function(mods, sat, cross, ...) {
  out <- lapply(mods, allSums, sat, cross, ...)
  return(out)
}

genMods <- function(dat.select, mods, dat, pass.fun, cv, sat, ...) {
  dat <- selectData(dat[dat.select == 1,], cv)
  mods <- lapply(mods, pass.fun, dat$est, sat, ...)
  out <- lapply(mods, genModHelp, sat, dat$cross, ...)
  return(out)
}

trimData <- function(sat, trim.list, mods, dat, pass.fun, cv, ...) {
  return(lapply(trim.list, genMods, mods, dat, pass.fun, cv, sat, ...))
}

satMods <- function(trash, sat.list, trim.list, mods, dat, pass.fun, cv, ...) {
  return(lapply(sat.list, trimData, trim.list, mods, dat, pass.fun, cv, ...))  
}

# commands

full.race <- satMods(NA, list(FALSE), race.select, race.mods, dat2, genModRace, FALSE, chains = chains, gens = gens)
half.race <- mclapply(1:reps, satMods, sat.list, list(race.select[[1]]), race.mods, dat2, genModRace, TRUE, chains = chains, gens = gens, mc.cores = mccount)
full.coh <- satMods(NA, list(FALSE), list(rep(1, nrow(dat1))), coh.mods, dat1, genModCoh, FALSE, chains = chains, gens = gens)
half.coh <- mclapply(1:reps, satMods, list(FALSE), list(rep(1, nrow(dat1))), coh.mods, dat1, genModCoh, TRUE, chains = chains, gens = gens, mc.cores = mccount)

##### munging

# table1

temp <- melt(full.race)
temp[-1] <- lapply(temp[-1], factor)
levels(temp$L1) <- c("false", "true")
temp$value <- round(temp$value, digits = 3)

tab1 <- temp[temp$L5 == "fits",]
tab1$measure <- c("bict", "bic", "tau", "r", "k")
table1 <- acast(tab1, L1 + L4 + L2 ~ measure)[,c(2,1,5,4,3)]

# table2

tab2 <- temp[temp$L5 == "coefs" & temp$L1 == "false",]
tab2$coef.name <- attr(terms(race.mods[[1]]), "term.labels")
table2 <- acast(tab2,L4 + L2 ~ coef.name)

# table3

temp <- melt(half.race)
temp[-1] <- lapply(temp[-1], factor)


tab3 <- temp[temp$L6 == "rec" & temp$L2 == "false",]
tab3$coefs <- factor(colnames(attr(terms(race.mods[[1]]), "factors")))
table3 <- acast(tab3, L5 ~ coefs, function(x) sum(x, na.rm = TRUE)/reps)


# table4

tab4 <- temp[temp$L6 == "fits",]
tab4$measure <- c("bict", "bic", "tau", "r", "k", "cross.tau", "cross.r")
table4 <- acast(tab4, L2 + L5 ~ measure, mean, na.rm = TRUE)[,c(2,1,7,6,5,4,3)]

# table5

temp <- melt(full.coh)
temp[-1] <- lapply(temp[-1], factor)
temp$value <- round(temp$value, digits = 3)

tab5 <- temp[temp$L5 == "fits",]
tab5$measure <- c("bict", "bic", "tau", "r", "k")
table5 <- acast(tab5, L4 + L3 ~ measure)[,c(2,1,5,4,3)]

# table6

temp <- melt(half.coh)
temp[-1] <- lapply(temp[-1], factor)

tab6 <- temp[temp$L6 == "rec",]
tab6$coefs <- factor(colnames(attr(terms(coh.mods[[1]]), "factors")))
table6 <- acast(tab6, L5 + L4 ~ coefs, function(x) sum(x, na.rm = TRUE)/reps)

# table7

tab7 <- temp[temp$L6 == "fits",]
tab7$measure <- c("bict", "bic", "tau", "r", "k", "cross.tau", "cross.r")
table7 <- acast(tab7, L4 + L5 ~ measure, mean, na.rm = TRUE)[,c(2,1,7,6,5,4,3)]
table7 <- apply(table7, c(1,2), round, digits = 3)

# bootstrap standard errors

race.mods2 <- list(mod1 = atbmean ~ raceampdiff + emsmean + imsmean)

gemmSE <- function(dat.select, mod, dat) {
  dat <- dat[dat.select == 1,]
  new.dat <- dat[sample(1:nrow(dat), nrow(dat), replace = TRUE),]
  g.mod <- gemm(mod, new.dat, n.chains = 10, n.gens = 10)
  return(coef(g.mod))
}

apGemm <- function(trash, select.list, mod, dat) {
  return(lapply(select.list, gemmSE, mod, dat))
}

gemm.se <- mclapply(1:reps, apGemm, race.select, race.mods2[[1]], dat2, mc.cores = mccount)

temp <- melt(gemm.se)
gemm.sterr <- ddply(temp[temp$Var1 == 1,], .(Var2, L2), summarise, means = mean(value, na.rm = TRUE), sds = sd(value, na.rm = TRUE), q.025 = quantile(value, probs = .025, na.rm = TRUE), q.0975 = quantile(value, probs = .975, na.rm = TRUE))

allSE <- function(dat.select, mod, dat) {
  dat <- dat[dat.select == 1,]
  new.dat <- dat[sample(1:nrow(dat), nrow(dat), replace = TRUE),]
  ols <- lm(mod, new.dat)
  out <- list(ols = ols,
    ord.log = polr(update(mod, factor(.) ~ .), new.dat, Hess = TRUE),
    huber = rlm(mod, new.dat, psi = psi.huber),
    bisquare = rlm(mod, new.dat, psi = psi.bisquare),
    hampel = rlm(mod, new.dat, psi = psi.hampel),
    ols.bic = bestglm(cbind(ols$model[,-1], ols$model[,1,drop=FALSE])))
  out$ols.bic$BestModel$ols <- ols
  real.out <- lapply(out, getCoef)
  return(real.out)
}

apAll <- function(trash, select.list, mod, dat) {
  return(lapply(select.list, allSE, mod, dat))
}

all.se <- lapply(1:reps, apAll, race.select, race.mods[[1]], dat2)

temp <- melt(all.se)
temp$coefs <- factor(colnames(attr(terms(race.mods[[1]]), "factors")))
all.sterr <- ddply(temp, .(coefs, L2, L3), summarise, ses = sd(value))
table.se <- acast(all.sterr, L3 + L2 ~ coefs)


### figures

# figure 1: bivariate scatter for race data with outcome

fig1.dat <- melt(dat2[names(dat2) %in% c("atbmean", "pollampdiff", "RaceIAT", "poliatDmean", "raceampdiff", "imsmean", "emsmean", "pamean", "Rstroopeffect", "ssrt")], id.vars = "atbmean")
levels(fig1.dat$variable) <- c("Political AMP", "Race IAT", "Political IAT", "Race AMP", "Internal MCP", "External MCP", "Political Attitude", "Stroop", "Stop Signal")
fig1 <- ggplot(fig1.dat, aes(x = value, y = atbmean)) + geom_point() + facet_wrap(~variable) + ylab("Attitude Towards Blacks") + xlab("") + theme_classic()

# figure 2: univariate hist for race data

fig2.dat <- melt(dat2[names(dat2) %in% c("atbmean", "pollampdiff", "RaceIAT", "poliatDmean", "raceampdiff", "imsmean", "emsmean", "pamean", "Rstroopeffect", "ssrt")])
levels(fig2.dat$variable) <- c("Attitude Towards Blacks", "Political AMP", "Race IAT", "Political IAT", "Race AMP", "Internal MCP", "External MCP", "Political Attitude", "Stroop", "Stop Signal")
fig2 <- ggplot(fig2.dat, aes(x = value)) + geom_histogram() + facet_wrap(~variable) + theme_classic()

# figure 3: bivariate scatter for CoH, normal and log outcome

fig3.dat <- melt(dat1[names(dat1) %in% c("murderrate", "percent_pastures", "gini_usethis", "GNPpercapita")], id.vars = "murderrate")
levels(fig3.dat$variable) <- c("Pastureland", "Gini", "Per Capita GDP")
fig3.dat <- rename(fig3.dat, c(murderrate = "Homicide"))
fig3a <- ggplot(fig3.dat, aes(x = value, y = Homicide)) + geom_point() + facet_wrap(~variable, scales = "free_x") + xlab("") + theme_classic()
fig3b <- ggplot(fig3.dat, aes(x = value, y = Homicide)) + geom_point() + facet_wrap(~variable, scales = "free_x") + ylab("log(Homicide)") + xlab("") + scale_y_log10() + theme_classic()

# figure 4: univariate hist for CoH

fig4.dat <- melt(dat1[names(dat1) %in% c("murderrate", "percent_pastures", "gini_usethis", "GNPpercapita")])
levels(fig4.dat$variable) <- c("Homicide", "Pastureland", "Gini", "Per Capita GDP")
fig4 <- ggplot(fig4.dat, aes(x = value)) + geom_histogram() + facet_wrap(~variable, scales = "free") + theme_classic()

# figure 5: half-sample precover for gemm and OLS mods (CoH), half-sample precover for gemm and GLM (CoH)

fig5.dat <- ddply(tab6, .(L5, L4, coefs), summarise, p.rec = sum(value)/length(value))
levels(fig5.dat$coefs) <- c("Gini", "GNP Per Capita", "Pastureland")
fig5a <- ggplot(fig5.dat[fig5.dat$L5 %in% c("gemmr", "ols", "ols.bic"),], aes(x = coefs, y = p.rec, fill = L5)) + geom_bar(stat = "identity", position = "dodge") + ylab(c("p(recovery)")) + facet_wrap(~L4) + theme_classic() + scale_fill_grey()
fig5b <- ggplot(fig5.dat[fig5.dat$L5 %in% c("gemmr", "poiss", "qpoiss", "nb") & fig5.dat$L4 == "full",], aes(x = coefs, y = p.rec, fill = L5)) + geom_bar(stat = "identity", position = "dodge") + ylab(c("p(recovery)")) + theme_classic() + scale_fill_grey()

fig5balt <- ggplot(fig5.dat[fig5.dat$L5 %in% c("gemmr", "poiss", "qpoiss", "nb") & fig5.dat$L4 == "full",], aes(x = L5, y = p.rec)) + geom_bar(stat = "identity", position = "dodge") + facet_grid(~coefs) + ylab(c("p(recovery)")) + theme_classic() + scale_fill_grey()

save(list = ls(), file = "meth_sims.RData")