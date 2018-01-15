library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(rstan)
library(psy)
library(rstanarm)

fesc <- read.csv("fescdata.csv")
fesc2 <- read.csv("fesc2data.csv")

# fesc cleaning

fesc$precort <- fesc$PRECORT
fesc$postcort <- fesc$POSTCORT
fesc$prestai <- fesc$PRESTAI
fesc$poststai <- fesc$POSTSTAI
fesc$cortchange <- fesc$postcort - fesc$precort
fesc$staichange <- fesc$poststai - fesc$prestai

# fesc2 cleaning

# STAI in study 2 was coded as sum rather than mean, thus the correction for number of items

fesc2$precort <- fesc2$PRECORT
fesc2$postcort <- fesc2$POSTCORT
fesc2$prestai <- fesc2$PRESTAITOT/10
fesc2$poststai <- fesc2$POSTSTAITOT/10
fesc2$cortchange <- fesc2$postcort - fesc2$precort
fesc2$staichange <- fesc2$poststai - fesc2$prestai

fesc$study <- "Study 1"
fesc2$study <- "Study 2"
get.var <- c("prestai", "poststai", "precort", "postcort", "staichange", "cortchange", "study")
graph.data <- rbind(fesc[,get.var], fesc2[,get.var])
graph.data$study <- factor(graph.data$study)
graph.data <- graph.data[complete.cases(graph.data),]

###### results

# figure 1

fescscatter <- ggplot(graph.data, aes(staichange, cortchange), color = "black") +
  geom_point() +
  geom_smooth(method = "lm", color = "black") +
  facet_wrap(~study) +
  xlab("Change in STAI") +
  ylab("Change in cortisol") +
  theme_bw()

# table 1

fes1_mod <- '
data {
  int<lower=1> N;
  vector[N] y;
  vector[N] x;
}
parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
  real<lower=1> nu;
}
model {
  beta ~ double_exponential(0,1);
  nu ~ gamma(2, 0.1);
  y ~ student_t(nu, alpha + beta * x, sigma);
}
'

fesc <- fesc[complete.cases(fesc$staichange, fesc$cortchange),]

fes1_dat <- list(N = nrow(fesc),
                 y = fesc$cortchange,
                 x = fesc$staichange)

fit.1 <- stan(model_code = fes1_mod, data = fes1_dat)

# table 2

tab2 <- graph.data %>%
  gather(variable, value, -study) %>%
  group_by(study, variable) %>%
  summarise(means = mean(value, na.rm = TRUE),
            sds   = sd(value, na.rm = TRUE))

# table 3

a.1 <- cor.test(~ prestai + poststai, data = fesc)
a.2 <- cor.test(~ prestai + poststai, data = fesc2)
a.c <- cor.test(~ prestai + poststai, data = graph.data)

s.1 <- cor.test(~ precort + postcort, data = fesc)
s.2 <- cor.test(~ precort + postcort, data = fesc2)
s.c <- cor.test(~ precort + postcort, data = graph.data)

c.1 <- cor.test(~ staichange + cortchange, data = fesc)
c.2 <- cor.test(~ staichange + cortchange, data = fesc2)
c.c <- cor.test(~ staichange + cortchange, data = graph.data)

# table 4

get.p <- as.array(fit.1)

prior.mean <- median(get.p[,,"beta"])
prior.var <- var(c(get.p[,,"beta"]))

fes2_mod <- '
data {
  int<lower=1> N;
  vector[N] y;
  vector[N] x;
  real prior_mean;
  real prior_var;
}
parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
  real<lower=1> nu;
}
model {
  beta ~ normal(prior_mean,prior_var);
  nu ~ gamma(2, 0.1);
  y ~ student_t(nu, alpha + beta * x, sigma);
}
'

fesc2 <- fesc2[complete.cases(fesc2$staichange, fesc2$cortchange),]

fes2_dat <- list(N = nrow(fesc2),
                 y = fesc2$cortchange,
                 x = fesc2$staichange,
                 prior_mean = prior.mean,
                 prior_var = prior.var)

fit.2 <- stan(model_code = fes2_mod, data = fes2_dat)

# multilevel model & vis

fit.3 <- stan_lmer(cortchange ~ staichange + (1 + staichange | study), data = graph.data)

graph.data$Var2 <- seq(nrow(graph.data))

g2 <- posterior_predict(fit.3, draws = 20) %>%
  reshape2::melt(.) %>%
  mutate(Var2 = unclass(factor(Var2))) %>%
  left_join(dplyr::select(graph.data, staichange, study, Var2))

predictscatter <- ggplot(g2, aes(staichange, value)) +
  geom_smooth(aes(group = Var1), se = FALSE, method = "lm", alpha = .2, color = "gray") +
  geom_point(data = graph.data, aes(staichange, cortchange)) +
  facet_wrap(~study) +
  theme_classic() +
  labs(y = "Change in cortisol", x = "Change in STAI")

# in-text

cb1 <- cronbach(fesc[,8:17])$alpha
cb2 <- cronbach(fesc2[,7:16])$alpha

cort.age <- lm(precort ~ AGE, data = fesc)
postcort.age <- lm(postcort ~ AGE, data = fesc)

stai.change.1 <- with(fesc, t.test(prestai, poststai, paired = TRUE))
cort.change.1 <- with(fesc, t.test(precort, postcort, paired = TRUE))
change.1 <- cor.test(~ staichange + cortchange, data = fesc)

cort.sex.desc <- fesc2 %>%
  select(GENDER, precort, postcort) %>%
  gather(variable, value, -GENDER) %>%
  group_by(GENDER, variable) %>%
  summarise(means = mean(value, na.rm = TRUE),
            sds = sd(value, na.rm = TRUE))
cort.sex.2 <- t.test(precort ~ GENDER, data = fesc2, var.equal = TRUE)

cort.sex.2p <- t.test(postcort ~ GENDER, data = fesc2, var.equal = TRUE)

stai.change.2 <- with(fesc2, t.test(prestai, poststai, paired = TRUE))
cort.change.2 <- with(fesc2, t.test(precort, postcort, paired = TRUE))
change.2 <- cor.test(~ staichange + cortchange, data = fesc2)

save.image("fesc.RData")