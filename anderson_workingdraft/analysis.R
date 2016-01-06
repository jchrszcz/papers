# required packages
library(dplyr)
library(ggplot2)
library(reshape2)
library(MASS)
library(tidyr)
library(lme4)
library(rstan)

##### Read in data and format for modeling
all.dat <- read.csv("data.csv")
names(all.dat) <- tolower(names(all.dat))

# factor and scale
all.dat$mite <- factor(all.dat$mindintheeyes)
all.dat$age <- scale(all.dat$age)[,1]
all.dat$fsiq <- scale(all.dat$fsiq)[,1]
all.dat$split.age <- factor(ifelse(all.dat$age < 0, 4, 6))

dat <- all.dat %>%
  filter(outlying_fa == 0, outlying_motion == 0) %>%
  mutate(unc_left = scale(lh.unc_avg_weight)[,1],
    unc_right = scale(rh.unc_avg_weight)[,1],
    ilf_left = scale(lh.ilf_avg_weight)[,1],
    ilf_right = scale(rh.ilf_avg_weight)[,1],
    cst_left = scale(lh.cst_avg_weight)[,1],
    cst_right = scale(rh.cst_avg_weight)[,1],
    leftamy = scale(leftamy)[,1],
    rightamy = scale(rightamy)[,1],
    unc_right_md = scale(rh.unc_md) [,1],
    unc_left_md =scale(lh.unc_md) [,1],
    ilf_left_md =scale(lh.ilf_md) [,1],
    ilf_right_md = scale(rh.ilf_md) [,1],
    cst_left_md =scale(lh.cst.md) [,1],
    cst_right_md = scale(rh.cst.md) [,1],
    rotation = scale(avg_rotation) [,1],
    translation = scale(avg_translation) [,1]
  )

null <- polr(mite ~ 1, data = dat)

fit1 <- polr(mite ~ fsiq + sex + translation + rotation + age * unc_left + cst_left + leftamy + cst_left:age + leftamy:age, data = dat)
fit2 <- polr(mite ~ fsiq + sex + translation + rotation + age * unc_right + cst_right + rightamy + cst_right:age + rightamy:age, data = dat)
fit3 <- polr(mite ~ fsiq + sex + translation + rotation + age * ilf_left + cst_left + leftamy + cst_left:age + leftamy:age, data = dat)
fit4 <- polr(mite ~ fsiq + sex + translation + rotation + age * ilf_right + cst_right + rightamy + cst_right:age + rightamy:age, data = dat)

# pseudo R^2

1 - logLik(fit1) / logLik(null)
1 - logLik(fit2) / logLik(null)
1 - logLik(fit3) / logLik(null)
1 - logLik(fit4) / logLik(null)

1 - logLik(review) / logLik(null)

1 - logLik(fit5) / logLik(null)
1 - logLik(fit6) / logLik(null)

##### Checking the continuous multilevel models

mdat <- dat %>%
  dplyr::select(id,
                mite,
                unc_left,
                age,
                split.age,
                fsiq,
                cst_left,
                unc_right,
                cst_right,
                ilf_left,
                ilf_right) %>%
  gather(area, value, -mite, -age, -fsiq, -id, -split.age) %>%
  separate(area, into = c("area", "hemisphere"), sep = "\\_")

mdat$hemisphere <- factor(mdat$hemisphere)
mdat$area <- factor(mdat$area)

mlm1 <- mdat %>%
  filter(hemisphere == "left", area != "ilf") %>%
  lmer(data = ., value ~ fsiq + age * mite * area + (1|id))

mlm2 <- mdat %>%
  filter(hemisphere == "right", area != "ilf") %>%
  lmer(data = ., value ~ fsiq + age * mite * area + (1|id))

mlm3 <- mdat %>%
  filter(hemisphere == "left", area != "unc") %>%
  lmer(data = ., value ~ fsiq + age * mite * area + (1|id))

mlm4 <- mdat %>%
  filter(hemisphere == "right", area != "unc") %>%
  lmer(data = ., value ~ fsiq + age * mite * area + (1|id))


mlm.fit <- lmer(value ~ hemisphere * area * unclass(mite) * age + (1 | id), data = mdat)

##### Followup ordinal correlation for subsets

with(dat[dat$split.age == "4",], cor.test(mindintheeyes, unc_left, method = "kendall"))

with(dat[dat$split.age == "6",], cor.test(mindintheeyes, unc_left, method = "kendall"))

##### Plots

ggplot(data = mdat, aes(x = value, y = mite, color = split.age)) +
  geom_jitter(position = position_jitter(height = .1, width = 0)) +
  facet_grid(area ~ hemisphere)

fig.dat <- all.dat %>%
  filter(outlying_fa == 0, outlying_motion == 0) %>%
  mutate(unc_left = lh.unc_avg_weight,
    unc_right = rh.unc_avg_weight,
    ilf_left = lh.ilf_avg_weight,
    ilf_right = rh.ilf_avg_weight,
    new.mite = as.numeric(as.character(dat$mite))
  ) %>%
  dplyr::select(id,
                new.mite,
                unc_left,
                split.age,
                fsiq,
                unc_right,
                ilf_left,
                ilf_right) %>%
  gather(area, value, -new.mite, -split.age, -fsiq, -id) %>%
  separate(area, into = c("area", "hemisphere"), sep = "\\_")

fig.dat$hemisphere <- as.factor(fig.dat$hemisphere)
fig.dat$area <- as.factor(fig.dat$area)

levels(fig.dat$hemisphere) <- c("Left", "Right")
levels(fig.dat$area) <- c("ILF", "UNC")

fig2 <- ggplot(fig.dat, aes(split.age, value, fill = split.age)) +
  geom_boxplot() +
  facet_grid(area ~ hemisphere) +
  scale_fill_brewer(type = "qual", palette = 2) +
  xlab("Age") +
  ylab("FA Weighted Average") +
  theme_classic() +
  guides(fill = FALSE)

fig3 <- ggplot(fig.dat, aes(x = value, y = new.mite, color = split.age)) +
  geom_jitter(position = position_jitter(height = .1, width = 0)) +
  geom_smooth(method = "lm") +
  facet_grid(area ~ hemisphere) +
  theme_classic() +
  xlab("FA Weighted Average") +
  ylab("SERT") +
  scale_color_brewer(type = "qual", palette = 2) +
  guides(color=guide_legend(title="Age"))

ggsave(filename = "~/Dropbox/DTI shared/Writing/Figures/fig2.png", plot = fig2, width = 4, height = 4)

ggsave(filename = "~/Dropbox/DTI shared/Writing/Figures/fig3.png", plot = fig3, width = 6, height = 4)

##### Bayesian analysis

# create predictor matrix
x <- model.matrix(~ fsiq + rotation + translation + sex * age * unc_left + cst_left + cst_left:age + cst_left:sex + cst_left:age:sex + leftamy + leftamy:age + leftamy:sex + leftamy:age:sex + unc_right + unc_right:sex + unc_right:age + unc_right:age:sex, dat)

bdat <- list(
  N = nrow(dat),
  P = ncol(x),
  K = length(levels(dat$mite)),
  y = unclass(dat$mite),
  x = x
)

# Stan code for ordered logistic model
bmod <- "
data {
  int<lower=1> N;
  int<lower=1> K;
  int<lower=1> P;
  int<lower=1,upper=K> y[N];
  row_vector[P] x[N];
}
parameters {
  vector[P] beta;
  ordered[K - 1] c;
}
model {
  beta ~ student_t(4.0, 0.0, 5.0/4.0);
  for (n in 1:N)
    y[n] ~ ordered_logistic(x[n] * beta, c);
}
"

bfit <- stan(model_code = bmod, data = bdat)