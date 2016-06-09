library("lme4")
library("ggplot2")
library("dplyr")
library("tidyr")
library("reshape2")
library("rstan")

###### First training study data #####

fig1 <- read.csv("Training20130204.csv", na.strings = ".") %>%
  mutate(idx = 1:n()) %>%
  gather(measure, value, -idx) %>%
  mutate(measure = gsub("([0-9]+)","~\\1", measure)) %>%
  separate(measure, into = c("measure", "session"), sep = "~") %>%
  group_by(measure, session) %>%
  summarise(value = mean(value, na.rm = TRUE)) %>%
  ggplot(aes(session, value)) +
    geom_point() +
    geom_line(size = .75, group = 1) +
    facet_wrap(~measure, scale = "free_y", ncol = 2) +
    theme_bw(base_size = 14) +
    guides(size = F) +
    ylab("Average Level of Achievement") +
    xlab("Session") +
    theme(panel.grid.major.x = element_blank(), 
          panel.grid.major.y = element_blank())

# ggsave(plot = fig1, file = "study1.png", height = 8, width = 8)

##### Study 2 #####

condFun <- function(x) {
  if (grepl("6", x)) {
    return("combo")
  } else {
      if (grepl("13", x)) {
        if(grepl("floop|nback", x)) {
          return("interference")
        } else {
          return("spatial")
        }
      }
      "control"
  }
}

readFun <- function(x) {
  tmp <- read.csv(x)
  tmp$condition <- condFun(x)
  tmp$task <- strsplit(x, split = "\\.|6|1")[[1]][1]
  tmp$subject <- factor(tmp$subject)
  melt(tmp)
}

dat <- do.call("rbind", lapply(list("floop6.csv",
                                    "floop13.csv",
                                    "nback6.csv",
                                    "nback13.csv",
                                    "mem6.csv",
                                    "mem13.csv",
                                    "shape6.csv",
                                    "shape13.csv",
                                    "follow.csv",
                                    "remember.csv"),
                                    readFun))

dat$scale.value <- unlist(with(dat, tapply(value, list(condition, task),
  FUN = function(x) (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE))))
dat$task <- factor(dat$task, labels = c("floop","follow me","memnosyne","n-back","remember me","shapebuilder"))
dat$variable <- factor(dat$variable, labels = as.character(1:length(unique(dat$variable))))

##### Hours in training #####

time <- dat %>%
  filter(!is.na(value)) %>%
  dplyr::count(subject, condition) %>%
  dplyr::mutate(time = n * ifelse(condition == "combo", 6.5/60, 13/60)) %>%
  group_by(condition) %>%
  summarise(trained = mean(time))

##### Plotting #####

# Figure 2

fig2 <- dat %>%
  filter(!variable %in% as.character(36:41)) %>%
  mutate(variable = as.numeric(variable)) %>%
  group_by(condition, task, variable) %>%
  summarise(value = mean(value, na.rm = TRUE)) %>%
  ggplot(aes(variable, value)) +
    geom_point() +
    geom_line(size = .75, group = 1) +
    facet_wrap(~ condition + task, scale = "free_y", ncol = 2,
      labeller = label_wrap_gen(multi_line=FALSE)) +
    ylab("Average Level of Achievement") +
    xlab("Session") +
    guides(size = F) +
    theme_bw(base_size = 14) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank())

# ggsave(plot = fig2, filename = "all.png", height = 8, width = 6)

##### Pre/post data #####

pdat <- read.csv("CogAssessmentDatasheetJan28_2012.csv",
            na.strings = ".")
names(pdat) <- tolower(names(pdat))
pdat <- pdat %>%
          dplyr::select(participant, completed, condition, nfcscore,
                      rspanscore, t2rspanscore, rspantotal, t2rspantotal,
                      shapescore, t2shapescore, ravens, t2ravens) %>%
          gather(measure, value, -participant, -completed, -condition, -nfcscore) %>%
          mutate(time = ifelse(grepl("t2", measure), "2", "1"),
                 measure = gsub("t2", "", measure)) %>%
          spread(measure, value)

##### Bayesian Multilevel Model #####

smod <- "
  data {
    int N;
    int K;
    int P;
    vector[N] y;
    matrix[N,P] x;
    int idx[N];
  }
  parameters {
    real alpha;
    vector[K] alpha_idx;
    vector[P] beta;
    real<lower=0> sigma_y;
    real<lower=0> sigma_idx;
  }
  transformed parameters {
    vector[N] yhat;
    for (n in 1:N)
      yhat[n] <- alpha_idx[idx[n]] + dot_product(beta, x[n]);
  }
  model {
    alpha ~ normal(0,1);
    beta ~ normal(0,1);
    sigma_y ~ inv_chi_square(2);
    sigma_idx ~ inv_chi_square(2);
    alpha_idx ~ normal(alpha, sigma_idx);
    y ~ student_t(4, yhat, sigma_y);
  }
"

tmp <- subset(pdat, completed == "1" & complete.cases(ravens))
tmp$condition <- factor(tmp$condition)

x <- model.matrix(~ -1 + condition + condition:time, data = tmp)

sdat <- list(
  y = scale(tmp$ravens)[,1],
  x = x,
  idx = unclass(factor(tmp$participant)),
  N = nrow(tmp),
  K = length(unique(tmp$participant)),
  P = ncol(x)
)

fit <- stan(model_code = smod, data = sdat)
# print(fit, probs = c(.025, .975),
#   pars = c("alpha", "beta", "sigma_idx", "sigma_y"))