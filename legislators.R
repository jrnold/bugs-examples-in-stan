## ----legislators_setup,message=FALSE-------------------------------------
library("pscl")
library("tidyverse")
library("forcats")
library("stringr")
library("rstan")
library("sn")

## ----legislators_identification------------------------------------------
xi <- c(-1, -0.5, 0.5, 1)
alpha <- c(1, 0, -1)
beta <- c(-0.5, 0, 0.5)
y <- matrix(c(1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1), 3, 4)
k <- 1

list(sum(plogis(y - (alpha + beta %o% xi))),
     sum(plogis(y - (alpha + -beta %o% -xi))),
     sum(plogis(y - ((alpha - beta * k) + beta %o% (xi + k)))),
     sum(plogis(y - ((alpha + (beta / k) %o% (xi * k))))))

## ------------------------------------------------------------------------
data("s109", package = "pscl")

## ------------------------------------------------------------------------
s109

## ------------------------------------------------------------------------
s109_vote_data <- as.data.frame(s109$vote.data) %>%
  mutate(rollcall = paste(session, number, sep = "-"),
         passed = result %in% c("Confirmed", "Agreed To", "Passed"),
         votestotal = yeatotal + naytotal,
         yea_pct = yeatotal / (yeatotal + naytotal),
         unanimous = yea_pct %in% c(0, 1),
         close = yea_pct < 0.35 | yea_pct > 0.65,
         lopsided = yea_pct < 0.025 | yea_pct > 0.975) %>%
  filter(!unanimous) %>%
  select(-unanimous) %>%
  mutate(.rollcall_id = row_number())

s109_legis_data <- as.data.frame(s109$legis.data) %>%
  rownames_to_column("legislator") %>%
  mutate(.legis_id = row_number(),
         party = fct_recode(party,
                            "Democratic" = "D",
                            "Republican" = "R",
                            "Independent" = "Indep"))

s109_votes <- s109$votes %>%
  as.data.frame() %>%
  rownames_to_column("legislator") %>%
  gather(rollcall, vote, -legislator) %>%
  # recode to Yea (TRUE), Nay (FALSE), or missing
  mutate(yea = NA,
         yea = if_else(vote %in% c(1, 2, 3), TRUE, yea),
         yea = if_else(vote %in% c(4, 5, 6), FALSE, yea)
         ) %>%
  filter(!is.na(yea)) %>%
  inner_join(dplyr::select(s109_vote_data, rollcall, .rollcall_id), by = "rollcall") %>%
  inner_join(dplyr::select(s109_legis_data, legislator, party, .legis_id), by = "legislator")

partyline <-
  s109_votes %>%
  group_by(.rollcall_id, party) %>%
  summarise(yea = mean(yea)) %>%
  spread(party, yea) %>%
  ungroup() %>%
  mutate(partyline = NA_character_,
         partyline = if_else(Republican < 0.1 & Democratic > 0.9,
                             "Democratic", partyline),
         partyline = if_else(Republican > 0.9 & Democratic < 0.1,
                             "Republican", partyline)) %>%
  rename(pct_yea_D = Democratic, pct_yea_R = Republican) %>%
  select(-Independent)

s109_vote_data <-
  left_join(s109_vote_data, partyline, by = ".rollcall_id")

## ----legislators_mod_ideal_point_1,results='hide',cache.extra=tools::md5sum("stan/ideal_point_1.stan")----
mod_ideal_point_1 <- stan_model("stan/ideal_point_1.stan")

## ------------------------------------------------------------------------
mod_ideal_point_1

## ----legislators_xi_init_1-----------------------------------------------
xi_1 <-
  s109_legis_data %>%
  mutate(
         xi = if_else(legislator == "FRIST (R TN)", 1,
                 if_else(legislator == "REID (D NV)", -1, NA_real_)),
         init = if_else(party == "Republican", 1,
                        if_else(party == "Democratic", -1, 0)))


## ----legislators_data_1--------------------------------------------------
legislators_data_1 <-
  within(list(), {
    y <- as.integer(s109_votes$yea)
    y_idx_leg <- as.integer(s109_votes$.legis_id)
    y_idx_vote <- as.integer(s109_votes$.rollcall_id)
    Y_obs <- length(y)
    N <- max(s109_votes$.legis_id)
    K <- max(s109_votes$.rollcall_id)
    # priors
    alpha_loc <- 0
    alpha_scale <- 5
    beta_loc <- 0
    beta_scale <- 2.5
    N_xi_obs <- sum(!is.na(xi_1$xi))
    idx_xi_obs <- which(!is.na(xi_1$xi))
    xi_obs <- xi_1$xi[!is.na(xi_1$xi)]
    N_xi_param <- sum(is.na(xi_1$xi))
    idx_xi_param <- which(is.na(xi_1$xi))
    tau_scale <- 5
    zeta_loc <- 0
    zeta_scale <- 10
  })

## ----legislators_init_1--------------------------------------------------
legislators_init_1 <- list(
  list(xi_param = xi_1$init[is.na(xi_1$xi)])
)

## ----message=FALSE,results='hide',warning=FALSE--------------------------
mod_ideal_point_1 <- stan_model("stan/ideal_point_1.stan")

## ----legislators_fit_1,results='hide'------------------------------------
legislators_fit_1 <-
  sampling(mod_ideal_point_1, data = legislators_data_1,
           chains = 1, iter = 500,
           init = legislators_init_1,
           refresh = 100,
           pars = c("alpha", "beta", "xi"))

## ----legislator_summary_1------------------------------------------------
legislator_summary_1 <-
  bind_cols(s109_legis_data,
           as_tibble(summary(legislators_fit_1, par = "xi")$summary)) %>%
  mutate(legislator = fct_reorder(legislator, mean))

## ----legislator_plot_1,fig.height=8,fig.width=4,fig.cap="Estimated Ideal Points of the Senators of the 109th Congress"----
ggplot(legislator_summary_1,
       aes(x = legislator, y = mean,
           ymin = `2.5%`, ymax = `97.5%`, colour = party)) +
  geom_pointrange() +
  coord_flip() +
  scale_color_manual(values = c(Democratic = "blue", Independent = "gray", Republican = "red")) +
  labs(y = expression(xi[i]), x = "", colour = "Party") +
  theme(legend.position = "bottom")

## ------------------------------------------------------------------------
map_df(c(-50, 0, 50),
        function(alpha) {
          tibble(x = seq(-4, 4, by = 0.1),
                 density = dsn(x, alpha = alpha),
                 alpha = alpha)
        }) %>%
  ggplot(aes(x = x, y = density, colour = factor(alpha))) +
  geom_line()

## ----legislators_mod_ideal_point_3,results='hide',cache.extra=tools::md5sum("stan/ideal_point_3.stan")----
mod_ideal_point_3 <- stan_model("stan/ideal_point_3.stan")

## ------------------------------------------------------------------------
mod_ideal_point_3

## ----legislators_data_2--------------------------------------------------
legislators_data_2 <-
  within(list(), {
    y <- as.integer(s109_votes$yea)
    y_idx_leg <- as.integer(s109_votes$.legis_id)
    y_idx_vote <- as.integer(s109_votes$.rollcall_id)
    Y_obs <- length(y)
    N <- max(s109_votes$.legis_id)
    K <- max(s109_votes$.rollcall_id)
    # priors
    alpha_loc <- 0
    alpha_scale <- 5
    beta_loc <- 0
    beta_scale <- 2.5
    xi_skew <- if_else(s109_legis_data$legislator == "FRIST (R TN)", 50, 0)
  })

## ----legislators_init_2--------------------------------------------------
legislators_init_2 <- function(chain_id) {
  list(xi_raw = if_else(s109_legis_data$party == "Republican", 1,
                    if_else(s109_legis_data$party == "Democratic", -1, 0)))
}

## ----legislators_fit_2,results='hide'------------------------------------
legislators_fit_2 <- sampling(mod_ideal_point_3,
                              data = legislators_data_2,
                              init = legislators_init_2,
                              chains = 1, iter = 500,
                              pars = c("alpha", "beta", "xi"))

## ----legislators_mod_ideal_point_2,results='hide',cache.extra=tools::md5sum("stan/ideal_point_2.stan")----
mod_ideal_point_2 <- stan_model("stan/ideal_point_2.stan")

## ------------------------------------------------------------------------
mod_ideal_point_2

## ----legislators_data_3--------------------------------------------------
legislators_data_3 <-
  within(list(), {
    y <- as.integer(s109_votes$yea)
    y_idx_leg <- as.integer(s109_votes$.legis_id)
    y_idx_vote <- as.integer(s109_votes$.rollcall_id)
    Y_obs <- length(y)
    N <- max(s109_votes$.legis_id)
    K <- max(s109_votes$.rollcall_id)
    # priors
    alpha_loc <- 0
    alpha_scale <- 5
    beta_loc <- rep(0, K)
    beta_scale <- rep(2.5, K)
    beta_skew <- if_else(s109_vote_data$rollcall == "2-169", 50, 0)
  })

## ----legislators_init_3--------------------------------------------------
legislators_init_3 <- function(chain_id) {
  list(beta = if_else(s109_vote_data$partyline %in% "Republican", 1,
                      if_else(s109_vote_data$partyline %in% "Democratic", -1,
                              0)),
       alpha = plogis(s109_vote_data$yea_pct),
       xi_raw = if_else(s109_legis_data$party == "Republican", 1,
                                      if_else(s109_legis_data$party == "Democratic", -1, 0)))
}

## ----legislators_fit_3,results='hide'------------------------------------
legislators_fit_3 <- sampling(mod_ideal_point_2,
                              data = legislators_data_3,
                              init = legislators_init_3,
                              chains = 1, iter = 500,
                              pars = c("alpha", "beta", "xi"))

