################################################################################
#
# Simulate longitudinal dataset with:
#
#   - gender          : Confounder
#   - SES             : Confounder
#   - PFOS/PFOA/PFHxS/PFNA : Exposure
#   - LOD per stof     : elke PFAS-stof heeft een vaste, vooraf gekozen detectie-
#                        limiet (PFOS_lod = 0.2, enz.). Voor elke stof <naam> komen
#                        er drie kolommen in de dataset:
#                          <naam>      : ruwe meting, met -3 als sentinelwaarde
#                                        voor metingen onder de LOD (in plaats van
#                                        de werkelijke -- onbekende -- waarde)
#                          <naam>_lod  : constante kolom met de gehanteerde LOD
#                          <naam>_imp  : STOCHASTISCHE imputatie op basis van
#                                        een gefitte, links-gecensureerde log-
#                                        normale verdeling (zie stap 4), i.p.v.
#                                        de -3 sentinel
#   - random intercept : subject-specifieke random intercept b_i ~ N(0, sigma_b^2)
#                        met sigma_b groot t.o.v. de residuele SD, voor een sterke
#                        within-subject correlatie (hoge ICC) tussen de 3 metingen
#   - missing data     : enkel op de UITKOMST (neutrofielenaantal), met VASTE
#                        steekproefgroottes per golf (301 op jaar 0, 220 op
#                        jaar 7, 200 op jaar 14), monotone uitval, en wie het
#                        eerst afhaakt gestuurd door een SES-gelinkt dropoutrisico
#
# Output: een long-format data.frame (1 rij per kind per meetmoment).

set.seed(20260821)


### Design parameters
#--------------------
n_subjects <- 301            
ages       <- c(0, 7, 14)
n_times    <- length(ages)

## PFAS-panel: vier blootstellingen, elk met een vaste, vooraf gekozen LOD.
exposure_names <- c("PFOS", "PFOA", "PFHxS", "PFNA")
LOD <- c(PFOS = 0.1, PFOA = 0.1, PFHxS = 0.2, PFNA = 0.2)

## SES: categorisch met 3 niveaus, met aandeel van de cohorte per niveau
SES_levels <- c("laag", "midden", "hoog")
SES_probs  <- c(0.30, 0.40, 0.30)

## Ware regressieparameters voor de uitkomst (neutrofielenaantal, 10^9/L)
beta0        <-  4.50   # intercept (gemiddeld neutrofielenaantal bij referentie)
beta_age     <- -0.03   # licht dalend neutrofielenaantal met leeftijd
beta_gender  <-  0.90   # STERK effect van gender op de uitkomst (confounding-pad 2)

## Effect per PFAS-stof op de uitkomst (log-schaal) -- bewust verschillend sterk:
beta_logPFOS  <- 0.65   # STERK effect
beta_logPFOA  <- 0.28   # ZWAK effect
beta_logPFHxS <- 0.00   # GEEN effect
beta_logPFNA  <- 0.00   # GEEN effect

## Effect van SES-niveau op de uitkomst, t.o.v. referentieniveau "laag"
beta_SES_midden <- -0.10
beta_SES_hoog   <- -0.25

## Random-intercept / residuele variantie -> sterke within-subject correlatie
sigma_b   <- 0.90       # SD van de random intercept (subject-niveau)
sigma_eps <- 0.45       # SD van het meetfoutresidu (meetmoment-niveau)
icc_true  <- sigma_b^2 / (sigma_b^2 + sigma_eps^2)

## Parameters voor het PFAS-blootstellingsmodel (lognormaal, per stof).
## gender is hier de confounder: jongens hebben systematisch hogere PFAS-
## waarden dan meisjes, voor alle vier de stoffen (met verschillende sterkte).
## De vier stoffen delen bovendien een subject-niveau "blootstellingsburden"
## (b_exposure), zodat ze onderling gecorreleerd zijn -- zoals typisch is
## voor een PFAS-mengsel -- met een stofspecifieke lading (burden_loading).
mu_logexp_female    <- c(PFOS = log(4.0), PFOA = log(2.5), PFHxS = log(0.2), PFNA = log(0.3))
gender_logexp_effect <- c(PFOS = 0.45,    PFOA = 0.30,     PFHxS = 0.25,     PFNA = 0.20)
age_logexp_effect    <- c(PFOS = -0.02,   PFOA = -0.015,   PFHxS = -0.01,    PFNA = -0.01)
sigma_logexp_time    <- c(PFOS = 0.25,    PFOA = 0.30,     PFHxS = 0.35,     PFNA = 0.40)
burden_loading       <- c(PFOS = 1.00,    PFOA = 0.80,     PFHxS = 0.60,     PFNA = 0.50)
sigma_burden         <- 0.35            # SD van de gedeelde blootstellingsburden

cat(sprintf("Ware ICC (random intercept / totale variantie) = %.3f\n", icc_true))

## ---- 2. Subject-niveau covariaten (tijdsonafhankelijk) ---------------------

subject_id <- seq_len(n_subjects)
gender     <- rbinom(n_subjects, 1, prob = 0.51)          # 1 = jongen, 0 = meisje
SES        <- factor(
  sample(SES_levels, n_subjects, replace = TRUE, prob = SES_probs),
  levels = SES_levels
)

## Subject-specifieke random intercepts (zorgen voor sterke correlatie tussen
## de 3 herhaalde metingen van dezelfde persoon, en tussen de 4 PFAS-stoffen)
b_outcome  <- rnorm(n_subjects, mean = 0, sd = sigma_b)
b_exposure <- rnorm(n_subjects, mean = 0, sd = sigma_burden)

subjects <- data.frame(
  id          = subject_id,
  gender      = gender,
  SES         = SES,
  b_outcome   = b_outcome,
  b_exposure  = b_exposure
)

## ---- 3. Long-format dataset: 1 rij per subject per meetmoment --------------

dat <- subjects[rep(seq_len(n_subjects), each = n_times), ]
dat$age <- rep(ages, times = n_subjects)
rownames(dat) <- NULL

## ---- 4. Simuleer de ware PFAS-blootstellingen (lognormaal, per stof) -------
## Voor elke stof in exposure_names: log-lineair model met gender- en
## leeftijdseffect plus de gedeelde blootstellingsburden (voor correlatie
## tussen de stoffen), gevolgd door LOD-censurering met de vaste, vooraf
## gekozen LOD (zie stap 1). Voor elke stof <naam> ontstaan drie kolommen:
##   <naam>      : ruwe meting, met -3 als sentinelwaarde voor metingen onder
##                 de LOD (in plaats van de werkelijke, onbekende waarde)
##   <naam>_lod  : constante kolom met de gehanteerde LOD voor die stof
##   <naam>_imp  : STOCHASTISCHE imputatie op basis van de (gefitte) lognormale
##                 verdeling: een links-gecensureerde normale MLE (mu, sigma)
##                 wordt geschat op log-schaal met enkel de boven-LOD-metingen
##                 en de kennis "waarde < LOD" voor de rest (analoog aan hoe
##                 een analist dit in de praktijk zou doen, zonder de ware
##                 waarden te kennen); de below-LOD-waarden worden vervolgens
##                 getrokken uit die gefitte verdeling, GETRUNKEERD op (0, LOD)
##                 -- dus een aselecte trekking i.p.v. de vaste LOD/sqrt(2)-
##                 substitutie van voorheen.
## below_LOD_<naam> blijft daarnaast beschikbaar als expliciete 0/1-indicator.

## Negatieve log-likelihood van een links-gecensureerde normale verdeling op
## log-schaal (par = c(mu, log(sigma)), sigma via log-parametrisatie zodat
## optim() niet buiten het toegelaten gebied (sigma > 0) kan optimaliseren).
.censored_lognormal_nll <- function(par, log_obs, below_flag, log_lod) {
  mu <- par[1]
  sigma <- exp(par[2])
  ll <- ifelse(
    below_flag == 0,
    dnorm(log_obs, mean = mu, sd = sigma, log = TRUE),
    pnorm(log_lod, mean = mu, sd = sigma, log.p = TRUE)
  )
  -sum(ll)
}

imputation_params <- list()   # bewaart de gefitte (mu, sigma) per stof, ter controle

for (expo in exposure_names) {
  log_mean <- mu_logexp_female[expo] +
    gender_logexp_effect[expo] * dat$gender +
    age_logexp_effect[expo]    * dat$age +
    burden_loading[expo]       * dat$b_exposure
  
  log_true <- log_mean + rnorm(nrow(dat), mean = 0, sd = sigma_logexp_time[expo])
  true_val <- exp(log_true)
  
  below_flag <- as.integer(true_val < LOD[expo])
  raw_val    <- ifelse(below_flag == 1, -3, true_val)              # sentinel = -3
  
  ## "Geobserveerde" log-waarde: gekend en gelijk aan log(true_val) boven de
  ## LOD, onbekend (NA) onder de LOD -- exact wat een analist in de praktijk
  ## zou hebben, gebruikt om (mu, sigma) van de onderliggende lognormale
  ## verdeling te schatten via censored-likelihood MLE.
  log_obs_for_fit <- ifelse(below_flag == 0, log_true, NA_real_)
  log_lod <- log(LOD[expo])
  
  start <- c(
    mean(log_obs_for_fit, na.rm = TRUE),
    log(sd(log_obs_for_fit, na.rm = TRUE))
  )
  fit <- optim(
    par = start, fn = .censored_lognormal_nll, method = "BFGS",
    log_obs = log_obs_for_fit, below_flag = below_flag, log_lod = log_lod
  )
  mu_hat    <- fit$par[1]
  sigma_hat <- exp(fit$par[2])
  imputation_params[[expo]] <- c(mu_hat = mu_hat, sigma_hat = sigma_hat)
  
  ## Stochastische imputatie: trek voor elke below-LOD waarneming een waarde
  ## uit N(mu_hat, sigma_hat), getrunkeerd op (-Inf, log_lod), via de inverse-
  ## CDF-methode; boven de LOD blijft de imputatie gelijk aan de (ware/
  ## geobserveerde) waarde.
  imp_val <- true_val
  n_below <- sum(below_flag == 1)
  if (n_below > 0) {
    u <- runif(n_below) * pnorm(log_lod, mean = mu_hat, sd = sigma_hat)
    log_imp_below <- qnorm(u, mean = mu_hat, sd = sigma_hat)
    imp_val[below_flag == 1] <- exp(log_imp_below)
  }
  
  dat[[paste0(expo, "_true")]]      <- true_val   # ware (ongekende) concentratie
  dat[[expo]]                       <- raw_val    # ruw, -3 = onder LOD
  dat[[paste0(expo, "_lod")]]       <- LOD[expo]   # constante LOD-kolom
  dat[[paste0(expo, "_imp")]]       <- imp_val     # geimputeerde versie
  dat[[paste0("below_LOD_", expo)]] <- below_flag
}

## ---- 5. Simuleer de uitkomst (neutrofielenaantal) via het mixed model ------
## Merk op: de vier stoffen worden hier ALLE VIER expliciet in het lineaire
## predictor opgenomen, met bewust verschillende sterkte -- PFOS heeft een
## sterk effect, PFOA een zwak effect, PFHxS en PFNA hebben coefficient 0 (dus
## per constructie geen effect, ook al staan ze wel in de formule). De
## uitkomst wordt gegenereerd op basis van de WARE (ongecensureerde)
## blootstelling, zoals in werkelijkheid het geval is -- de onderzoeker
## observeert enkel de (deels gecensureerde) <naam>/<naam>_imp-kolommen.

linpred <- beta0 +
  beta_age      * dat$age +
  beta_logPFOS  * log(dat$PFOS_true) +
  beta_logPFOA  * log(dat$PFOA_true) +
  beta_logPFHxS * log(dat$PFHxS_true) +
  beta_logPFNA  * log(dat$PFNA_true) +
  beta_gender   * dat$gender +
  beta_SES_midden * (dat$SES == "midden") +
  beta_SES_hoog   * (dat$SES == "hoog") +
  dat$b_outcome

dat$neutrophil_true <- linpred + rnorm(nrow(dat), mean = 0, sd = sigma_eps)
dat$neutrophil_true <- pmax(dat$neutrophil_true, 0.1)   # neutrofielenaantal kan niet negatief zijn

## ---- 6. Missing data -- ENKEL op de uitkomst, met VASTE steekproefgroottes -
## In plaats van een kans-gebaseerd MAR/MCAR-mechanisme wordt de uitval nu
## vastgelegd op exacte aantallen per golf: op jaar 7 blijven nog maar 220 van
## de 301 deelnemers over, op jaar 14 nog maar 200. De uitval is monotoon --
## wie al ontbreekt op jaar 7 blijft ook ontbreken op jaar 14 -- en wie het
## eerst afhaakt wordt bepaald door een dropout-risicoscore (lager SES-niveau
## + ruis -> hoger risico), zodat het MAR-verband met SES behouden blijft ook
## al liggen de steekproefgroottes per golf nu vast. Op jaar 0 (baseline) is
## er geen missingness: alle 301 deelnemers hebben er een uitkomstmeting.

n_wave2 <- 220   # aantal deelnemers met een uitkomstmeting op jaar 7
n_wave3 <- 200   # aantal deelnemers met een uitkomstmeting op jaar 14
stopifnot(n_wave2 <= n_subjects, n_wave3 <= n_wave2)

dropout_risk <- -0.50 * (subjects$SES == "midden") +
  -1.00 * (subjects$SES == "hoog") +
  rnorm(n_subjects)   # lager SES-niveau + ruis -> hoger dropoutrisico

## Rangschik subjecten van hoog naar laag risico: wie het hoogst scoort haakt
## het eerst af.
dropout_order <- order(dropout_risk, decreasing = TRUE)

n_drop_wave2       <- n_subjects - n_wave2   # al weg tegen jaar 7
n_drop_wave3_extra <- n_wave2 - n_wave3      # extra weg tussen jaar 7 en jaar 14

id_dropped_wave2 <- subject_id[dropout_order[seq_len(n_drop_wave2)]]
id_dropped_extra <- subject_id[dropout_order[seq(n_drop_wave2 + 1L, n_drop_wave2 + n_drop_wave3_extra)]]
id_dropped_wave3 <- c(id_dropped_wave2, id_dropped_extra)   # monotoon: wave2-uitval telt mee

dat$missing_outcome <- 0L
dat$missing_outcome[dat$age == 7  & dat$id %in% id_dropped_wave2] <- 1L
dat$missing_outcome[dat$age == 14 & dat$id %in% id_dropped_wave3] <- 1L

dat$neutrophil <- dat$neutrophil_true
dat$neutrophil[dat$missing_outcome == 1] <- NA

## ---- 7. Opkuisen en dataset klaarzetten ------------------------------------

dat$gender_f <- factor(dat$gender, levels = c(0, 1), labels = c("meisje", "jongen"))
dat$time     <- factor(dat$age, levels = ages, labels = c("0 jaar", "7 jaar", "14 jaar"))

expo_raw_cols   <- exposure_names
expo_lod_cols   <- paste0(exposure_names, "_lod")
expo_imp_cols   <- paste0(exposure_names, "_imp")
expo_below_cols <- paste0("below_LOD_", exposure_names)

dataset <- dat[, c("id", "time", "age", "gender_f", "SES",
                   expo_raw_cols, expo_lod_cols, expo_imp_cols, expo_below_cols,
                   "neutrophil")]
names(dataset)[names(dataset) == "gender_f"] <- "gender"

cat(sprintf("\nAantal rijen (subjecten x meetmomenten): %d\n", nrow(dataset)))
cat("SES-niveaus (aantal subjecten):\n")
print(table(subjects$SES))
cat("\n% metingen onder LOD, per stof (-3 = sentinel in de ruwe kolom):\n")
for (expo in exposure_names) {
  cat(sprintf("  %-6s (LOD = %.3f) : %.1f%%\n",
              expo, LOD[expo], 100 * mean(dat[[paste0("below_LOD_", expo)]])))
}
cat("\nGefitte censored-lognormal parameters (log-schaal), gebruikt voor de\n")
cat("stochastische imputatie van waarden onder de LOD:\n")
for (expo in exposure_names) {
  p <- imputation_params[[expo]]
  cat(sprintf("  %-6s mu_hat = %.3f, sigma_hat = %.3f\n", expo, p["mu_hat"], p["sigma_hat"]))
}
cat(sprintf("%% ontbrekende uitkomstmetingen             : %.1f%%\n",
            100 * mean(is.na(dataset$neutrophil))))
cat("Aantal deelnemers met een uitkomstmeting, per meetmoment (doel: 301 / 220 / 200):\n")
print(tapply(!is.na(dataset$neutrophil), dataset$time, sum))
cat("Ontbrekende uitkomst per meetmoment (%):\n")
print(round(100 * tapply(is.na(dataset$neutrophil), dataset$time, mean), 1))

## ---- 8. (optioneel) valideer het simulatie-DGP met lme4 --------------------
## Fit het model op de gesimuleerde data (complete-case op de uitkomst) om na
## te gaan of de geschatte parameters in de buurt van de ware waarden liggen.
## Vereist het lme4-package; wordt overgeslagen als dat niet beschikbaar is.

if (requireNamespace("lme4", quietly = TRUE)) {
  fit <- lme4::lmer(
    neutrophil ~ age + log(PFOS_imp) + log(PFOA_imp) + gender + SES + (1 | id),
    data = dataset
  )
  cat("\n--- lme4::lmer fit op de gesimuleerde (complete-case) data ---\n")
  print(summary(fit))
  
  vc      <- as.data.frame(lme4::VarCorr(fit))
  icc_est <- vc$vcov[vc$grp == "id"] / sum(vc$vcov)
  cat(sprintf("\nGeschatte ICC: %.3f (ware ICC was %.3f)\n", icc_est, icc_true))
} else {
  cat("\n(lme4 niet geinstalleerd -- installeer met install.packages('lme4') ",
      "om het mixed model te fitten en de simulatie te valideren)\n", sep = "")
}

## ---- 9. Wegschrijven --------------------------------------------------------

write.csv(dataset, "3xg_gesimuleerde_dataset.csv", row.names = FALSE)
cat("\nWeggeschreven naar: 3xg_gesimuleerde_dataset.csv\n")
