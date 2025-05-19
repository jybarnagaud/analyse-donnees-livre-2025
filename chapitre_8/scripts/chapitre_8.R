##----------------------------------------------------------------------------##

# JY. Barnagaud et O. Gimenez. Analyse de données pour les écologues. Biotope
# Editions, 2025.

#                       -------Chapitre 8-------                               #

## commandes R pour la reproduction des exemples du chapitre
## à jour 23 janvier 2025

## Il est recommandé d'utiliser ce script avec la version de R la plus récente,
# consulter https://cran.r-project.org/

## les packages requis en début de script doivent être installés avant de lancer
# les commandes. La disponibilité et les mises à jours de packages ne sont pas
# garanties.

## ouvrez ce script dans RStudio depuis le projet chapitre_2.Rproj pour plus de
# flexibilité. En cas d'utilisation hors du projet, spécifiez le répertoire
# courant à partir de la fonction setwd()

##----------------------------------------------------------------------------##

## packages --------------------------------------------------------------------

# utiliser cette fonction pour installer et charger automatiquement les packages

install.packages("pacman")

pacman::p_load(
  # génériques
  tidyverse,
  formatR,
  jtools,
  reshape2,
  
  # graphiques
  patchwork,
  viridis,
  ggplot2,
  cowplot,
  
  # modèles linéaires généralisés
  arm,
  # fonction invlogit
  MASS,
  # glm avec loi négative binomiale
  pscl,
  # zero inflated poisson
  glmmTMB,
  # autre codage de la zero inflated poisson
  
  # sorties graphiques de modèles
  effects,
  sjPlot,
  visreg,
  ggeffects,
  
  # interprétation et diagnostics de modèles
  DHARMa,
  biostat3,
  questionr,
  oddsratio,
  ggfortify
  
)

# pour installer manuellement chaque package : install.packages(nomdupackage)
# pour charger manuellement chaque package : library(nomdupackage)

## options ---------------------------------------------------------------------

# enlève fond gris, garde grille et grossit la taille des légendes

theme_set(theme_light(base_size = 16))

# reproductibilité des simulations

set.seed(150)

# choisir le chapitre (si vous passez par le projet analyse-donnees-livre-2025.Rproj)
# si vous passez par le projet chapitre_X.Rproj, ignorez cette commande

setwd("chapitre_8")

## jeux de données -------------------------------------------------------------

lomolino <- read.csv2("donnees/lomolino.csv", dec = ".", row.names = 1)
sterne <- read.csv("donnees/sterne_pierregarin_hyeres.csv")
tetras <- read.csv2("donnees/tetras_beaufortain.csv")
rifleman <- read.csv2("donnees/xenique.csv", dec = ".")

## 8.1 Modéliser des données de comptages---------------------------------------

## La question et les données : un test de la théorie des îles

# données

head(lomolino)
summary(lomolino)

# exploration graphique

p1 <- lomolino %>%
  ggplot() +
  aes(x = distance_source, y = n_especes) +
  geom_point(size = 3,
             alpha = 0.5,
             col = "#3b528b") +
  labs(x = "Distance à la source (km)", y = "Richesse spécifique") +
  theme_classic()

p2 <- lomolino %>%
  ggplot() +
  aes(x = latitude, y = n_especes) +
  geom_point(size = 3,
             alpha = 0.5,
             col = "#3b528b") +
  labs(x = "Latitude (°)", y = "Richesse spécifique") +
  theme_classic()

p3 <- lomolino %>%
  ggplot() +
  aes(x = surface, y = n_especes) +
  geom_point(size = 3,
             alpha = 0.5,
             col = "#3b528b") +
  labs(x = "Surface (km²)", y = "Richesse spécifique") +
  theme_classic()

p4 <- lomolino %>%
  ggplot() +
  aes(x = log(surface), y = n_especes) +
  geom_point(size = 3,
             alpha = 0.5,
             col = "#3b528b") +
  labs(x = "log(Surface (km²))", y = "Richesse spécifique") +
  theme_classic()

(p1 + p2) / (p3 + p4) + plot_annotation(tag_levels = 'A')

# régulariser la variable surface en la transformant en log :

lomolino$log.surface <- log(lomolino$surface)

## Les limites de la régression linéaire et des transformations

# histogramme de la variable de réponse

lomolino %>%
  ggplot() +
  aes(x = n_especes) +
  geom_histogram(bins = 9,
                 color = 'black',
                 fill = 'white') +
  labs(x = "Richesse spécifique", y = "Effectif") +
  theme_classic()

# transformation de la richesse spécifique en racine carrée

sq.n_especes <- sqrt(lomolino$n_especes)
lomolino$sq.n_especes <- sq.n_especes

# histogrammes

p1 <- lomolino %>%
  ggplot() +
  aes(x = sq.n_especes) +
  geom_histogram(bins = 8,
                 color = 'black',
                 fill = 'white') +
  labs(x = expression(sqrt("Richesse spécifique")), y = "Effectif") +
  theme_classic()

p2 <- lomolino %>%
  ggplot() +
  aes(x = n_especes, y = sq.n_especes) +
  geom_point(size = 2,
             alpha = 0.5,
             col = "#3b528b") +
  labs(x = "Richesse spécifique", y = expression(sqrt("Richesse spécifique"))) +
  theme_classic()

p1 + p2

# transformation log

lg.n_especes <- log(lomolino$n_especes + 1)
lomolino$log.n_especes <- lg.n_especes

# effet de la transformation log

p1 <- lomolino %>%
  ggplot() +
  aes(x = log.n_especes) +
  geom_histogram(bins = 8,
                 color = 'black',
                 fill = 'white') +
  labs(x = "log(Richesse spécifique + 1)", y = "Effectif") +
  theme_classic()

p2 <- lomolino %>%
  ggplot() +
  aes(x = n_especes, y = log.n_especes) +
  geom_point(size = 2,
             alpha = 0.5,
             col = "#3b528b") +
  labs(x = "Richesse spécifique", y = "log(Richesse spécifique + 1)") +
  theme_classic()

p1 + p2

# test de Shapiro sur la richesse log transformée

shapiro.test(lomolino$log.n_especes)
shapiro.test(lomolino$sq.n_especes)

# modèle avec NSPEC transformée en racine carrée

sq.lomo <- lm(sq.n_especes ~ distance_source + latitude + log.surface, data = lomolino)

# résidus

par(mfrow = c(2, 2))
plot(sq.lomo)

# résultats

summary(sq.lomo)

# représentation graphique avec la librairie ggeffects
# NB : plusieurs fonctions sont possibles avec cette librairie. Elles ont des comportements légèrement différents, mais dans la majorité des cas simples les rendus seront similaires. Consultez la documentation de la librairie afin de bien comprendre les différences.

p1 <- plot(
  ggeffect(sq.lomo, terms = "distance_source"),
  show_residuals = T,
  show_title = F,
  colors = "#3b528b",
  dot_alpha = 0.5
) +
  labs(x = "Distance à la source (km)", y = expression(sqrt("Richesse spécifique"))) +
  theme_classic()
p1

# prédiction

x.pred <- list(
  distance_source = 700,
  latitude = mean(lomolino$latitude),
  log.surface = mean(lomolino$log.surface)
)
predict(sq.lomo, newdata = x.pred)

# le même modèle avec une transformation log

log.lomo <-
  lm(log.n_especes ~ distance_source + latitude + log.surface, data = lomolino)

x.pred <- list(
  distance_source = 700,
  latitude = mean(lomolino$latitude),
  log.surface = mean(lomolino$log.surface)
)
predict(log.lomo, newdata = x.pred)

# le même modèle sans transformation:

lomo <-
  lm(n_especes ~ distance_source + latitude + log.surface, data = lomolino)
x.pred <- list(
  distance_source = 700,
  latitude = mean(lomolino$latitude),
  log.surface = mean(lomolino$log.surface)
)
predict(lomo, newdata = x.pred)

## Le GLM de Poisson (modèle linéaire généralisé de Poisson ou régression de Poisson)

# simuler une loi de Poisson vs une loi de Poisson surdispersée

df <- data.frame(
  poisson = rpois(1000, mean(lomolino$n_especes)),
  surdpoisson = rpois(1000, mean(lomolino$n_especes)) * rbinom(1000, 1, 0.6)
)

# histogrammes des simulations :

p1 <- df %>%
  ggplot() +
  aes(x = poisson) +
  geom_histogram(bins = 10,
                 color = 'black',
                 fill = 'white') +
  labs(x = "Loi de Poisson", y = "Effectif") +
  theme_classic()

p2 <- df %>%
  ggplot() +
  aes(x = surdpoisson) +
  geom_histogram(bins = 10,
                 color = 'black',
                 fill = 'white') +
  labs(x = "Loi de Poisson surdispersée", y = "Effectif") +
  theme_classic()

p1 + p2

# moyenne et variance des données de Lomolino

mean(lomolino$n_especes)
var(lomolino$n_especes)

# simulation effet du lien log

xtest <- rnorm(1000, 0, 1)
lxtest <- exp(xtest)

df <- data.frame(xtest = xtest, lxtest = lxtest)

p1 <- df %>%
  ggplot() +
  aes(x = xtest, y = lxtest) +
  geom_point(size = 2, alpha = 0.5) +
  labs(x = "x", y = "exp(x)") +
  theme_classic()

p2 <- df %>%
  ggplot() +
  aes(x = lxtest, y = log(lxtest)) +
  geom_point(size = 2, alpha = 0.5) +
  labs(x = "exp(x)", y = "log(x)") +
  theme_classic()

p1 + p2

# Le GLM de Poisson sous R

lomo.pois <- glm(
  n_especes ~ distance_source + latitude + log.surface,
  family = "poisson",
  data = lomolino
)

# Contrôler les conditions d'application

par(mfrow = c(2, 2))
plot(lomo.pois)

# résultats

summary(lomo.pois)

# Mesurer l'ajustement d'un GLM de Poisson : pourcentage de déviance expliquée

100 * (lomo.pois$null.deviance - lomo.pois$deviance) / lomo.pois$null.deviance

# Vérifier la surdispersion résiduelle (méthode approximative, pratique pour sa rapidité)

res.lomo.pois <- residuals(lomo.pois, type = "pearson")
N <- nrow(lomolino)
k <- length(coef(lomo.pois))
sum(res.lomo.pois ^ 2) / (N - k)

# tests de surdispersion avec la librairie DHARma (à préférer à la précédente pour des résultats définitifs)

testDispersion(lomo.pois)

# Interpréter les effets des variables dans un GLM de Poisson

summary(lomo.pois)

## Les incidence-rate ratios : faciliter l’interprétation biologique des paramètres d’un GLM de Poisson

# incidence rate ratios avec la librairie biostat3

eform(lomo.pois)

# avec la librairie oddsratio, on peut choisir l'intervalle sur lequel calculer l'incidence rate ratio. Ici, on veut l'IRR pour un incrément de distance à la source de 1, 50, 100km, et pour l'étendue totale de la variable.

or_glm(lomolino,
       lomo.pois,
       incr = list(
         distance_source = 1,
         latitude = 0,
         log.surface = 0
       ))

or_glm(lomolino,
       lomo.pois,
       incr = list(
         distance_source = 50,
         latitude = 0,
         log.surface = 0
       ))

or_glm(lomolino,
       lomo.pois,
       incr = list(
         distance_source = 100,
         latitude = 0,
         log.surface = 0
       ))

range.distance <-
  max(lomolino$distance_source)  - min(lomolino$distance_source)

or_glm(
  lomolino,
  lomo.pois,
  incr = list(
    distance_source = range.distance,
    latitude = 0,
    log.surface = 0
  )
)

## Représentation graphique des effets d’un GLM de Poisson

# avec ggeffects, sur l'échelle naturelle

p1 <- plot(
  ggeffect(lomo.pois, terms = "distance_source"),
  show_residuals = T,
  show_title = F
) +
  labs(x = "Distance à la source (km)", y = "Richesse spécifique marginale") +
  theme_classic()

p2 <- plot(
  ggeffect(lomo.pois, terms = "log.surface"),
  show_residuals = T,
  show_title = F
) +
  xlab("log(Surface (km²))") +
  ylab("Richesse spécifique marginale") +
  theme_classic()

p1 + p2

# à la main, sur l'échelle du lien

preds <- predict.glm(
  object = lomo.pois,
  newdata = data.frame(
    distance_source = lomolino$distance_source,
    latitude = mean(lomolino$latitude),
    log.surface = mean(lomolino$log.surface)
  ),
  type = 'link',
  se.fit = TRUE
)

upper <- preds$fit + 1.96 * preds$se.fit
lower <- preds$fit - 1.96 * preds$se.fit
predframe <- data.frame(
  lwr = lower,
  upr = upper,
  fit = preds$fit,
  distance_source = lomolino$distance_source,
  n_especes = lomolino$n_especes
)

p3 <- ggplot(predframe, aes(distance_source, log(n_especes))) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = 'grey90') +
  geom_point(color = 'grey60') +
  geom_line(aes(distance_source, fit),
            col = 'black',
            linewidth = .7) +
  labs(x = "Distance à la source (km)", y = "Richesse spécifique (échelle log)") +
  theme_classic()

preds <- predict.glm(
  object = lomo.pois,
  newdata = data.frame(
    distance_source = mean(lomolino$distance_source),
    latitude = mean(lomolino$latitude),
    log.surface = lomolino$log.surface
  ),
  type = 'link',
  se.fit = TRUE
)

upper <- preds$fit + 1.96 * preds$se.fit
lower <- preds$fit - 1.96 * preds$se.fit

predframe <- data.frame(
  lwr = lower,
  upr = upper,
  fit = preds$fit,
  log.surface = lomolino$log.surface,
  n_especes = lomolino$n_especes
)

p4 <- ggplot(predframe, aes(log.surface, log(n_especes))) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = 'grey90') +
  geom_point(color = 'grey60') +
  geom_line(aes(log.surface, fit),
            col = 'black',
            linewidth = .7) +
  labs(x = "log(Surface (km²))", y = "Richesse spécifique (échelle log)") +
  theme_classic()

# Figure 8.12:

(p1 + p2) / (p3 + p4) + plot_annotation(tag_levels = 'A')

# nb : pour choisir une autre valeur de référence que la moyenne des covariables, par exemple ici en les annulant : utiliser ggpredict - voir ci-dessous

# attention, ce n'est que pour exemple : cette prédiction n'est pas pertinente parce que nous n'avons ni une latitude de 0, ni une log.surf=0 dans nos données

p3 <- plot(
  ggpredict(
    lomo.pois,
    terms = "distance_source",
    condition = c(latitude = 0, log.surface = 0)
  ),
  show_residuals = T,
  show_title = F
) +
  labs(x = "Distance à la source (km)", y = "Richesse spécifique marginale") +
  theme_classic()

# effet distance aux minimum et maximum de latitude

p4 <- plot(
  ggpredict(
    lomo.pois,
    terms = "distance_source",
    condition = c(latitude = min(lomolino$latitude))
  ),
  show_residuals = T,
  show_title = F
) +
  labs(x = "Distance à la source (km)", y = "Richesse spécifique marginale", title = "marge sud, log(surface moyenne)") +
  theme_classic()

p5 <- plot(
  ggpredict(lomo.pois, terms = "distance_source", condition = c(LATI = max(lomolino$latitude))),
  show_residuals = T,
  show_title = F
) +
  labs(x = "Distance à la source (km)", y = "Richesse spécifique marginale", title = "marge nord, log(surface moyenne)") +
  theme_classic()

p1b <- p1 + labs(title = "latitude moyenne, log(surface moyenne)")

p1b + p4 + p5

## 8.2 Modéliser des données binaires-------------------------------------------

## La question et les données : variations temporelles d’occupation d’îlots chez la sterne pierregarin

# les données

head(sterne)
summary(sterne)

# comptages par années : plan d'échantillonnage parfaitement équilibré

table(sterne[, c("ilot", "annee")])
tapply(sterne$ilot, INDEX = sterne$annee, FUN = "length")

# représentation graphique

p0 <- ggplot(sterne) +
  aes(x = annee, y = presence) +
  geom_bar(stat = "identity") +
  labs(x = "Années", y = "Nombre d'îlots occupés") +
  theme_classic()

# exprimé en ratio - par ans

synthese.annees <-
  aggregate(
    sterne$presence,
    by = list(sterne$annee),
    FUN = function(x) {
      100 * sum(x) / 16
    }
  )

synthese.annees <- rbind(synthese.annees, c(2019, NA))
colnames(synthese.annees) = c("annee", "proportion")

head(synthese.annees)

# exprimé en ratio - par colonies

synthese.colonies <-
  aggregate(sterne$presence, by = list(sterne$ilot), FUN = sum)
colnames(synthese.colonies) = c("ilot", "duree")

head(synthese.annees)

# représentations graphiques des synthèses

p1 <- sterne %>%
  mutate(ilot = as.factor(ilot)) %>%
  group_by(ilot) %>%
  mutate(synthese.colonies = sum(presence)) %>%
  ungroup() %>%
  mutate(ilot = fct_reorder(ilot, synthese.colonies)) %>%
  ggplot(aes(x = factor(annee, levels = c(2012:2021)), y = ilot)) +
  geom_tile(aes(fill = as.factor(presence)), color = "gray70") +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual(values = c("white", "black")) +
  labs(x = "Années", y = "Ilots", fill = "Occupation") +
  theme_classic()

p2 <- ggplot(synthese.annees) +
  aes(x = annee, y = proportion) +
  geom_point(size = 2) +
  geom_line(size = 1) +
  labs(x = "Années", y = "Pourcentage d'îlots occupés") +
  theme_classic()

(p0 + plot_spacer()) / (p1 + p2) + plot_annotation(tag_levels = 'A')

## Le modèle linéaire généralisé binomial avec des données binaires

mod.sternes0 <- lm(proportion ~ annee, data = synthese.annees)

# résidus

par(mfrow = c(2, 2))
plot(mod.sternes0)

# résultats

summary(mod.sternes0)

# prédictions : interpolation. Pourquoi pas, mais on agrège les colonies souvent occupées et celles qui ne le sont pas

predict(mod.sternes0, newdata = list(annee = 2019))

# extrapolation --> proportion d'occupation inférieure à 0

predict(mod.sternes0, newdata = list(annee = 1990))

# extrapolation --> proportion d'occupation au-delà de 100%

predict(mod.sternes0, newdata = list(annee = 2052))

# Le lien logit

df.psi <- data.frame(logit.psi <- rnorm(100, 0, 3), psi <- invlogit(logit.psi))

ggplot(df.psi) +
  aes(x = psi, y = logit.psi) +
  geom_line(size = 1) +
  labs(x = "psi", y = "logit(psi)") +
  theme_classic()

# le GLM binomial avec une loi de Bernoulli

sterne$ilot <- factor(sterne$ilot)

sternes.bern <- glm(presence ~ annee + ilot, family = "binomial", data =
                      sterne)

# résidus

par(mfrow = c(2, 2))
plot(sternes.bern)

# résultats

summary(sternes.bern)

# odds ratios avec la librairie questionr

odds.ratio(sternes.bern)

# odds ratios avec la librairie biostat3

eform(sternes.bern)

# odds ratios avec la librairie oddsratio (incréments d'un an, trois ans, cinq ans)

or_glm(sterne, sternes.bern, incr = list(annee = 1))
or_glm(sterne, sternes.bern, incr = list(annee = 3))
or_glm(sterne, sternes.bern, incr = list(annee = 5))

# prédictions

predict(sternes.bern, newdata = list(annee = 2019, ilot = "A"))

predict(
  sternes.bern,
  newdata = list(annee = 2019, ilot = "A"),
  type = "response",
  se.fit = T
)

predict(
  sternes.bern,
  newdata = list(annee = 2019, ilot = "MR_IF1"),
  type = "response",
  se.fit = T
)

# graphiques avec ggeffects

p3 <- plot(ggpredict(sternes.bern, terms = "annee"), show_residuals = T) +
  labs(x = "Années", y = "Probabilité de colonisation", title = "") +
  theme_classic()

p4 <- plot(ggpredict(sternes.bern, terms = "ilot"), show_residuals = T) +
  labs(x = "Ilots", y = "Probabilité de colonisation", title = "") +
  theme(axis.text.x = element_text(angle = 90)) +
  theme_classic() +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  ))

p3 + p4

## 8.3 Modéliser des proportions------------------------------------------------

## La question et les données : le succès reproducteur du tétras lyre

# données sur les tétras

tetras$UN <- factor(tetras$UN)
head(tetras)
summary(tetras)

# données brutes

ggplot(tetras) +
  aes(x = UN, y = 100 * Nichees / Poules) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  labs(x = "Sites", y = "Taux de reproduction") +
  theme_classic()

# modèle

repro.tetras <- glm(cbind(Nichees, Poules - Nichees) ~ UN,
                    data = tetras,
                    family = binomial)

# résidus

par(mfrow = c(2, 2))
plot(repro.tetras)

# résultats

summary(repro.tetras)

# surdispersion (la surdispersion n'existe pas avec une loi de Bernoulli mais doit être contrôlée avec une loi binomiale)

res.repro.tetras <- residuals(repro.tetras, type = "pearson")
N <- nrow(tetras)
k <- length(coef(repro.tetras))
sum(res.repro.tetras ^ 2) / (N - k)

# vérification avec la librairie DHARMa

testDispersion(repro.tetras)

# odds ratios

odds.ratio(repro.tetras)

# représentation graphique des résultats

probs <- ggpredict(repro.tetras, terms = "UN")

plot(probs, show_residuals = T) +
  labs(x = "Sites", y = "Taux de reproduction", title = "") +
  theme_classic()

## 8.4 Les données de comptage surdispersées------------------------------------

## La question et les données : en forêt avec le Titipounamu

# les données

rifleman$region <- factor(rifleman$region)
head(rifleman)
summary(rifleman)

## Détecter un excès de zéros

#explorer les données

p1 <- rifleman %>%
  ggplot() +
  aes(x = PROPFOR, y = total) +
  geom_point(size = 2, alpha = 0.5) +
  labs(x = "% de forêt dans le paysage", y = "Comptages de Xéniques") +
  theme_classic()

p2 <- rifleman %>%
  ggplot() +
  aes(x = PROPNATFOR, y = total) +
  geom_point(size = 2, alpha = 0.5) +
  labs(x = "% de forêt native par parcelle", y = "") +
  theme_classic()

p3 <- rifleman %>%
  ggplot() +
  aes(x = SHA500, y = total) +
  geom_point(size = 2, alpha = 0.5) +
  labs(x = "Hétérogénéité du paysage", y = "") +
  theme_classic()

p4 <- rifleman %>%
  ggplot() +
  aes(x = total) +
  geom_histogram(bins = 10,
                 color = 'black',
                 fill = 'white') +
  labs(x = "nombre de Xéniques", y = "Effectif") +
  theme_classic()

(p1 + p2) / (p3 + p4) + plot_annotation(tag_levels = "A")

# surdispersion liée aux zéros

mean(rifleman$total)
var(rifleman$total)
sum(rifleman$total == 0) / nrow(rifleman)  # proportion de zéros dans le jeu de données

# données simulées

rif.sim <- rpois(1000, mean(rifleman$total))

rif.sim %>%
  as.data.frame() %>%
  ggplot() +
  aes(x = rif.sim) +
  geom_histogram(bins = 10,
                 color = 'black',
                 fill = 'white') +
  labs(x = "Données simulées", y = "Effectif") +
  theme_classic()

mean(rif.sim)
var(rif.sim)
sum(rif.sim == 0) / length(rif.sim)

## Commençons simple : un GLM de Poisson

# GLM de Poisson

rif.pois <- glm(total ~ PROPFOR + PROPNATFOR + SHA500,
                family = poisson,
                data = rifleman)

# résidus

par(mfrow = c(2, 2))
plot(rif.pois)

# surdispersion

res.rif <- residuals(rif.pois, type = "pearson")
N  <- nrow(rifleman)
p  <- length(coef(rif.pois))
sum(res.rif ^ 2) / (N - p)

# avec DHARMa :

testDispersion(rif.pois)

## Une première solution : le modèle quasi-Poisson

# modèle quasi-Poisson

rif.qpois <- glm(total ~ PROPFOR + PROPNATFOR + SHA500,
                 family = quasipoisson,
                 data = rifleman)
# résidus

par(mfrow = c(2, 2))
plot(rif.qpois)

# résultats

summary(rif.qpois)

# surdispersion

res.qrif <- residuals(rif.qpois, type = "pearson")
N  <- nrow(rifleman)
p  <- length(coef(rif.qpois))
sum(res.qrif ^ 2) / (N - p)

## Le GLM négatif-binomial

# données simulées (librairie MASS pour la fonction rnegbin)

lambda.sim.pois <- rpois(1000, mean(rifleman$total))
lambda.sim.nb <-
  rnegbin(1000, mu = mean(rifleman$total), theta = 0.3)

p1 <-
  lambda.sim.pois %>%
  as.data.frame() %>%
  ggplot() +
  aes(x = lambda.sim.pois) +
  geom_histogram(bins = 10,
                 color = 'black',
                 fill = 'white') +
  labs(x = "Loi de Poisson, espérance = 0.17", y = "Effectif")  +
  theme_classic()

breaks <-
  unique(unlist(ggplot_build(p1)$data[[1]][, c("xmin", 'xmax')]))

p2 <-
  lambda.sim.nb %>%
  as.data.frame() %>%
  ggplot() +
  aes(x = lambda.sim.nb) +
  geom_histogram(
    bins = 10,
    color = 'black',
    fill = 'white',
    breaks = breaks
  ) +
  labs(x = "Loi Négative Binomiale, espérance = 0.17, dispersion = 0.3", y = "Effectif") +
  theme_classic()

p1 / p2 + plot_annotation(tag_levels = "A")

## modèle négatif binomial (librairie MASS)

rif.negbin <- glm.nb(total ~ PROPFOR + PROPNATFOR + SHA500, data = rifleman)

# résidus

par(mfrow = c(2, 2))
plot(rif.negbin)

# résultats

summary(rif.negbin)

# surdispersion

res.rif.nb <- residuals(rif.negbin, type = "pearson")
N  <- nrow(rifleman)
p  <- length(coef(rif.negbin))
sum(res.rif.nb ^ 2) / (N - p)

# avec DHARMa:

testDispersion(rif.negbin)

# résultats :

summary(rif.negbin)

# odds ratios :

eform(rif.negbin)

# effets marginaux :

p1 <- plot(ggeffect(rif.negbin, terms = "PROPFOR"), show_residuals = T) +
  theme_classic()

p2 <- plot(ggeffect(rif.negbin, terms = "PROPNATFOR"), show_residuals = T) +
  theme_classic()

p3 <- plot(ggeffect(rif.negbin, terms = "SHA500"), show_residuals = T) +
  theme_classic()

p4 <- ggplot(rif.negbin, aes(x = .fitted, y = .resid)) +
  geom_point(size = 2, alpha = 0.2) +
  labs(x = "Valeurs prédites", y = "Résidus") +
  theme_classic()

(p1 + p2) / (p3 + p4)

## Le modèle ZIP (Zero Inflated Poisson)

# variation régionale des occurrences, des fréquences et des effectifs lorsque l'espèce est présente:

# nombre de points d'occurrence

rif.occ <-
  aggregate(
    rifleman$total,
    by = list(rifleman$region),
    FUN = function(x) {
      sum(x > 0)
    }
  )
colnames(rif.occ) = c("region", "tot.occ")

# fréquence relative

n.region <- table(rifleman$region)
rif.occ$freq <- rif.occ$tot.occ / n.region

# graphiques

p0 <- ggplot(rif.occ) +
  aes(x = fct_reorder(region, -tot.occ), y = tot.occ) +
  geom_col(fill = "steelblue", col = "white") +
  theme_classic() +
  labs(x = "Région", y = "Nombre de points d'occurrence")

p1 <- ggplot(rif.occ) +
  aes(x = fct_reorder(region, -tot.occ), y = freq) +
  geom_col(fill = "steelblue", col = "white") +
  theme_classic() +
  labs(x = "Région", y = "Fréquence d'occurrence relative")

rif.tot.pres <- droplevels(subset(rifleman, total > 0))
rif.tot.pres$region <-
  factor(rif.tot.pres$region, levels = c("Mont", "Banks", "Foot_hills"))

p2 <- ggplot(rif.tot.pres) +
  aes(x = region, y = total) +
  geom_boxplot() +
  theme_classic() +
  labs(x = "Région", y = "Comptage de xéniques par point de présence")

(p0 + p1) / (p2 + plot_spacer())

# moyenne et variance des comptages sur les points de présence :

mean(rif.tot.pres$total)
var(rif.tot.pres$total)

# plots données brutes (les mêmes qu'en début de partie)

ggplot(rifleman) +
  aes(x = total) +
  geom_histogram(bins = 10,
                 color = 'black',
                 fill = "white") +
  labs(x = "Comptage de xéniques", y = "Effectif") +
  theme_classic()

# modèle avec la librairie pscl

rif.zip <- zeroinfl(total ~ PROPFOR + PROPNATFOR + SHA500 | region,
                    dist = "poisson",
                    data = rifleman)

# résidus

df <- data.frame(pred = predict(rif.zip),
                 res = residuals(rif.zip, type = "pearson"))
ggplot(df) +
  aes(x = pred, y = res) +
  geom_point(size = 2, alpha = 0.5) +
  labs(x = "Valeurs prédites", y = "Résidus") +
  theme_classic()

#  surdispersion (approximée)

res.rif.zip <- residuals(rif.zip, type = "pearson")
N  <- nrow(rifleman)
p  <- length(coef(rif.zip))
sum(res.rif.zip ^ 2) / (N - p)

# FIGURE 8.36: surdispersion avec DHARMa 
# (nécessite de recoder la ZIP avec la librairie glmmTMB, voir : https://github.com/florianhartig/DHARMa/issues/16)

rif.zip.2 <-
  glmmTMB(
    total ~ PROPFOR + PROPNATFOR + SHA500,
    ziformula = ~ region ,
    family = "poisson",
    data = rifleman
  )
testDispersion(rif.zip.2)

# résultats

summary(rif.zip)

# odds ratios :

eform(rif.zip)

# simulation pour expliciter le comportement de la fonction zeroinfl : la partie "Zero-inflation" du summary exprime les paramètres comme une probabilité d'avoir un zéro (et non une probabilité d'avoir un 1, comme c'est le cas dans un GLM binomial habituel) (cette simulation ne figure pas dans l'ouvrage)

A <- rbinom(60, 1, 0.9)
B <- rbinom(60, 1, 0.5)
C <- rbinom(60, 1, 0.2)

df = data.frame(bloc = rep(c("A", "B", "C"), each = 60), pres = c(A, B, C))

df$count = 0

n.count <- length(which(df$pres == 1))
rand.count <- rpois(n.count, 10)

df[which(df$pres == 1), "count"] = rand.count

# GLM binomial

glb <- glm(pres ~ bloc, family = "binomial", data = df)
summary(glb)

# ZIP (pas de covariable sur la couche des comptages)

zip <- zeroinfl(count ~ 1 | bloc, dist = "poisson", data = df)
summary(zip)

# On voit bien que les paramètres sont exactement opposés dans le GLM binomial (proba de présence) et dans la ZIP (proba d'absence)

## Encart 8.1 : le modèle Hurdle, une variante du modèle ZIP--------------------

rif.hur <- hurdle(total ~ PROPFOR + PROPNATFOR + SHA500 | region,
                  dist = "poisson",
                  data = rifleman)

# surdispersion résiduelle approximée

res.rif.hur <- residuals(rif.hur, type = "pearson")
N  <- nrow(rifleman)
p  <- length(coef(rif.hur))
sum(res.rif.hur ^ 2) / (N - p)

# résidus

df <- data.frame(pred = predict(rif.hur),
                 res = residuals(rif.hur, type = "pearson"))
df %>%
  ggplot() +
  aes(x = pred, y = res) +
  geom_point(size = 3, alpha = 0.5) +
  labs(x = "Valeurs prédites", y = "Valeurs résiduelles")

# résultats

summary(rif.hur)

# comparaison des AIC

AIC(rif.zip, rif.hur)

## Traiter les excès de zéros : qu’est-ce que ça change vraiment ?

# Table des AIC

AIC(rif.pois, rif.negbin, rif.zip, rif.hur)

# intervalles de confiance des différents modèles

conf.pois <- confint(rif.pois)
conf.qpois <- confint(rif.qpois)
conf.negbin <- confint(rif.negbin)
conf.zip <- confint(rif.zip)
conf.hur <- confint(rif.hur)

# graphique

df <- as.data.frame(rbind(conf.pois, conf.qpois, conf.negbin, conf.zip, conf.hur))

colnames(df) <- c('low2.5', 'up2.5')

df$est <- c(coef(rif.pois),
            coef(rif.qpois),
            coef(rif.negbin),
            coef(rif.zip),
            coef(rif.hur))

df$param <- c(
  row.names(conf.pois),
  row.names(conf.qpois),
  row.names(conf.negbin),
  row.names(conf.zip),
  row.names(conf.hur)
)

row.names(df) <- NULL

df$model <- rep(
  c(
    "GLM Poisson",
    "GLM quasi-Poisson",
    "GLM nég. binom.",
    "ZIP",
    "Hurdle"
  ),
  c(
    nrow(conf.pois),
    nrow(conf.qpois),
    nrow(conf.negbin),
    nrow(conf.zip),
    nrow(conf.hur)
  )
)

df %>%
  mutate(
    param = as_factor(param),
    param = fct_recode(
      param,
      Foothills = "zero_regionFoot_hills",
      Mont = "zero_regionMont",
      Plain = "zero_regionPlain",
      PROPFOR = "count_PROPFOR",
      PROPNATFOR = "count_PROPNATFOR",
      SHA500 = "count_SHA500"
    )
  ) %>%
  filter(param %in% c('PROPFOR', 'PROPNATFOR', 'SHA500', 'Foothills', 'Mont', 'Plain')) %>%
  ggplot() +
  aes(x = param, y = est, color = model) +
  geom_pointrange(aes(ymax = up2.5, ymin = low2.5), position = position_dodge(width = 0.4)) +
  scale_x_discrete("") +
  coord_cartesian(ylim = c(-4, 4)) + # contraint axe y sans supprimer des données comme avec ylim(-4,4)
  geom_hline(yintercept = 0, lty = 2) +
  labs(color = NULL, x = "Variable", y = "Paramètre - IC 95%") +
  scale_color_viridis_d() +
  theme_classic() +
  geom_vline(xintercept = 3.5) +
  annotate(
    "text",
    x = c(2, 5),
    y = -4,
    label = c("couche poissonienne", "couche binomiale")
  )