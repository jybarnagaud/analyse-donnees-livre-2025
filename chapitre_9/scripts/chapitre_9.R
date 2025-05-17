##----------------------------------------------------------------------------##

# JY. Barnagaud et O. Gimenez. Analyse de données pour les écologues. Biotope
# Editions, 2025.

#                       -------Chapitre 9-------                               #

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

# install.packages("pacman") # si vous n'avez pas encore le package {pacman}

pacman::p_load(
  # génériques
  tidyverse,
  reshape2,
  jtools,
  patchwork,
  cowplot,
  viridis,
  mapview,
  
  # modèles
  boot,
  polynom,
  mgcv,
  investr,
  
  # sélection de modèles
  lmtest,
  MuMIn,
  
  # représentations graphiques de modèles
  effects,
  sjPlot,
  ggeffects,
  visreg
  
)

# pour installer manuellement chaque package : install.packages(nomdupackage)
# pour charger manuellement chaque package : library(nomdupackage)

## options ---------------------------------------------------------------------

# enlève fond gris, garde grille et grossit la taille des légendes

theme_set(theme_light(base_size = 16))

# reproductibilité des simulations

set.seed(2020)

# choisir le chapitre (si vous passez par le projet analyse-donnees-livre-2025.Rproj)
# si vous passez par le projet chapitre_X.Rproj, ignorez cette commande

setwd("chapitre_9")

## jeux de données -------------------------------------------------------------

lomolino <- read.csv2("donnees/Lomolino.csv", dec = ".")
brg <- read.csv2("donnees/phenologie_arbres.csv")
clim <- read.csv2("donnees/climat_LR.csv", dec = ".")
mulette <- read.csv2("donnees/mulette_perliere.csv", dec = ".")
tetras <- read.csv2("donnees/tetras_lyre.csv", dec = ".")
amphi <- read.csv2("donnees/amphi_tendances.csv", dec = ".")
calais <- read.csv2("donnees/vegetation_pas_de_calais.csv", dec = ".")

## 9.1 Sélectionner un modèle parmi plusieurs-----------------------------------

# recharger les données (cf script du chapitre 8)

lomolino$log.surface <- log(lomolino$surface)
summary(lomolino)

## deux hypothèses emboitées:

# hypothèse spatiale

m.spatial <-
  glm(
    n_especes ~ log.surface + latitude + distance_source,
    family = "poisson",
    data = lomolino
  )

# hypothèse sur la surface

m.surface <-
  glm(n_especes ~ log.surface, family = "poisson", data = lomolino)

# on contrôle les résidus

par(mfrow = c(2, 2))
plot(m.spatial)

# l' autre modèle:

par(mfrow = c(2, 2))
plot(m.surface)

# ratio de vraisemblances (librairie lmtest)

lrtest(m.spatial, m.surface)

## comparer deux hypothèses non emboîtées

# modèle

m.ecol <-
  glm(
    n_especes ~ dens_foret + dens_prairie + dist_pleistocene + latitude,
    family = "poisson",
    data = lomolino
  )

# résidus

par(mfrow = c(2, 2))
plot(m.ecol)

# comparer les AIC

AIC(m.spatial, m.surface, m.ecol)

# exemple fictif AIC

L1 <- 10 # vraisemblance du modèle 1
L2 <- 10 # vraisemblance du modèle 2

2 * 3 - 2 * log(10)
2 * 4 - 2 * log(10)


## Sélectionner le modèle le plus parcimonieux

# Note : nous ne recommandons pas l'utilisation de cette méthode de réduction
# de modèle dans le cas général. Elle est éventuellement utile  pour certains usages
# spécifiques, mais que vous avez peu de chances de rencontrer. Nous la présentons
# uniquement parce que vous la rencontrerez dans la littérature, mais les
# recommandations des biostatisticiens vont actuellement plutôt dans le sens de
# ne construire qu'un seul modèle ou un jeu restreint de modèles tous biologiquement
# pertinents.

# données climatiques

summary(clim)

# on prépare les données

clim2 <-
  aggregate(clim[, "temperature_moyenne"], FUN = "mean", by = list(clim$POSTE, clim$annee))
colnames(clim2) <- c('POSTE', "annee", "temperature_moyenne")

clim3 <-
  aggregate(clim[, "cumul_precip"], FUN = "sum", by = list(clim$POSTE, clim$annee))
colnames(clim3) <- c('POSTE', "annee", "cumul_precip")

clim4 <- merge(clim2, clim3, by = c("POSTE", "annee"))
clim4$POSTE <- factor(clim4$POSTE)

summary(clim4)

# représentation graphique

ggplot(clim4) +
  aes(x = annee, y = temperature_moyenne, col = POSTE) +
  geom_line() +
  geom_point() +
  labs(x = "Années", y = "Température moyenne, °C") +
  theme_classic() +
  scale_colour_viridis_d()

# modèle maximum

mod.clim.max <-
  lm(
    temperature_moyenne ~ annee + POSTE + cumul_precip,
    data = clim4,
    na.action = na.fail
  )

# résidus

par(mfrow = c(2, 2))
plot(mod.clim.max)

# tous les modèles possibles avec AIC

sel.clim <- dredge(mod.clim.max, rank = "AIC")
sel.clim

# tous les modèles possibles avec  AICc

sel.clim.aicc <- dredge(mod.clim.max, rank = "AICc")
sel.clim.aicc

## Que faire si deux modèles sont équivalents ?

# choisir quelques modèles selon un critère donné

keep.clim <- get.models(sel.clim, delta <= 2)
keep.clim

keep.clim.alt = get.models(sel.clim, cumsum(weight) < 0.95)
keep.clim.alt

# moyenner les modèles

avg.clim <- model.avg(keep.clim)
summary(avg.clim)

# intervalle de confiance

ci.cond <- confint(avg.clim)
ci.cond

# autre version

ci.full <- confint(avg.clim, full = T)
ci.full

# représentation graphique :

par(mfrow = c(1, 2))
plot(avg.clim)
plot(avg.clim, intercept = F)

# FIGURE 9.4 - autre représentation possible :

p1 <- plot_model(avg.clim,
                 type = "est",
                 title = "",
                 colors = "viridis") +
  labs(x = "Coefficient de pente", y = "Ville") + font_size(axis_title.x = 12, axis_title.y = 12)

p2 <- plot_model(
  avg.clim,
  terms = "annee",
  type = "emm",
  title = "",
  colors = "viridis"
) +
  theme_classic() +
  labs(x = "Années", y = "Température moyenne annuelle prédite")

p3 <- plot_model(
  avg.clim,
  terms = "POSTE",
  type = "emm",
  title = "",
  colors = "viridis"
) +
  theme_classic() +
  labs(x = "Villes", y = "Température moyenne annuelle prédite")

p4 <- plot_model(
  avg.clim,
  terms = "cumul_precip",
  type = "emm",
  title = "",
  colors = "viridis"
) +
  theme_classic() +
  labs(x = "Cumul de précipiations", y = "Température moyenne annuelle prédite")

(p1 + p2) / (p3 + p4) +
  theme_sjplot(base_size = 12, base_family = "") +
  plot_annotation(tag_levels = "A")

## 9.2  Modéliser des interactions entre variables------------------------------

## Une interaction variable qualitative - variable quantitative continue

# les données : bourgeonnement pour 2 espèces (scores de 1 à 9)

brg$code_arbre <- factor(brg$code_arbre)
brg$essence <- factor(brg$essence)
brg$parcelle <- factor(brg$parcelle)

summary(brg)

# on recale la date

brg$date_rel_recal <- brg$date_rel - min(brg$date_rel) + 1

summary(brg)
dim(brg)

# représentation graphique des données brutes

brg %>%
  mutate(essence = fct_recode(essence, "Chêne" = "CHS", "Pin" = "PS")) %>%
  ggplot() +
  aes(x = date_rel_recal, y = bmean, color = essence) +
  geom_point(size = 3, alpha = 0.5) +
  scale_color_manual(values = c("#5ec962", "#440154")) +
  labs(x = "Date relative", y = "Degré de débourrement", color = NULL) +
  theme_classic()

# anova sur les espèces

bg.anova <- lm(bmean ~ essence, data = brg)
summary(bg.anova)

# résidus

par(mfrow = c(2, 2))
plot(bg.anova)

# histogramme de la variable de réponse

p1 <- ggplot(brg) +
  aes(x = bmean) +
  geom_histogram(bins = 9,
                 color = 'black',
                 fill = 'white') +
  labs(x = "Degré de débourrement", y = "Fréquence") +
  theme_classic()

p2 <- ggplot(brg) +
  aes(x = bmean) +
  geom_histogram(bins = 9,
                 color = 'black',
                 fill = 'white') +
  labs(x = "Degré de débourrement", y = "Fréquence") +
  theme_classic() +
  facet_wrap( ~ essence)

(p1 + plot_spacer()) / p2 +
  plot_annotation(tag_levels = "A")

# modèle à 2 variables

bg.mod <- lm(bmean ~ date_rel_recal + essence, data = brg)
summary(bg.mod)

# modèle à 2 variables + interaction

bg.mod.int <-
  lm(bmean ~ date_rel_recal + essence + (date_rel_recal:essence),
     data = brg)

bg.mod.int <- lm(bmean ~ date_rel_recal * essence, data = brg) # les deux écritures sont équivalentes

# résidus

par(mfrow = c(2, 2))
plot(bg.mod.int)

# on choisit la structure (étape facultative, en général si vous avez introduit
# une interaction dans votre modèle, vous savez pourquoi sur le plan biologique
# et n'avez donc pas de raison de la supprimer)

AIC(bg.mod, bg.mod.int)

# le summary

summary(bg.mod.int)

# représentation graphique, à la main

date.pred <- seq(
  from = min(brg$date_rel_recal),
  to = max(brg$date_rel_recal),
  length.out = 100
)

xpred.ch <- list(date_rel_recal = date.pred, essence = rep("CHS", 100))

ypred.ch <- predict(bg.mod.int, newdata = xpred.ch, interval = "confidence")


xpred.ps <- list(date_rel_recal = date.pred, essence = rep("PS", 100))

ypred.ps <- predict(bg.mod.int, newdata = xpred.ps, interval = "confidence")

df <- data.frame(
  date.pred = rep(date.pred, 2),
  essence = c(rep('Chêne', length(date.pred)), rep('Pin', length(date.pred))),
  xpred = c(xpred.ch$date_rel_recal, xpred.ps$date_rel_recal),
  ypred = rbind(ypred.ch, ypred.ps)
)

brg %>%
  mutate(essence = fct_recode(essence, "Chêne" = "CHS", "Pin" = "PS")) %>%
  ggplot() +
  aes(x = date_rel_recal, y = bmean, color = essence) +
  geom_point(size = 3, alpha = 0.5) +
  geom_line(data = df, aes(x = date.pred, y = ypred.fit, color = essence)) +
  geom_line(data = df,
            aes(x = date.pred, y = ypred.lwr, color = essence),
            lty = 2) +
  geom_line(data = df,
            aes(x = date.pred, y = ypred.upr, color = essence),
            lty = 2) +
  scale_color_manual(values = alpha(c("#5ec962", "#440154"), 0.5)) +
  ylim(0, 9) +
  expand_limits(x = 0, y = 0) +
  labs(x = "Date relative", y = "Degré de débourrement", color = NULL) +
  theme_classic()

# avec ggeffect (comme dans le livre)

plot(ggeffect(bg.mod.int, terms = c("date_rel_recal", "essence")), show_residuals = T) +
  labs(x = "Date relative", y = "Débourrement", col = "Essence") +
  theme_classic() +
  scale_colour_manual(
    aesthetics = c("fill", "colour"),
    values = c("#5ec962", "#440154"),
    labels = c("Chêne", "Pin")
  )

# avec sjPlot

plot_model(bg.mod.int, type = "int")

## Une interaction entre variables quantitatives

# charger les données

mulette$station <- factor(mulette$station)

head(mulette)
summary(mulette)

# relations avec les variables

p1 <- ggplot(mulette) +
  aes(x = factor(presence),
      y = pente,
      fill = factor(presence)) +
  geom_boxplot(alpha = 0.5) +
  labs(x = "Présence de la mulette", y = "Pente (%)") +
  theme_classic() +
  scale_fill_manual(values = c("#fde725", "#3b528b")) +
  theme(legend.position = "none")

p2 <- ggplot(mulette) +
  aes(x = factor(presence),
      y = coniferes,
      fill = factor(presence)) +
  geom_boxplot(alpha = 0.5) +
  labs(x = "Présence de la mulette", y = "Conifères (%)") +
  theme_classic() +
  scale_fill_manual(values = c("#fde725", "#3b528b")) +
  theme(legend.position = "none")

p3 <- ggplot(mulette) +
  aes(x = coniferes,
      y = pente,
      col = factor(presence)) +
  geom_point(size = 2, alpha = 0.5) +
  labs(y = "Pente(%)", x = "Conifères (%)", col = "Mulette") +
  theme_classic() +
  scale_color_manual(values = c("#fde725", "#3b528b"))

(p1 + p2) / (p3 + plot_spacer()) + plot_annotation(tag_levels = "A")

# modèle

mod.mulette <-
  glm(presence ~ pente * coniferes, family = binomial, data = mulette)

# résidus

par(mfrow = c(2, 2))
plot(mod.mulette)

# résumé du modèle

summary(mod.mulette)

# FIGURE 9.14 : représentation des effets

p1 <-
  plot(ggpredict(mod.mulette, terms = c("coniferes", "pente[0,1,100]")), show_residuals = TRUE) +
  labs(x = "Conifères (%)", y = "Probabilité de présence", title = "(A) valeurs spécifiques") +
  scale_y_continuous(labels = scales::percent, lim = c(0, 1)) +
  scale_color_viridis_d(aesthetics = c("color", "fill"))

p2 <-
  plot(ggpredict(mod.mulette, terms = c("coniferes", "pente[meansd]")), show_residuals = TRUE) +
  labs(x = "Conifères (%)", y = "Probabilité de présence", title = "(B) moyenne ± écart-type") +
  scale_y_continuous(labels = scales::percent, lim = c(0, 1)) +
  scale_color_viridis_d(aesthetics = c("color", "fill"))

p3 <-
  plot(ggpredict(mod.mulette, terms = c("coniferes", "pente[quart2]")), show_residuals = TRUE) +
  labs(x = "Conifères (%)", y = "Probabilité de présence", title = "(C) quartiles") +
  scale_y_continuous(labels = scales::percent, lim = c(0, 1)) +
  scale_color_viridis_d(aesthetics = c("color", "fill"))

p4 <-
  plot(ggpredict(mod.mulette, terms = c("coniferes", "pente[minmax]")), show_residuals = TRUE) +
  labs(x = "Conifères (%)", y = "Probabilité de présence", title = "(D) minimum - maximum") +
  scale_y_continuous(labels = scales::percent, lim = c(0, 1)) +
  scale_color_viridis_d(aesthetics = c("color", "fill"))

(p1 + p2) / (p3 + p4)

## Une interaction entre deux variables qualitatives

# données

calais$secteur <- factor(calais$secteur)
calais$quadrat <- factor(calais$quadrat)
calais$periode <- factor(calais$periode)
calais$periode <- relevel(calais$periode, ref = "fin_90")

summary(calais)

# visualisation de l'interaction sur données brutes

ggplot(calais) +
  aes(x = periode, y = richesse) +
  geom_boxplot() +
  facet_wrap(vars(secteur)) +
  theme_classic()

# le modèle

flore.calais <- lm(hveg ~ periode + secteur + (periode:secteur), data = calais)
flore.calais <- lm(hveg ~ periode * secteur, data = calais) # les deux codes sont équivalents

# résidus

par(mfrow = c(2, 2))
plot(flore.calais)

# résumé

summary(flore.calais)

## 9.3 Modéliser un optimum de réponse------------------------------------------

# simulation d'une courbe quadratique

x <- seq(-5, 5, length = 200)
a <- -1
b <- 0
c <- 1
y <- a * x ^ 2 + b * x + c

df <- data.frame(y = y, x = x)

ggplot(df) +
  geom_line(aes(x = x, y = y), linewidth = 2) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = "succès reproducteur", x = "climat") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  geom_segment(
    aes(
      x = -5,
      xend = 5,
      y = -25,
      yend = -25
    ),
    linewidth = 1,
    arrow = arrow(length = unit(0.6, "cm"))
  ) +
  geom_segment(
    aes(
      x = -5,
      xend = -5,
      y = -25,
      yend = 5
    ),
    linewidth = 1,
    arrow = arrow(length = unit(0.6, "cm"))
  ) +
  coord_cartesian(clip = "off") +
  annotate(
    "text",
    label = "froid-humide",
    x = -2.5,
    y = -24,
    color = "#3b528b",
    size = 6
  ) +
  annotate(
    "text",
    label = "chaud-sec",
    x = 2.5,
    y = -24,
    color = "#3b528b",
    size = 6
  )

# données sur la reproduction des tétras lyre dans les Préalpes

tetras$prop <- tetras$Nichees / tetras$Poules
tetras$UN <- factor(tetras$UN)
tetras$RN <- factor(tetras$RN)

summary(tetras)

# cartographie

carte.tetras <-
  mapview(
    tetras,
    xcol = "long_UN",
    ycol = "lat_UN",
    zcol = "prop",
    crs = "epsg:4326",
    map.types = "OpenStreetMap"
  )

# représentation graphique des variations de la NAO

synth.nao <- cbind(unique(tetras$Annee),
                   tapply(tetras$NAOdjfm, INDEX = tetras$Annee, FUN = "mean"))
colnames(synth.nao) <- c("annee", "moy")
synth.nao <- as.data.frame(synth.nao)
synth.nao <- synth.nao[order(synth.nao$annee), ]

ggplot(synth.nao) +
  aes(x = annee, y = moy) +
  geom_line() +
  geom_point(size = 5, colour = "white") +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = "Années", y = "NA0 (décembre - mars)") +
  theme_classic()

## Représenter un optimum dans un modèle linéaire

# forme quadratique hypothétique (librairie boot pour la fonction inv.logit)

nao.sim <- rnorm(1000, mean(synth.nao[, 2]), sd(synth.nao[, 2]))
y <- 0.5 + 0.1 * nao.sim - 0.3 * (nao.sim ^ 2)
psi <- inv.logit(y)
z <- cbind(nao.sim, psi)
z <- z[order(z[, 1]), ]
z <- as.data.frame(z)

ggplot(z) +
  aes(x = nao.sim, y = psi) +
  geom_line(linewidth = 1.5) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "NAO décembre - mars", y = "Taux de reproduction") +
  theme_classic()

## FIGURE 9.23 : formes des courbes produites par un polynôme de second degré

# formes quadratiques simulées inspirées de https://en.wikipedia.org/wiki/Quadratic_equation

# B = 0
n <- 200
x <- seq(-2, 2, length = n)
a <- c(-2, -1, 0, 1, 2)
grid <- expand.grid(x = x, a = a)
b <- 0
c <- 0
y <- grid$a * grid$x ^ 2 + b * grid$x + c

df <- data.frame(y = y,
                 x = x,
                 shape = c(
                   rep("a = -2", n),
                   rep("a = -1", n),
                   rep("a = 0", n),
                   rep("a = 1", n),
                   rep("a = 2", n)
                 ))

p1 <- ggplot(df) +
  geom_line(aes(
    x = x,
    y = y,
    group = shape,
    color = shape
  ), linewidth = 1.5) +
  labs(color = NULL, y = expression (y = ax ^ 2 + bx + c)) +
  theme_classic() +
  labs(subtitle = "(A)") +
  scale_color_viridis_d()

# B = 2
b <- 1
y <- grid$a * grid$x ^ 2 + b * grid$x + c

df <- data.frame(y = y,
                 x = x,
                 shape = c(
                   rep("a = -2", n),
                   rep("a = -1", n),
                   rep("a = 0", n),
                   rep("a = 1", n),
                   rep("a = 2", n)
                 ))

p2 <- ggplot(df) +
  geom_line(aes(
    x = x,
    y = y,
    group = shape,
    color = shape
  ), linewidth = 1.5) +
  labs(color = NULL, y = expression (y = ax ^ 2 + bx + c)) +
  theme_classic() +
  labs(subtitle = "(B)") +
  scale_color_viridis_d()

# B = 1
b <- (-2)
y <- grid$a * grid$x ^ 2 + b * grid$x + c

df <- data.frame(y = y,
                 x = x,
                 shape = c(
                   rep("a = -2", n),
                   rep("a = -1", n),
                   rep("a = 0", n),
                   rep("a = 1", n),
                   rep("a = 2", n)
                 ))

p3 <- ggplot(df) +
  geom_line(aes(
    x = x,
    y = y,
    group = shape,
    color = shape
  ), linewidth = 1.5) +
  labs(color = NULL, y = expression (y = ax ^ 2 + bx + c)) +
  theme_classic() +
  labs(subtitle = "(C)") +
  scale_color_viridis_d()

p1 + p2 + p3

## exploration graphique des données

p1 <- ggplot(tetras) +
  aes(x = NAOdjfm, y = prop) +
  geom_point(size = 2, alpha = 0.5) +
  labs(x = "NAO décembre - mars", y = "Succès reproducteur ") +
  labs(subtitle = "(A)") +
  theme_classic()

p2 <- ggplot(tetras) +
  aes(x = RN, y = prop) +
  geom_boxplot() +
  labs(x = "Région", y = "Succès reproducteur ") +
  coord_flip() +
  labs(subtitle = "(B)") +
  theme_classic()

p1 + p2

## différents codes possibles pour le même modèle

tetras.mod.quad.v1 <- glm(
  cbind(Nichees, Poules - Nichees) ~ poly(NAOdjfm, 2) + RN,
  family = "binomial",
  data = tetras
)

# alternative (voir le texte du livre pour l'explication de raw = T)

tetras.mod.quad <- glm(
  cbind(Nichees, Poules - Nichees) ~ poly(NAOdjfm, 2, raw = T) + RN,
  family = "binomial",
  data = tetras
)

# alternative avec la fonction I, qui équivaut à la précédente

tetras.mod.quad <- glm(
  cbind(Nichees, Poules - Nichees) ~ NAOdjfm + I(NAOdjfm ^ 2) + RN,
  family = "binomial",
  data = tetras
)

# on vérifie les résidus

par(mfrow = c(2, 2))
plot(tetras.mod.quad)

# modèle linéaire pour comparaison

tetras.mod.lin <- glm(cbind(Nichees, Poules - Nichees) ~ NAOdjfm + RN,
                      family = "binomial",
                      data = tetras)

# le terme quadratique est-il justifié?

AIC(tetras.mod.lin, tetras.mod.quad)

# résumé

summary(tetras.mod.quad)

# déviance expliquée

(276.84 - 187.52) / 276.84

# effet marginal

visreg(tetras.mod.quad, xvar = "NAOdjfm", scale = "response")

# avec ggeffect

plot(ggeffect(tetras.mod.quad, terms = "NAOdjfm"), show_residuals = T) +
  labs(title = "", x = "NAO décembre-mars", y = "Succès reproducteur marginal") +
  theme_classic()

## Calculer l’optimum d’un modèle quadratique

# la fonction :

y <- function(x) {
  coef(tetras.mod.quad)[1] +
    coef(tetras.mod.quad)[2] * x +
    coef(tetras.mod.quad)[3] * x * x
}

# calcul de l'optimum

xmax <-
  optimize(y, interval = range(tetras$NAOdjfm), maximum = TRUE)
xmax

# autre manière de faire avec la librairie polynom

p0 <- polynomial(c(
  coef(tetras.mod.quad)[1],
  coef(tetras.mod.quad)[2],
  coef(tetras.mod.quad)[3]
))

p1 <- deriv(p0)
p2 <- deriv(p1) # dérivées 1 et 2

xm <- solve(p1) # max et min de p0

xmax <- xm[predict(p2, xm) < 0] # selectionner les maxima
xmax

# la variance d'un ratio peut etre obtenue par la delta-method, voir eq (20) sur https://www.stat.cmu.edu/~hseltman/files/ratio.pdf

muR <- coef(tetras.mod.quad)[2]
muS <- coef(tetras.mod.quad)[3]
varR <- vcov(tetras.mod.quad)[2, 2]
varS <- vcov(tetras.mod.quad)[3, 3]
covRS <- vcov(tetras.mod.quad)[2, 3]

var_max <-
  muR ^ 2 / muS ^ 2 * (varR / muR ^ 2 - 2 * covRS / (muR * muS) +
                         varS / muS ^ 2)
var_max

# plot marginal avec l'optimum

plot(ggeffect(tetras.mod.quad, terms = "NAOdjfm"), show_residuals = T) +
  labs(title = "", x = "NAO décembre-mars", y = "Succès reproducteur marginal") +
  geom_vline(
    xintercept = xmax$maximum,
    col = "#3b528b",
    linetype = "dashed"
  ) +
  theme_classic() +
  annotate(
    "rect",
    xmin = xmax$maximum - 1.96 * sqrt(var_max),
    xmax = xmax$maximum + 1.96 * sqrt(var_max),
    ymin = -Inf,
    ymax = Inf,
    fill = "#3b528b",
    alpha = 0.2
  )

## Une interaction entre un terme quadratique et une variable qualitative : variation de l’optimum climatique des tétras par massifs

tetras.mod.quad.int <-
  glm(
    cbind(Nichees, Poules - Nichees) ~ poly(NAOdjfm, 2, raw = TRUE) * RN,
    family = "binomial",
    data = tetras
  )

# résidus

par(mfrow = c(2, 2))
plot(tetras.mod.quad.int)

# effets marginaux

visreg(tetras.mod.quad.int, xvar = "NAOdjfm", by = "RN")

# avec ggeffect

plot(ggeffect(tetras.mod.quad.int, terms = c("NAOdjfm", "RN")),
     show_residuals = T,
     facet = T) +
  labs(title = "", x = "NAO décembre-mars", y = "Succès reproducteur marginal") +
  scale_color_viridis_d(aesthetics = c("color", "fill")) +
  theme_classic()

# attention néanmoins à ne pas abuser des interactions :

AIC(tetras.mod.quad.int, tetras.mod.quad)

## 9.4 Les modèles additifs généralisés-----------------------------------------

## La question et les données : modéliser la phénologie de la grenouille rousse

# données de l'article de Prodon et al, Global Change Biology, 2017

amphi$codesp <- factor(amphi$codesp)
summary(amphi)

# modèle linéaire sur la Grenouille rousse

rantem <- subset(amphi, codesp == "RANTEM")

# représentation graphique données brutes

ggplot(rantem) +
  aes(x = an, y = date_resc) +
  geom_point(size = 3, alpha = 0.5) +
  labs(title = "Grenouille rousse", x = "Années", y = "Date de sortie moyenne") +
  theme_classic()

# modèle

rantem.mod <- lm(date_resc ~ an, data = rantem)

# résidus

par(mfrow = c(2, 2))
plot(rantem.mod)

# modèle quadratique

rantem.mod.quad <- lm(date_resc ~ poly(an, 2, raw = TRUE), data = rantem)

# résidus

par(mfrow = c(2, 2))
plot(rantem.mod.quad)

# lequel est le meilleur? aucun.

AIC(rantem.mod, rantem.mod.quad)

# représentation graphique

visreg(rantem.mod.quad, xvar = "an")

# avec ggeffect

plot(ggeffect(rantem.mod.quad, terms = "an"), show_residuals = T) +
  labs(title = "Grenouille rousse", x = "Années", y = "Date de sortie") +
  theme_classic()

## La spline, faire passer une fonction au plus près des données

# GAM sur la grenouille rousse avec la librairie mgcv

rantem.mod.gam <- gam(date_resc ~ s(an), data = rantem)

# checks résidus

gam.check(rantem.mod.gam, old.style = T)

# résumé

summary(rantem.mod.gam)

# comparaison graphique entre le GAM et le modèle quadratique (en R de base, plus facile à exploiter avec la librairie mgcv)

plot(
  rantem.mod.gam,
  shift = coef(rantem.mod.gam)[1],
  se = T,
  ylim = range(rantem$date_resc),
  lwd = 2,
  rug = F,
  shade = T,
  xlab = "Années",
  ylab = "Date de sortie",
  main = ""
)
points(
  rantem$an,
  rantem$date_resc,
  pch = 21,
  bg = "gray80",
  col = "gray80",
  cex = 1
)

p.quad <- predict(rantem.mod.quad,
                  newdata = list(an = seq(min(rantem$an), max(rantem$an), by = 1)),
                  se.fit = TRUE)

lines(seq(min(rantem$an), max(rantem$an), by = 1),
      p.quad$fit,
      col = "#3b528b",
      lwd = 2)
low <- p.quad$fit - 1.96 * p.quad$se.fit
high <- p.quad$fit + 1.96 * p.quad$se.fit
lines(
  seq(min(rantem$an), max(rantem$an), by = 1),
  low,
  lwd = 2,
  lty = 3,
  col = "#3b528b"
)
lines(
  seq(min(rantem$an), max(rantem$an), by = 1),
  high,
  lwd = 2,
  lty = 3,
  col = "#3b528b"
)
legend(
  "topright",
  lwd = 2,
  col = c("gray80", "#3b528b"),
  legend = c("GAM", "modèle quadratique")
)

## Le compromis ajustement – parcimonie avec un GAM

# contraintes sur les k

rantem.mod.gam.k3 <- gam(date_resc ~ s(an, k = 3), data = rantem)
rantem.mod.gam.k4 <- gam(date_resc ~ s(an, k = 4), data = rantem)
rantem.mod.gam.k5 <- gam(date_resc ~ s(an, k = 5), data = rantem)
rantem.mod.gam.k6 <- gam(date_resc ~ s(an, k = 6), data = rantem)

# plot des contraintes sur k

par(mfrow = c(2, 2))
plot(
  rantem.mod.gam.k3,
  shift = coef(rantem.mod.gam)[1],
  show_residuals = T,
  cex = 5,
  rug = F,
  shade = T,
  ylab = "Date de sortie",
  xlab = "Années"
)
legend("topright", legend = paste(c("k max = ", "edf ="), c(3, round(
  summary(rantem.mod.gam.k3)$edf, 2
)), sep = " "), bty = "n")
plot(
  rantem.mod.gam.k4,
  shift = coef(rantem.mod.gam)[1],
  show_residuals = T,
  cex = 5,
  rug = F,
  shade = T,
  ylab = "Date de sortie",
  xlab = "Années"
)
legend("topright", legend = paste(c("k max = ", "edf ="), c(4, round(
  summary(rantem.mod.gam.k4)$edf, 2
)), sep = " "), bty = "n")
plot(
  rantem.mod.gam.k5,
  shift = coef(rantem.mod.gam)[1],
  show_residuals = T,
  cex = 5,
  rug = F,
  shade = T,
  ylab = "Date de sortie",
  xlab = "Années"
)
legend("topright", legend = paste(c("k max = ", "edf ="), c(5, round(
  summary(rantem.mod.gam.k5)$edf, 2
)), sep = " "), bty = "n")
plot(
  rantem.mod.gam.k6,
  shift = coef(rantem.mod.gam)[1],
  show_residuals = T,
  cex = 5,
  rug = F,
  shade = T,
  ylab = "Date de sortie",
  xlab = "Années"
)
legend("topright", legend = paste(c("k max = ", "edf ="), c(6, round(
  summary(rantem.mod.gam.k6)$edf, 2
)), sep = " "), bty = "n")

## La phénologie de plusieurs espèces : un GAM avec une interaction entre une spline et une variable qualitative

# GAM avec toutes les espèces

all.sp.gam <- gam(date_resc ~ codesp + s(an, by = codesp), data = amphi)

# courbes correspondantes

titres.plots <- data.frame(
  codes = levels(amphi$codesp),
  especes = c(
    "Alyte accoucheur",
    "Crapaud commun",
    "Rainette méridionale",
    "Triton palmé",
    "Pélodyte ponctué",
    "Grenouille rousse",
    "Triton marbré"
  )
)

par(mfrow = c(3, 3))
for (i in 1:nlevels(amphi$codesp)) {
  sub_data <- amphi[amphi$codesp == i, ]
  plot_title <- titres.plots[i, "especes"]
  
  plot.gam(all.sp.gam, select = i, main = plot_title)
  
}

## Encart 9.1. Utiliser un GAM pour l’interpolation : une carte de chaleur

# le gam comme outil spatial

head(tetras[, c("Nichees", "Poules", "long_UN", "lat_UN")])

# ici le gam  va nous servir à créer une surface de réponse

tetras.mod.gam.spatial <-
  gam(
    cbind(Nichees, Poules - Nichees) ~ s(long_UN, lat_UN, k = 16),
    family = "binomial",
    data = tetras
  )

# une carte de chaleur :

vis.gam(
  tetras.mod.gam.spatial,
  view = c("lat_UN", "long_UN"),
  plot.type = "contour"
)

# variations de la carte de chaleur selon la complexité du modèle

k.seq <- c(4, 8, 12, 16)

par(mfrow = c(2, 2))

for (i in 1:length(k.seq)) {
  tetras.mod.gam.spatial <-
    gam(
      cbind(Nichees, Poules - Nichees) ~ s(long_UN, lat_UN, k = k.seq[i]),
      family = "binomial",
      data = tetras
    )
  
  vis.gam(
    tetras.mod.gam.spatial,
    view = c("lat_UN", "long_UN"),
    plot.type = "contour",
    main = paste("k = ", k.seq[i], sep = "")
  )
}

## 9.5 Une très brève introduction aux modèles non linéaires--------------------

## on recharge les données de débourrement et on sépare les espèces

# données brutes

pp <- brg %>%
  mutate(essence = fct_recode(
    essence,
    'Pin sylvestre' = 'PS',
    'Chêne sessile' = 'CHS'
  )) %>%
  ggplot() +
  aes(x = as_factor(date_rel_recal),
      y = bmean,
      fill = essence) +
  geom_boxplot() +
  labs(x = "Date", y = "Degré de débourrement", fill = NULL) +
  scale_fill_manual(values = c("#5ec962", "#440154")) +
  theme_classic() +
  theme(legend.position = "top")
pp

# simulation d'une sigmoïde

sigmoid <- function(x, a, b, c, d) {
  a + (b - a) / (1 + exp((c - x) / d))
}

x <- rnorm(100, 0, 1)
a <- 1
b <- 1.2
c <- 0
d <- 0.3

y <- sigmoid(x, a, b, c, d)
df <- data.frame(x = x, y = y)

ggplot(df) +
  aes(x = x, y = y) +
  geom_line(size = 1.5) +
  geom_hline(
    yintercept = a,
    linetype = "dashed",
    col = "#3b528b",
    linewidth = 1.2
  ) +
  annotate(
    "text",
    x = -2,
    y = a + 0.01,
    label = "a",
    col = "#3b528b",
    size = 8
  ) +
  geom_hline(
    yintercept = b,
    linetype = "dashed",
    col = "#3b528b",
    linewidth = 1.2
  ) +
  annotate(
    "text",
    x = -2,
    y = b - 0.01,
    label = "b",
    col = "#3b528b",
    size = 8
  ) +
  geom_vline(
    xintercept = c,
    linetype = "dashed",
    col = "#21918c",
    linewidth = 1.2
  ) +
  annotate(
    "text",
    x = c - 0.2,
    y = 1.15,
    label = "c",
    col = "#21918c",
    linewidth = 8
  ) +
  geom_segment(
    x = -1,
    xend = 1,
    y = sigmoid(-1, a, b, c, d),
    yend = sigmoid(1, a, b, c, d),
    linetype = "dashed",
    col = "tan2",
    linewidth = 1.2
  ) +
  annotate(
    "text",
    x = 0.8,
    y = 1.15,
    label = "d",
    col = "tan2",
    linewidth = 8
  ) +
  theme_classic()

# d'abord les pins ! sigmoïde à 4 param

brg.pin <- subset(brg, essence == "PS")

# le modèle

deb.pin.mod <-
  nls(bmean ~ SSfpl(date_rel_recal, a, b, c, d), data = brg.pin)

# vérification des résidus et ajustement

df <- data.frame(pred = predict(deb.pin.mod), res = residuals(deb.pin.mod))

# plot

p1 <- ggplot(df) +
  aes(x = pred, y = res) +
  geom_point(size = 2, alpha = 0.5) +
  labs(x = "Valeurs prédites", y = "Résidus") +
  theme_classic()

df2 <- data.frame(bmean = brg.pin$bmean, pred = predict(deb.pin.mod))

p2 <- ggplot(df2) +
  aes(x = bmean, y = pred) +
  geom_point(size = 3, alpha = 0.5) +
  geom_abline(
    aes(intercept = 0, slope = 1),
    lty = 2,
    color = '#21918c',
    size = 1.5
  ) +
  labs(x = "Degrés de débourrement réels", y = "Degrés de débourrement prédits") +
  theme_classic()

p1 + p2 +
  plot_annotation(tag_levels = "A")

# résumé

summary(deb.pin.mod)

# intervalles de confiance

confint(deb.pin.mod)

# puis les chênes

brg.ch <- subset(brg, essence == "CHS")

# le modèle

deb.ch.mod <-
  nls(bmean ~ SSfpl(date_rel_recal, a, b, c, d), data = brg.ch)
summary(deb.ch.mod)
confint(deb.ch.mod)

# plot des deux modèles sur les données brutes (librairie investr pour l'intervalle de confiance)

prd.pin <-
  predFit(deb.pin.mod, interval = "prediction", level = 0.95)

z <- data.frame(
  prd.pin = prd.pin[, 'fit'],
  low = prd.pin[, 'lwr'],
  high = prd.pin[, 'upr'],
  date = brg.pin$date_rel_recal,
  fdate = as_factor(brg.pin$date_rel_recal)
)

z <- z[order(z[, 'fdate']), ]

prd.ch <- predFit(deb.ch.mod, interval = "prediction", level = 0.95)

z2 <- data.frame(
  prd.ch = prd.ch[, 'fit'],
  low = prd.ch[, 'lwr'],
  high = prd.ch[, 'upr'],
  date = brg.ch$date_rel_recal,
  fdate = as_factor(brg.ch$date_rel_recal)
)

z2 <- z2[order(z2[, 'fdate']), ]

pp +
  geom_line(
    data = z,
    aes(x = as.numeric(fdate), y = prd.pin),
    inherit.aes = FALSE,
    size = 1,
    color = "#440154"
  ) +
  geom_ribbon(
    data = z,
    aes(
      x = as.numeric(fdate),
      ymin = low,
      ymax = high
    ),
    inherit.aes = FALSE,
    alpha = 0.2,
    fill = "#440154"
  ) +
  geom_line(
    data = z2,
    aes(x = as.numeric(fdate), y = prd.ch),
    inherit.aes = FALSE,
    size = 1,
    color = "#5ec962"
  ) +
  geom_ribbon(
    data = z2,
    aes(
      x = as.numeric(fdate),
      ymin = low,
      ymax = high
    ),
    inherit.aes = FALSE,
    alpha = 0.2,
    fill = "#5ec962"
  )
