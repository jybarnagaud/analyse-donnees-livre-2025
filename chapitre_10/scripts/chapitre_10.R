##----------------------------------------------------------------------------##

# JY. Barnagaud et O. Gimenez. Analyse de données pour les écologues. Biotope
# Editions, 2025.

#                       -------Chapitre 10-------                               #

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

# install.packages("pacman")

pacman::p_load(
  # génériques
  reshape2,
  tidyr,
  ggplot2,
  patchwork,
  viridis,
  sf,
  mapview,
  leaflet,
  lwgeom,
  rnaturalearth,
  rnaturalearthdata,
  formatR,
  
  # modèles mixtes
  lme4,
  
  # représentations graphiques
  PerformanceAnalytics,
  ggeffects,
  sjPlot,
  
  # diagnostics et post-analyses
  DHARMa,
  questionr,
  broom.mixed,
  # odds ratios
  MuMIn,
  performance,
  see
)

# pour installer manuellement chaque package : install.packages(nomdupackage)
# pour charger manuellement chaque package : library(nomdupackage)

## options ---------------------------------------------------------------------

# enlève fond gris, garde grille et grossit la taille des légendes

theme_set(theme_light(base_size = 16))

# reproductibilité des simulations

set.seed(2020)

# options pour les cartes

mapviewOptions(fgb = FALSE)

## jeux de données -------------------------------------------------------------

aqui <-
  read.csv2("donnees/chenille_processionnaire.csv", encoding = "UTF-8")

tetras <- read.csv2("donnees/tetras_lyre.csv", dec = ".")
couesnon <- read.csv2("donnees/avifaune_couesnon.csv")
geojson <- st_read("donnees/mailles_landbio.geojson")

## 10.1 Le problème des strates-------------------------------------------------

# La question et les données : tendance temporelle d’infestation de la chenille processionnaire

aqui$placette <- factor(aqui$placette)
aqui$region <- factor(aqui$region)
aqui$prop_attaq <- aqui$nbpinsattac / aqui$nbpins
summary(aqui)

# Carte des placettes (données brutes)

regions <- st_read("donnees/regions.geojson")

ggplot(regions) +
  geom_sf() +
  geom_point(
    data = aqui,
    mapping = aes(x = longitude, y = latitude, colour = region),
    size = 2
  ) +
  scale_color_viridis(discrete = TRUE) +
  labs(color = "Régions") +
  theme_classic()

# Séries temporelles d'infestation par régions :

ggplot(aqui) +
  aes(
    x = factor(annee),
    y = prop_attaq,
    group = placette,
    color = region
  ) +
  geom_line(show.legend = FALSE) +
  labs(x = "Années", y = "Proportion de pins infestés") +
  facet_wrap(~ region, ncol = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_color_viridis_d()

# On change l'échelle des années par commodité d'interprétation :

aqui$annee_resc <- aqui$annee - min(aqui$annee) + 1

## Un premier modèle

# Le GLM :

mod0 <- glm(
  cbind(nbpinsattac, nbpins - nbpinsattac) ~ annee_resc,
  family = binomial,
  data = aqui
)

# Interprétation :

summary(mod0)

# Odds ratios :

odds.ratio(mod0)

# Résidus :

par(mfrow = c(2, 2))
plot(mod0)

# cartographie des résidus:

aqui$res.mod0 <- residuals(mod0)

p0 <- plot(ggpredict(mod0, terms = c("annee_resc")), show_residuals = T)

p1a <- ggplot(regions) +
  geom_sf() +
  geom_point(
    data = aqui,
    mapping = aes(x = longitude, y = latitude, colour = res.mod0),
    size = 2
  ) +
  scale_color_viridis(discrete = F) +
  labs(color = "Résidus du \n GLM binomial") +
  theme_classic()

p1b <- ggplot(aqui) +
  aes(x = placette, y = res.mod0) +
  geom_boxplot() +
  labs(x = "Placettes", y = "Résidus du GLM binomial") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))

p1a + p1b

## 10.2 Représenter une stratification par un modèle mixte----------------------

## Cadre d’application et limites de l’approche

# Carte : illustration de la notion d'effets aléatoires. Les placettes échantillonnées sont représentées comme des triangles rouges. Les placettes non échantillonnées sont les croix blanches.

rpoints <- st_sample(regions, 700) %>%
  st_sf() %>%
  st_transform(4326)

ggplot(regions) +
  geom_sf() +
  geom_sf(data = subset(regions, jointure_g %in% c("A", "B", "F", "G", "J")), aes(fill = jointure_g)) +
  geom_point(
    data = aqui,
    mapping = aes(x = longitude, y = latitude),
    colour = "red",
    pch = 17,
    size = 2
  ) +
  geom_sf(
    data = rpoints,
    pch = 3,
    col = 'white',
    alpha = 0.67
  ) +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  labs(color = "Régions") +
  theme_classic() +
  theme(legend.position = "none")

## L’implémentation sous R

# recodage des niveaux de facteurs par commodité d'affichage (en particulier pour la figure 17), sans impact sur le modèle.

aqui$placette2 <- NA

for (i in 1:nlevels(aqui$region)) {
  aqui[which(aqui$region == levels(aqui$region)[i]), "placette2"] <-
    paste(substr(levels(aqui$region)[i], 1, 1), aqui[which(aqui$region == levels(aqui$region)[i]), "placette"], sep =
            "")
}

aqui$placette2 <- factor(aqui$placette2)

# Un codage possible du modèle, avec la librairie lme4. On peut aussi utiliser la librairie glmmPQL (ou la librairie mgcv, qui réexploite cette dernière). En revanche, la librairie nlme, principale alternative à lme4 pour les modèles mixtes, n'ajuste pas de GLM.

mod2 <- glmer(
  cbind(nbpinsattac, nbpins - nbpinsattac) ~ annee_resc + (1 |
                                                             region / placette2),
  family = binomial,
  data = aqui
)

# Diagnostic basique des résidus :

plot(mod2, type = c("p", "smooth"), col.line = 1)

# QQ plot :

plot(mod2, xlab = "Valeurs prédites", ylab = "Résidus de Pearson du GLMM")
qqnorm(scale(residuals(mod2)))
abline(0, 1)

# hétéroscédasticité :

plot(
  mod2,
  sqrt(abs(resid(.))) ~ fitted(.),
  type = c("p", "smooth"),
  col.line = 1,
  xlab = "Valeurs prédites",
  ylab = "sqrt(abs(résidus))"
)

# surdispersion :

simulationOutput <- simulateResiduals(fittedModel = mod2)
testDispersion(simulationOutput)

## Interpréter un modèle mixte

# Le résumé :

summary(mod2)

# Les odds-ratios, avec un intervalle de confiance de Wald :

tidy(
  mod2,
  conf.int = TRUE,
  exponentiate = TRUE,
  effects = "fixed",
  conf.method = "Wald"
)

# Les odds-ratios, avec un intervalle de confiance par profilage :

tidy(
  mod2,
  conf.int = TRUE,
  exponentiate = TRUE,
  effects = "fixed",
  conf.method = "profile"
)

# Effet prédit des années (l'argument jitter = 0 évite un sur-étalement artefactuel des points prévu par ggpredict pour éviter les empilements de données sur le graphique) :

plot(ggpredict(mod2, terms = "annee_resc"),
     residuals = T,
     jitter = 0) + labs(x = "Années", y = "Taux d'infestation", title = "")

# répartition des placettes par région

res.aqui <- unique(aqui[, c("placette", "region")])
tapply(res.aqui$placette, INDEX = res.aqui$region, FUN = "length")

# Effet de l'incertitude sur les effets aléatoires :

p1 <- ggpredict(mod2, type = "re")
plot(p1, colors = "bw", show_title = F)

# Impact des effets aléatoires sur la tendance temporelle:

plot(ggpredict(
  mod2,
  terms = c("annee_resc", "region"),
  type = "re"
),
show_residuals = T,
jitter = F) +
  scale_color_viridis_d() +
  theme_classic() +
  labs(x = "Année relative", y = "Probabilité d'infestation", title = "")

# Effets aléatoires et leur incertitude - il s'agit d'une prédiction a posteriori :

p1 <-
  plot_model(mod2,
             type = "re",
             show.values = TRUE,
             title = "Effet aléatoire Placettes")[[1]]
p2 <-
  plot_model(mod2,
             type = "re",
             show.values = TRUE,
             title = "Effet aléatoire Régions")[[2]]

p1 + p2 + plot_annotation(tag_levels = "A")

# Estimation du R² :

r.squaredGLMM(mod2)

## pentes aléatoires  : exemple des tétras

# Variabilité d’une courbe entre strates : variations inter-massifs de la reproduction du tétras lyre dans les Alpes

tetras$prop <- tetras$Nichees / tetras$Poules
tetras$UN <- factor(tetras$UN)
tetras$RN <- factor(tetras$RN)

# Le modèle du chapitre 9

tetras.mod.quad.int  <- glm(
  cbind(Nichees, Poules - Nichees) ~ poly(NAOdjfm, 2, raw = TRUE) *
    RN,
  family = "binomial",
  data = tetras
)

# Modèle avec effets aléatoires emboités sur l'ordonnée à l'origine : même type de modèle qu'avec les chenilles

glmm.tetras  <- glmer(
  cbind(Nichees, Poules - Nichees) ~ poly(NAOdjfm, 2, raw = TRUE) + (1 |
                                                                       RN / UN),
  family = "binomial",
  data = tetras
)

# le modèle à pentes aléatoires

glmm.tetras.slope  <-
  glmer(
    cbind(Nichees, Poules - Nichees) ~ poly(NAOdjfm, 2, raw = TRUE) + (poly(NAOdjfm, 2, raw = TRUE) |
                                                                         RN / UN),
    family = "binomial",
    data = tetras
  )

# résidus

plot(glmm.tetras.slope,
     type = c("p", "smooth"),
     col.line = 1)

# graphique suivant

plot(
  glmm.tetras.slope,
  sqrt(abs(resid(.))) ~ fitted(.),
  type = c("p", "smooth"),
  col.line = 1,
  xlab = "Valeurs prédites",
  ylab = "sqrt(abs(résidus))"
)

# graphique suivant

plot(glmm.tetras.slope, xlab = "Valeurs prédites", ylab = "Résidus de Pearson du GLMM")
qqnorm(scale(residuals(glmm.tetras.slope)))
abline(0, 1)

# surdispersion

simulationOutput <-
  simulateResiduals(fittedModel = glmm.tetras.slope)
testDispersion(simulationOutput)

# pour comparaison : sans interaction ni effet aléatoire

glm.tetras <-
  glm(
    cbind(Nichees, Poules - Nichees) ~ poly(NAOdjfm, 2, raw = TRUE) +
      RN,
    family = "binomial",
    data = tetras
  )

# pour comparaison : sans aucun effet hiérarchique

glm.null <-
  glm(
    cbind(Nichees, Poules - Nichees) ~ poly(NAOdjfm, 2, raw = TRUE),
    family = "binomial",
    data = tetras
  )

# interprétation

summary(glmm.tetras)

# avec pente aléatoire

summary(glmm.tetras.slope)

# représentation comparative

p7b <-
  plot(ggpredict(glm.null, terms = "NAOdjfm [all]"), show_residuals = T) +
  labs(x = "NAO décembre - mars", y = "Probabilité de reproduction", title = "sans effet hiérarchique") +
  scale_y_continuous(labels = scales::percent, lim = c(0, 1)) +
  scale_x_continuous(lim = c(-1.5, 2))

p7 <-
  plot(ggpredict(glm.tetras, terms = "NAOdjfm [all]"), show_residuals = T) +
  labs(x = "NAO décembre - mars", y = "Probabilité de reproduction", title = "effet fixe RN") +
  scale_y_continuous(labels = scales::percent, lim = c(0, 1)) +
  scale_x_continuous(lim = c(-1.5, 2))

p8 <-
  plot(ggpredict(glmm.tetras, terms = "NAOdjfm [all]"),
       show_residuals = T) +
  labs(x = "NAO décembre - mars", y = "Probabilité de reproduction", title = "effets aléatoires RN/UN \n sur l'ordonnée à l'origine") +
  scale_y_continuous(labels = scales::percent, lim = c(0, 1)) +
  scale_x_continuous(lim = c(-1.5, 2))

p9 <-
  plot(ggpredict(glmm.tetras.slope, terms = "NAOdjfm [all]"),
       show_residuals = T) +
  labs(x = "NAO décembre - mars", y = "Probabilité de reproduction", title = "effets aléatoires RN/UN \n sur l'ordonnée à l'origine et la pente") +
  scale_y_continuous(labels = scales::percent, lim = c(0, 1)) +
  scale_x_continuous(lim = c(-1.5, 2))

(p7b + p7) / (p8 + p9) + plot_annotation(tag_levels = "A")

# représentation des effets pentes

plot(
  ggpredict(
    glmm.tetras.slope,
    terms = c("NAOdjfm[all]", "RN"),
    type = "re"
  ),
  show_residuals = T,
  facet = T,
  jitter = F
) +
  labs(x = "NAO décembre - mars", y = "Probabilité de reproduction", title = "") +
  scale_color_viridis_d(aesthetics = c("colour", "fill")) +
  theme_classic()

# pour ne représenter que quelques régions d'intérêt :

plot(
  ggpredict(
    glmm.tetras.slope,
    terms = c("NAOdjfm [all]", "RN [Bauges, Devoluy]"),
    type = "re"
  ),
  show_residuals = T,
  facet = T,
  jitter = F
) +
  labs(x = "NAO décembre - mars", y = "Probabilité de reproduction", title = "effets aléatoires RN/UN \n sur l'ordonnée à l'origine et la pente") +
  scale_color_viridis_d(aesthetics = c("colour", "fill"))

## Encart 10.3. La complexité d’un modèle mixte---------------------------------

# pour répliquer la simulation à l'identique

set.seed(2023)

# on crée 20 quadrats

blocs <- paste("B", 1:20, sep = "")

# l'intercept de chaque quadrat

a <- abs(rnorm(20, 70, 10))

# moyenne et erreur standard sur les 20 blocs (= ce que le modèle doit retrouver)

mean(a)
sd(a) / sqrt(length(a))

# pente de chaque bloc (qu'on force à être négative)

b <- rnorm(20, -5, 10)

# pente que le modèle doit retrouver

mean(b)
sd(b) / sqrt(length(b))

# variable explicative tirée aléatoirement en chaque point (centrée réduite). On tire les mêmes valeurs dans tous les blocs

x <- scale(seq(0, 118, 4))

# création du jeu de données

y <- matrix(NA, nrow = 30, ncol = 20)
colnames(y) <- paste("A", 1:20, sep = "")
eps <- rnorm(20, 0, 20) # erreur intra-point
for (i in 1:20) {
  y[, i] <- a[i] + b[i] * x
}
y2 <- melt(y + eps)
y2$x <- rep(x, times = 20)
colnames(y2) <- c("point", "quadrat", "concentration", "profondeur")

# Les données :

ggplot(y2) +
  aes(x = profondeur, y = concentration) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ quadrat)

# Les trois modèles :

m1 <- lm(concentration ~ profondeur, data = y2) # modèle linéaire sans effet bloc
m2 <- lmer(concentration ~ profondeur + (1 | quadrat), data = y2)
m3 <- lmer(concentration ~ profondeur + (profondeur |
                                           quadrat), data = y2) # forte sensibilité dès que les pentes sont trop hétérogènes

# paramètres :

summary(m1)
summary(m2)
summary(m3)

# intervalles de confiance :

confint(m1)
confint(m2)
confint(m3)

# courbes prédites :

p1 <- plot(
  ggpredict(m1, terms = "profondeur"),
  show_data = T,
  colors = "darkblue",
  dot_alpha = 0.2
) +
  labs(title = "sans effet quadrat")

p2 <- plot(
  ggpredict(m2, terms = "profondeur"),
  show_data = T,
  colors = "darkblue",
  dot_alpha = 0.2
) +
  labs(title = "effet aléatoire quadrat sur l'intercept")

p3 <- plot(
  ggpredict(m3, terms = "profondeur"),
  show_data = T,
  colors = "darkblue",
  dot_alpha = 0.2
) +
  labs(title = "effet aléatoire quadrat sur l'intercept et pente")

(p1 + p2) / (p3 + plot_spacer())

## 10.5 Le cas des échantillonnages à stratifications croisées------------------

## Un cas de stratification croisée : modéliser les relations abondance – paysage chez plusieurs espèces

# les données de points d'écoute :

couesnon$maille <- factor(couesnon$maille)
couesnon$espece <- factor(couesnon$espece)
couesnon$point <- factor(couesnon$point)

summary(couesnon)

# la liste des espèces :

levels(couesnon$espece)

# carte des points :

p1 <-
  mapview(
    geojson,
    color = "black",
    col.regions = "#fde725",
    legend = F,
    map.types = "OpenStreetMap"
  )

# pas dans l'ouvrage car fond de carte non libre de droits : pour afficher un fond satellite :

villes <- data.frame(
  ville = c("Rennes", "Fougères", "Combourg", "Dol de Bretagne"),
  latitude = c(48.1173, 48.3538, 48.4110, 48.5486),
  longitude = c(-1.6778, -1.2045, -1.7464, -1.7499)
)

villes_sp <- sf::st_as_sf(villes, coords = c("longitude", "latitude"))

p1 <-
  mapview(
    geojson,
    color = "black",
    col.regions = "#fde725",
    legend = F,
    map.types = "OpenStreetMap"
  ) + mapview(
    villes_sp,
    color = "white",
    col.regions = "black",
    legend = F,
    alpha = 1,
    label = villes_sp$ville,
    labelOptions = labelOptions(noHide = TRUE, textOnly = TRUE)
  )

# corrélation des variables explicatives (échelle mailles)

couesnon.maille <- unique(couesnon[, c("maille",
                                       "heterogeneite",
                                       "configuration",
                                       "prairies",
                                       "haies")])
couesnon.cor <- couesnon.maille[, -1]
chart.Correlation(couesnon.cor, histogram = T, pch = 19)

# moyenne et variance des nombres d'oiseaux par point :

mean(couesnon$total)
var(couesnon$total)

# nombre de points et de mailles d'occurrence de chaque espèce :

df.sp <- data.frame(occ = tapply(
  couesnon$total,
  INDEX = couesnon$espece,
  FUN = function(x) {
    sum(x > 0)
  }
))
df.sp$espece <- rownames(df.sp)

p1 <- ggplot(df.sp) +
  aes(x =  reorder(espece, -occ), y = occ) +
  geom_bar(stat = "identity", fill = "#3b528b") +
  theme_classic() +
  labs(x = "Espèces", y = "Nombre de points d'occurrence \n (total= 104 points)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(text = element_text(size = 11))

# nombre de points d'occurrence par maille

df.sp2 <- aggregate(
  couesnon$total,
  by = list(couesnon$maille, couesnon$espece),
  FUN = function(x) {
    sum(x > 0)
  }
)
colnames(df.sp2) <- c("maille", "espece", "n_points_mailles")

p2 <- ggplot(df.sp2) +
  aes(x = reorder(espece, -n_points_mailles), y = n_points_mailles) +
  geom_boxplot() +
  theme_classic() +
  labs(x = "Espèces", y = "Nombre de points d'occurrence par maille \n (4 points par maille") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(text = element_text(size = 11)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(text = element_text(size = 11))

# nombre de mailles d'occurrence

df.sp3 <-
  aggregate(couesnon$total,
            by = list(couesnon$espece, couesnon$maille),
            FUN = sum)
colnames(df.sp3) <- c("espece", "maille", "total")

df.sp4 <- data.frame(occ = tapply(
  df.sp3$total,
  INDEX = df.sp3$espece,
  FUN = function(x) {
    sum(x > 0)
  }
))
df.sp4$espece <- rownames(df.sp4)

p3 <- ggplot(df.sp4) +
  aes(x =  reorder(espece, -occ), y = occ) +
  geom_bar(stat = "identity", fill = "#3b528b") +
  theme_classic() +
  labs(x = "Espèces", y = "Nombre de mailles d'occurrence \n (total= 26 mailles)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(text = element_text(size = 11))

(p1 + p2) / (p3 + plot_spacer()) + plot_annotation(tag_levels = "A")

## Un modèle sans stratification

# Le GLM sans aucune stratification :

glm.couesnon <-
  glm(
    total ~ poly(heterogeneite, 2) + poly(configuration, 2) + poly(prairies, 2) + poly(haies, 2),
    family = "poisson",
    data = couesnon
  )

# Les résidus:

par(mfrow = c(3, 2))
plot(glm.couesnon)

sim.glm.couesnon2 <- simulateResiduals(fittedModel = glm.couesnon)
testDispersion(sim.glm.couesnon2, alternative = "greater")

# quelques points de résidus extrêmes : ils correspondent à des groupes d'étourneaux et aux trois seules données > 16 oiseaux

couesnon[c(102, 738, 892), ]
subset(couesnon, total > 16)

# les paramètres :

summary(glm.couesnon)

# on peut regarder les courbes :

p1 <- plot(ggpredict(glm.couesnon, terms = "heterogeneite[all]")) +
  labs(x = "Hétérogénéité de composition", y = "Effectif prédit", title =
         "") +
  theme_classic()

p2 <- plot(ggpredict(glm.couesnon, terms = "configuration[all]")) +
  labs(x = "Hétérogénéité de configuration", y = "Effectif prédit", title =
         "") +
  theme_classic()

p3 <- plot(ggpredict(glm.couesnon, terms = "prairies[all]")) +
  labs(x = "Surface de prairies", y = "Effectif prédit", title = "") +
  theme_classic()

p4 <- plot(ggpredict(glm.couesnon, terms = "haies[all]")) +
  labs(x = "Linéaire de haies", y = "Effectif prédit", title = "") +
  theme_classic()

(p1 + p2) / (p3 + p4) + plot_annotation(tag_levels = "A")

# analyse de sensibilité : le même modèle sans l'étourneau. C'est mieux, mais pas parfait : il doit y avoir d'autres espèces surabondantes ou des points ou mailles à fortes abondances. Nous pourrions poursuivre l'investigation (il le faudrait en cas réel), mais pour cet exemple nous avons surtout besoin de réaliser que l'identité de l'espèce est probablement en partie responsable de la surdispersion.

couesnon2 <- subset(couesnon, espece != "etourneau_sansonnet")

glm.couesnon2 <-
  glm(
    total ~ poly(heterogeneite, 2) + poly(configuration, 2) + poly(prairies, 2) + poly(haies, 2),
    family = poisson,
    data = couesnon2
  )

sim.glm.couesnon2 <- simulateResiduals(fittedModel = glm.couesnon2)

par(mfrow = c(3, 2))
plot(glm.couesnon2)
testDispersion(sim.glm.couesnon2, alternative = "greater")

# Variation résiduelle par mailles et par espèces

par(mfrow = c(2, 2))

boxplot(residuals(glm.couesnon) ~ couesnon$maille,
        xlab = "mailles",
        ylab = "résidus du GLM")

boxplot(residuals(glm.couesnon) ~ couesnon$espece,
        xlab = "espèces",
        ylab = "résidus du GLM")

## Un modèle mixte à effets aléatoires croisés : stratifications par mailles et par espèces

# visualiser la stratification croisée

ggplot(df.sp2) +
  aes(
    x = maille,
    y = espece,
    fill =  ifelse(n_points_mailles > 0, n_points_mailles, NA)
  ) +
  scale_fill_continuous(name = "Nombre de points d'occurrence \n par maille", na.value = 'white') +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# le modèle :

glmm.couesnon <-
  glmer(
    total ~ poly(heterogeneite, 2) + poly(configuration, 2) + poly(prairies, 2) + poly(haies, 2) +
      (1 | maille) + (1 | espece),
    family = "poisson",
    data = couesnon
  )

# résidus avec la librairie performance :

check_model(glmm.couesnon)

# on voit nettement l'impact de la stratification sur l'incertitude des paramètres

summary(glmm.couesnon)

# R² : il est très faible, suggérant que l'effet des variables paysagères est très bruité. Ce n'est pas forcément étonnant, car les variables explicatives ne sont pas mesurées à la même échelle que les abondances d'oiseaux. On pourrait réfléchir à un moyen de modéliser l'abondance des oiseaux à l'échelle de la maille, et estimer l'effet des variables paysagères sur ces abondances à la maille. C'est faisable à l'aide de modèles hiérarchiques mais nous n'irons pas jusque là.

r.squaredGLMM(glmm.couesnon)

# (Pas dans l'ouvrage) On peut néanmoins tenter une amélioration en structurant les effets des variables par espèces (pentes aléatoires). Attention, pas d'effet aléatoire maille sur les pentes : ça n'aurait pas de sens puisqu'il n'y a pas de variabilité des variables explicatives à l'intérieur des mailles.
# Ce modèle commence à être un peu long à tourner et va poser des problèmes de convergence à cause d'espèces insuffisamment échantillonnées le long des gradients d'intérêt : nous en resterons donc là.

p1 <- plot(ggpredict(glmm.couesnon, terms = "heterogeneite[all]"),
           show_residuals = T) +
  labs(x = "Hétérogénéité de composition", y = "Effectif prédit", title =
         "") +
  theme_classic() +
  ylim(0, 10)
p2 <- plot(ggpredict(glmm.couesnon, terms = "configuration[all]"),
           show_residuals = T) +
  labs(x = "Hétérogénéité de configuration", y = "Effectif prédit", title =
         "") +
  theme_classic() +
  ylim(0, 10)
p3 <- plot(ggpredict(glmm.couesnon, terms = "prairies[all]"),
           show_residuals = T) +
  labs(x = "Surface de prairies", y = "Effectif prédit", title = "") +
  theme_classic() +
  ylim(0, 10)
p4 <- plot(ggpredict(glmm.couesnon, terms = "haies[all]"),
           show_residuals = T) +
  labs(x = "Linéaire de haies", y = "Effectif prédit", title = "") +
  theme_classic() +
  ylim(0, 10)

(p1 + p2) / (p3 + p4) + plot_annotation(tag_levels = "A")

# Sans les résidus partiels : c'est plus lisible, mais ça ne doit pas faire oublier le bruit résiduel autour de ces courbes lié aux espèces mal prédites, principalement par excès de 0 ou de très grandes valeurs.

p1 <- plot(ggpredict(glmm.couesnon, terms = "heterogeneite[all]")) +
  labs(x = "Hétérogénéité de composition", y = "Effectif prédit", title =
         "") +
  theme_classic()
p2 <- plot(ggpredict(glmm.couesnon, terms = "configuration[all]")) +
  labs(x = "Hétérogénéité de configuration", y = "Effectif prédit", title =
         "") +
  theme_classic()
p3 <- plot(ggpredict(glmm.couesnon, terms = "prairies[all]")) +
  labs(x = "Surface de prairies", y = "Effectif prédit", title = "") +
  theme_classic()
p4 <- plot(ggpredict(glmm.couesnon, terms = "haies[all]")) +
  labs(x = "Linéaire de haies", y = "Effectif prédit", title = "") +
  theme_classic()

(p1 + p2) / (p3 + p4) + plot_annotation(tag_levels = "A")
