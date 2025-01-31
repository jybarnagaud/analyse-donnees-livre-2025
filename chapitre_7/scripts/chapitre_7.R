##----------------------------------------------------------------------------##

# JY. Barnagaud et O. Gimenez. Analyse de données pour les écologues. Biotope
# Editions, 2025.

#                       -------Chapitre 7-------                               #

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
  patchwork,
  viridis,
  webshot,
  fGarch,
  
  # tests post-hoc de Tukey
  multcomp,
  
  # sorties graphiques de modèles
  ggeffects,
  sjPlot,
  jtools,
  effects,
  visreg
  
)

# pour installer manuellement chaque package : install.packages(nomdupackage)
# pour charger manuellement chaque package : library(nomdupackage)

## options ---------------------------------------------------------------------

# enlève fond gris, garde grille et grossit la taille des légendes
theme_set(theme_light(base_size = 16))

# reproductibilité des simulations
set.seed(2020)

## jeux de données -------------------------------------------------------------

bufo <- read.csv2("donnees/pheno_bufobufo.csv", dec = ".")
biom.mesanges <- read.csv2("donnees/biometrie_mesanges_bleues.csv", dec =
                             ".")
cistudes <- read.csv2("donnees/cistudes_Lansargues.csv", dec = ".")
climat <- read.csv2("donnees/climat_LR.csv", dec = ".")
cistudes.pop <- read.csv2("donnees/cistudes_3pops.csv", dec = ".")
amphi <- read.csv2("donnees/amphibiens.csv", dec = ".")

## 7.2	Relation entre deux variables quantitatives continues-------------------

## Les données

bufo$date <- as.Date(bufo$date, format = "%d/%m/%Y")
summary(bufo)

# graphique recherché dans la suite du chapitre

bufo %>%
  ggplot() +
  aes(x = TmaxFev.Mars, y = date_resc) +
  geom_point() +
  geom_smooth(method = lm,
              color = 'red',
              se = FALSE) +
  labs(x = "Température (Février - Mars, °C)", y = "Date de sortie d'hibernation des crapauds") +
  theme_classic()

# histogramme des dates de sortie

ggplot(bufo) +
  aes(x = date_resc) +
  geom_histogram(bins = 10,
                 color = "black",
                 fill = "white") +
  labs(x = "Date de sortie", y = "Fréquence") +
  theme_classic()

# histogramme des températures

ggplot(bufo) +
  aes(x = TmaxFev.Mars) +
  geom_histogram(bins = 10,
                 color = "black",
                 fill = "white") +
  labs(x = "Température (février-mars, °C)", y = "Fréquence") +
  theme_classic()

# relation date~température

ggplot(bufo) +
  aes(x = TmaxFev.Mars, y = date_resc) +
  geom_point(size = 2) +
  labs(x = "Température (février-mars, °C)", y = "Date de sortie") +
  theme_classic()

## L’implémentation du modèle sous R

# Le modèle linéaire avec la fonction lm :

mod.bufo <- lm(date_resc ~ TmaxFev.Mars, data = bufo)

# graphiques de vérification des résidus. On peut aussi les afficher un par un en n'appelant pas la première ligne et en tapant "Entrée" entre chaque graphique.

par(mfrow = c(2, 2))
plot(mod.bufo)

# pour reproduire les deux premiers graphiques à  la main (non représenté dans l'ouvrage): on utilise la fonction *predict*, qui extrait les valeurs prédites du modèle, et la fonction *residuals*, qui extrait les résidus :

plot(
  predict(mod.bufo),
  residuals(mod.bufo),
  xlab = "valeurs prédites",
  ylab = "résidus",
  cex.lab = 1.5,
  bty = "n",
  pch = 21,
  bg = "black"
)

qqnorm(
  scale(residuals(mod.bufo)),
  pch = 21,
  bg = "black",
  bty = "n",
  main = "",
  xlab = "quantiles théoriques",
  ylab = "quantiles observés",
  cex.lab = 1.5
)
abline(0, 1, col = "red", lty = "dashed")

# Graphiques de résidus un par un :

plot(mod.bufo, which = 1)
bufo[c(6, 7, 32), ] # points extrêmes sur ce premier graphique :

plot(mod.bufo, which = 2)

plot(mod.bufo, which = 3)

plot(mod.bufo, which = 4)

## Encart 7.3-------------------------------------------------------------------

# Dans cette simulation, nous utilisons la librairie fGarch pour générer des résidus volontairement décalés à droite ou à gauche par rapport à une loi normale centrée.

# les coefficients qu'il faudra retrouver :

alpha.sim <- 2
beta.sim <- 0.5

## A chaque itération, on simule :

# -- une variable explicative x tirée dans une loi normale
# -- un résidu eps tiré dans une loi normale asymétrique (l'asymétrie reste toujours la même afin de nous focaliser sur l'effet de la réplication aléatoire de la simulation, mais vous pouvez la diminuer ou l'augmenter : argument xi dans rsnorm)
# -- une variable de réponse à partir de x, eps, alpha.sim, beta.sim
# -- on ré-estime les paramètres à partir d'une régression linéaire et on calcule les intervalles de confiance à 95%
# -- on étudie la différence entre coefficients simulés et coefficients estimés, en fonction de si la normalité des résidus a été rejetée ou non

df.lm <- data.frame()

set.seed(150)

for (i in 1:10) {
  x <- rnorm(30, 0, 1)
  eps <- rsnorm(30, mean = 0, sd = 1, xi = 10)
  
  y <- alpha.sim + beta.sim * x + eps
  
  lm1 <- lm(y ~ x)
  
  resid.lm1 <- residuals(lm1)
  
  
  sh.test <- shapiro.test(resid.lm1)[["p.value"]]
  
  df.lm <-
    rbind(df.lm, c(
      coef(lm1)[1],
      confint(lm1)[1, ],
      coef(lm1)[2],
      confint(lm1)[2, ],
      sh.test
    ))
}

colnames(df.lm) <-
  c("alpha",
    "alpha.low",
    "alpha.high",
    "beta",
    "beta.low",
    "beta.high",
    "pval.shapiro")

df.lm$sig <- "normalité non rejetée"
df.lm[which(df.lm$pval.shapiro < 0.05), "sig"] <-
  "normalité rejetée"

p.int <- ggplot(df.lm) +
  aes(x = 1:10, y = alpha) +
  geom_point(size = 3) +
  geom_segment(aes(
    x = 1:10,
    xend = 1:10,
    y = alpha.low,
    yend = alpha.high
  )) +
  geom_hline(yintercept = alpha.sim,
             col = "#440154",
             linetype = "dashed") +
  facet_wrap(~ sig) +
  labs(x = "coefficients simulés", y = "ordonnée à l'origine") +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

p.pente <- ggplot(df.lm) +
  aes(x = 1:10, y = beta) +
  geom_point(size = 3) +
  geom_segment(aes(
    x = 1:10,
    xend = 1:10,
    y = beta.low,
    yend = beta.high
  )) +
  geom_hline(yintercept = beta.sim,
             col = "#440154",
             linetype = "dashed") +
  facet_wrap(~ sig) +
  labs(x = "coefficients simulés", y = "pente") +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

p.int / p.pente + plot_annotation(tag_levels = 'A')

## L’ajustement et l’idée de corrélation

# Corrélation nulle:

plot(
  bufo$TmaxFev.Mars,
  bufo$date_resc,
  xlab = "Température (février-mars, °C)",
  ylab = "Date de première sortie",
  cex.lab = 1.5,
  bty = "n",
  type = "n"
)
points(
  bufo$TmaxFev.Mars,
  rnorm(nrow(bufo), mean(predict(mod.bufo)), sd(predict(mod.bufo))),
  pch = 21,
  cex = 1.5,
  bg = "black"
)
abline(mod.bufo)

# Situation observée :

plot(
  bufo$TmaxFev.Mars,
  bufo$date_resc
  ,
  xlab = "Température (février-mars, °C)",
  ylab = "Date de première sortie",
  cex.lab = 1.5,
  pch = 21,
  cex = 1.5,
  bg = "black",
  bty = "n"
)
abline(mod.bufo)

# Cas théorique 2 : corrélation exacte

plot(
  bufo$TmaxFev.Mars,
  bufo$date_resc,
  xlab = "Température (février-mars, °C)",
  ylab = "Date julienne de première sortie",
  cex.lab = 1.5,
  bty = "n",
  type = "n"
)
points(
  bufo$TmaxFev.Mars,
  predict(mod.bufo),
  pch = 21,
  cex = 1.5,
  bg = "black"
)
abline(mod.bufo)

# calcul du R² à la main

ecart.estime <- (predict(mod.bufo) - mean(bufo$date_resc)) ^ 2
ecart.estime
ecart.observe <- (bufo$date_resc - mean(bufo$date_resc)) ^ 2
ecart.observe
r.carre <- sum(ecart.estime) / sum(ecart.observe)
r.carre

## Encart 7.4-------------------------------------------------------------------

# jeu de données (partie du jeu de données mesures_nichoirs.csv, exploité dans les précédents chapitres : on ne conserve que les mésanges bleues, en peuplements de chênes, afin d'éviter des effets confondants):

summary(biom.mesanges)

# figure :

ggplot(biom.mesanges) +
  aes(x = tarse, y = poids) +
  geom_point(size = 2) +
  theme_classic() +
  labs(x = "Tarse (mm)", y = "Poids (g)")

# coefficient de corrélation :

cor.test(biom.mesanges$tarse, biom.mesanges$poids)

# exemple sur les données climatiques :

summary(climat)

# données brutes :

ggplot(climat) +
  aes(x = temperature_moyenne, y = cumul_precip) +
  geom_point(size = 2) +
  theme_classic() +
  labs(x = "Température moyenne mensuelle (°C)", y = "Cumul de précipitations mensuel (mm)")

# coefficient de corrélation :

cor.test(climat$temperature_moyenne, climat$cumul_precip)

# Modèle linéaire sur les mésanges :

lm.biom <- lm(poids ~ tarse, data = biom.mesanges)
summary(lm.biom)

## Comment interpréter ce modèle ?

# Tout ce que vous avez besoin de savoir sur le modèle de régression :

summary(mod.bufo)

# Représentation graphique:

ggplot(bufo) +
  aes(x = TmaxFev.Mars, y = date_resc) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "Température (février-mars, °C)", y = "Date de sortie") +
  theme_classic()

# intervalle de confiance à 95% sur la pente :

coefficients(summary(mod.bufo))[2, 1] - qt(0.975, nrow(bufo) - 1) * coefficients(summary(mod.bufo))[2, 2]
coefficients(summary(mod.bufo))[2, 1] + qt(0.975, nrow(bufo) - 1) * coefficients(summary(mod.bufo))[2, 2]

# intervalle de confiance à 95%  sur l'ordonnée à l'origine :

coefficients(summary(mod.bufo))[1, 1] - qt(0.975, nrow(bufo) - 1) * coefficients(summary(mod.bufo))[1, 2]
coefficients(summary(mod.bufo))[1, 1] + qt(0.975, nrow(bufo) - 1) * coefficients(summary(mod.bufo))[1, 2]

# intervalle de confiance sur la prédiction du modèle (qui permet d'établir la figure 7.15 - pas dans l'ouvrage) :

xforpredIC <- seq(
  from = min(bufo$TmaxFev.Mars),
  to = max(bufo$TmaxFev.Mars),
  length.out = 100
)
predIC <- predict(mod.bufo,
                  newdata = list(TmaxFev.Mars = xforpredIC),
                  interval = "confidence")
predIC <- as.data.frame(cbind(grille = xforpredIC, predIC))

predIC %>%
  ggplot() +
  aes(x = grille, y = fit) +
  geom_line() +
  geom_line(aes(x = grille, y = lwr), linetype = 2) +
  geom_line(aes(x = grille, y = upr), linetype = 2) +
  #  geom_ribbon(aes(ymin = lwr, ymax = upr),
  #              fill = "grey70",
  #              alpha = 0.3) +
  labs(x = "Température printanière °C", y = "Dates  prédites", title = "Intervalle de confiance à 95%")

## test sur les paramètres :

# calcul de t

t.stud <- (-5.873) / 1.461
t.stud

# p-value : probabilité de t dans une distribution théorique de t à 30 degrés de liberté :

2 * pt(-abs(t.stud), df = 30)

# quantiles de la loi de student pour un degré de liberté = 32-2

rand.student <- rt(100, df = 30)
qt.crit <- quantile(rand.student, p = c(0.025, 0.975))

# Figure 7.16

dens.student <- density(rand.student)
plot(
  dens.student,
  main = "",
  xlab = "Valeurs de t",
  ylab = "Densité de probabilité",
  bty = "n",
  xlim = c(-5, 5)
)
abline(
  v = qt.crit[1],
  lty = "dashed",
  lwd = 3,
  col = "darkblue"
)
abline(
  v = qt.crit[2],
  lty = "dashed",
  lwd = 3,
  col = "darkblue"
)

## Prédire de nouvelles données

# Générer un gradient sur lequel on veut prédire

temp.new <- seq(min(bufo$TmaxFev.Mars), max(bufo$TmaxFev.Mars), by = 0.1)
temp.new

# Prédictions à la main : on applique les coefficients du modèle à notre gradient de température :

dates.predites <- coef(mod.bufo)[1] + coef(mod.bufo)[2] * temp.new
dates.predites

# Avec predict :

dates.predites.ci <- predict(mod.bufo,
                             newdata = list(TmaxFev.Mars = temp.new),
                             interval = "confidence")
head(dates.predites.ci)

# Intervalle de confiance :

predict(mod.bufo,
        newdata = list(TmaxFev.Mars = 14.5),
        interval = "confidence")

# Intervalle de confiance de prédiction :

predict(mod.bufo,
        newdata = list(TmaxFev.Mars = 14.5),
        interval = "prediction")

# Figure : comparaison entre deux types d'intervalles de confiance

pr1 <- plot(
  ggpredict(mod.bufo, terms = "TmaxFev.Mars", interval = "confidence"),
  add.data = T,
  show.title = F
) +
  labs(x = "Température (février-mars, °C)", y = "Date de sortie") +
  ylim(0, 100)

pr2 <- plot(
  ggpredict(mod.bufo, terms = "TmaxFev.Mars", interval = "prediction"),
  add.data = T,
  show.title = F
) +
  labs(x = "Température (février-mars, °C)", y = "Date de sortie") +
  ylim(0, 100)

pr1 + pr2

# On peut prédire sur plusieurs températures d'un coup :

xpred <- predict(mod.bufo,
                 newdata = list(TmaxFev.Mars = c(12, 14, 15, 16, 17)),
                 interval = "confidence")
xpred

# sur une série de données :

dates.predites.ci <- predict(mod.bufo,
                             newdata = list(TmaxFev.Mars = temp.new),
                             interval = "confidence")
head(dates.predites.ci)

# L'intervalle de prédiction est plus conservatif :

dates.predites.ci <- predict(mod.bufo,
                             newdata = list(TmaxFev.Mars = temp.new),
                             interval = "prediction")
head(dates.predites.ci)

# la figure du modèle à la main, en utilisant la prédiction (cette figure n'est pas dans l'ouvrage)

xpred.bufo <- sort(bufo$TmaxFev.Mars)
pred.bufo <- predict(mod.bufo,
                     newdata = list(TmaxFev.Mars = xpred.bufo),
                     interval = "confidence")
df.pred.bufo <- data.frame(xpred.bufo, pred.bufo)

ggplot(df.pred.bufo) +
  geom_line(aes(x = xpred.bufo, y = fit), col = "blue") +
  geom_ribbon(aes(
    x = xpred.bufo,
    ymin = lwr,
    ymax = upr,
    color = NULL
  ), alpha = .15) +
  geom_point(aes(x = TmaxFev.Mars, y = date_resc), data = bufo) +
  labs(x = "Température maximale (février-mars, °C)", y = "Date de sortie")

## Prédiction ou extrapolation

# Figure : éviter l'extrapolation (le code ne figure pas dans l'ouvrage)

# on génère une séquence de températures qui va au-delà des données observées

xmax <- max(bufo$TmaxFev.Mars)
xseq <- seq(xmax, xmax + 5, length.out = 10)

# deux scénarios possibles d'évolution de la date de sortie d'hibernation sur les températures non échantillonnées

fit.max <- df.pred.bufo[which(df.pred.bufo$xpred.bufo == xmax), "fit"]
fit.sup <- c(fit.max, fit.max + rnorm(9, 1, 5))
df.sup <- data.frame(xseq, fit.sup)
df.sup$y2 <- predict(mod.bufo, newdata = list(TmaxFev.Mars = xseq))

# représentation graphique

ggplot(df.pred.bufo) +
  geom_line(aes(x = xpred.bufo, y = fit),
            col = "#440154",
            linewidth = 1.05) +
  geom_line(
    aes(x = xseq, y = fit.sup),
    col = "#fde725",
    data = df.sup,
    linewidth = 1.05
  ) +
  geom_line(
    aes(x = xseq, y = y2),
    col = "#1fa187",
    data = df.sup,
    linewidth = 1.05
  ) +
  geom_point(
    aes(x = TmaxFev.Mars, y = date_resc),
    data = bufo,
    size = 2,
    alpha = 0.5
  ) +
  labs(x = "Température (février-mars, °C)", y = "Date de sortie") +
  xlim(min(bufo$TmaxFev.Mars), xmax + 5) +
  theme_classic()

## Encart 7.5 (le script ne figure pas dans l'ouvrage) -------------------------

# On récupère la série de températures qui nous a servi au modèle sur les crapauds

clim <- unique(bufo[, c("annee", "TmaxFev.Mars")])

# on est en 2002 : on dispose donc d'une série de températures de 1983 à 2002.

clim.sub <- subset(clim, annee %in% 1983:2002)

# le modèle de régression qui décrit la tendance de température sur 1983 - 2002

lm.temp <- lm(TmaxFev.Mars ~ annee, data = clim.sub)

# représentation graphique
clim.plot <- clim.sub %>%
  ggplot() +
  aes(x = annee, y = TmaxFev.Mars) +
  geom_line(linewidth = 1.05, col = "gray50") +
  geom_abline(
    intercept = coef(lm.temp)[1],
    slope = coef(lm.temp)[2],
    color = "#440154",
    linewidth = 1.05
  ) +
  xlim(1983, 2014) +
  labs(x = "Années", y = "Température printanière (février-mars), °C")

# on va extrapoler cette tendance au futur 2002- 2015, en conservant l'écart type observé avant 2002:

set.seed(2010)
pred <- rnorm(
  n = 13,
  mean = predict(lm.temp, newdata = list(annee =
                                           2003:2015)),
  sd = sd(clim.sub$TmaxFev.Mars)
)

# extrapolation 1

pred1 <- data.frame(x = 2002:2015, y = c(clim.sub[nrow(clim.sub), 2], pred))

# on représente le modèle, la prédiction, et la vraie série de températures :

clim.sub2 <- subset(clim, annee >= 2002)

clim.plot +
  geom_line(
    data = pred1,
    aes(x = x, y = y),
    color = "#1fa187",
    linetype = 2,
    linewidth = 1.05
  ) +
  geom_line(
    data = clim.sub2,
    aes(x = annee, y = TmaxFev.Mars),
    color = "#fde725",
    linewidth = 1.05
  ) +
  theme_classic()

## A quoi sert l’ordonnée à l’origine ?

# Le modèle:

mod.bufo.si <- lm(date_resc ~ TmaxFev.Mars - 1, data = bufo)

# résidus

par(mfrow = c(2, 2))
plot(mod.bufo.si)

# paramètres estimés :

summary(mod.bufo.si)

## Décrire un modèle de régression dans un rapport ou une publication

# Table 7.3

# librairie ggeffect :

plot(ggeffect(mod.bufo), add.data = T)

# librairie sjPlot :

plot_model(mod.bufo, type = "pred", show.data = T)

# librairie jtools :

effect_plot(mod.bufo,
            pred = TmaxFev.Mars,
            interval = T,
            plot.points = T)

# librairie visreg :

visreg(mod.bufo, xvar = "TmaxFev.Mars")

# librairie effect :

plot(effect(
  term = "TmaxFev.Mars",
  mod = mod.bufo,
  partial.residuals = T
))

# Table 7.4 (avec sjPlot)

tab_model(
  mod.bufo,
  show.stat = T,
  show.df = T,
  pred.labels = c("Intercept", "Temperature"),
  file = "outputs/C7F21.html"
)

# 7.3	relation variable quantitative continue - variable qualitative binaire----

## La question et les données : dimorphisme sexuel chez les cistudes

#  explorer les données

cistudes$Sexe <- as.factor(cistudes$Sexe)
summary(cistudes)

# relation graphique

ggplot(cistudes) +
  aes(x = Sexe, y = Longueur) +
  geom_boxplot() +
  geom_point(alpha = 0.3,
             col = "darkblue",
             size = 2) +
  labs(y = "Longueur corporelle (mm)", x = "Sexe") +
  theme_classic()

# histogramme des longueurs

ggplot(cistudes) +
  aes(x = Longueur) +
  geom_histogram(bins = 10,
                 color = "black",
                 fill = "white") +
  labs(x = "Longueur corporelle (mm)", y = "Effectif") +
  theme_classic()

# modèle

mod.cistude <- lm(Longueur ~ Sexe, data = cistudes)

# résidus

par(mfrow = c(2, 2))
plot(mod.cistude)

# résumé du modèle

summary(mod.cistude)

# résultats formatés

tab_model(
  mod.cistude,
  show.stat = T,
  show.df = T,
  pred.labels = c("Intercept", "Temperature"),
  file = "outputs/C7F25.html"
)

# une sortie graphique parmi de nombreuses possibles

plot_model(mod.cistude, terms = "Sexe", type = "eff")

# une autre (représentation facile des points de données) (Figure 7.24)

plot(
  ggeffect(mod.cistude),
  add.data = T,
  dot.alpha = 0.2,
  show.title = F
)

## Encart 7.6. Modalité de variable qualitative manquante-----------------------

# simulation (le code pour simuler les données n'est pas dans l'ouvrage)

set.seed(2020)

# une variable x faite de A et de B

x <- as.factor(rep(c("A", "B"), 30))

# on la passe en numérique

x.num <- as.numeric(x) - 1
x.num

# valeurs des paramètres de régression

alpha <- 2
beta <- 0.5

# résidu

eps <- rnorm(30, 0, 5)

# jeu de données simulé

y <- alpha + beta * x.num + eps
df.sim <- data.frame(x, y)
head(df.sim)

# le modèle :

lm.sim <- lm(y ~ x , data = df.sim)
coef(lm.sim)

# valeur moyenne de A

mean(df.sim[which(x == "A"), "y"])

# valeur moyenne de B

mean(df.sim[which(x == "B"), "y"])

# ordonnée à l'origine + beta :

sum(coef(lm.sim))

# sans l'ordonnée à l'origine :

lm.sim.2 <- lm(y ~ -1 + x, data = df.sim)
coef(lm.sim.2)

## Le test de Student et le test de Welch

# calculer la statistique t de Student à la main

moy.cistude <- tapply(cistudes$Longueur, INDEX = cistudes$Sexe, FUN = "mean")

beta.cistude = moy.cistude[2] - moy.cistude[1]

effectif.sexes <-
  tapply(cistudes$Longueur, INDEX = cistudes$Sexe, FUN = "length")

var.intrasexes <-
  tapply(cistudes$Longueur, INDEX = cistudes$Sexe, FUN = "var")

t.cistude <-
  beta.cistude / sqrt((var.intrasexes[1] / effectif.sexes[1]) + (var.intrasexes[2] /
                                                                   effectif.sexes[2]))

# loi de student pour un degré de liberté = n-2

rand.student <- rt(100, df = 175)
dens.student <- density(rand.student)

plot(
  dens.student,
  main = "",
  xlab = "Valeurs de t",
  ylab = "Densité de probabilité",
  bty = "n"
)

qt.crit = quantile(rand.student, p = c(0.025, 0.975))

abline(
  v = qt.crit[1],
  lty = "dashed",
  col = "darkblue",
  lwd = 2
)

abline(
  v = qt.crit[2],
  lty = "dashed",
  col = "darkblue",
  lwd = 2
)

# calcul de la p-value (males différents que femelles = test bilatéral)

2 * pt(-abs(t.cistude), df = 175)

# calcul de la p-value (males plus petits que femelles = test unilatéral)

pt(-abs(t.cistude), df = 175)

# Faire le test t (ou le test de Welch) sous R : il y a deux possibilités.

# Soit on a deux échantillons séparés que l'on compare :

cis.F <- subset(cistudes, Sexe == "F")
cis.M <- subset(cistudes, Sexe == "M")

t.test(cis.M$Longueur, cis.F$Longueur) # test de Welch
t.test(cis.M$Longueur, cis.F$Longueur, var.equal = T) # test de Student

# Soit, comme dans l'ouvrage, on passe par une formule similaire à celle des modèles linéaires :

t.test(Longueur ~ Sexe, data = cistudes)
t.test(Longueur ~ Sexe, data = cistudes, var.equal = T)

## Encart 7.7-------------------------------------------------------------------

# densité de probabilité d'une loi de Student à 175 degrés de liberté

rand.student <- rt(100, df = 175)
dens.student <- density(rand.student)

# comparaison des seuils à 5% dans un test bilatéral et unilatéral

par(mfrow = c(1, 2))

# test bilatéral

plot(
  dens.student,
  main = "",
  xlab = "Valeurs de t",
  ylab = "Densité de probabilité",
  bty = "n"
)
qt.crit = quantile(rand.student, p = c(0.025, 0.975))
abline(
  v = qt.crit[1],
  lty = "dashed",
  col = "darkblue",
  lwd = 2
)
abline(
  v = qt.crit[2],
  lty = "dashed",
  col = "darkblue",
  lwd = 2
)

# test unilatéral

plot(
  dens.student,
  main = "",
  xlab = "Valeurs de t",
  ylab = "Densité de probabilité",
  bty = "n"
)
qt.crit.uni = quantile(rand.student, p = 0.05)
abline(
  v = qt.crit.uni,
  lty = "dashed",
  col = "orange",
  lwd = 2
)

# p-value d'un test bilatéral (nb : la p-value est un peu différente de la fonction t.test à cause d'un ajustement des degrés de liberté)

2 * pt(-abs(t.cistude), df = 175)

# p-value d'un test unilatéral

pt(-abs(t.cistude), df = 175)

# fonction pour le test t bilatéral

t.test(cis.M$Longueur, cis.F$Longueur)

# fonction pour le test t unilatéral

t.test(cis.M$Longueur, cis.F$Longueur, alternative = "less")

## L’analyse de variance (ANOVA)

ggplot(cistudes) +
  aes(x = Sexe, y = Longueur) +
  geom_boxplot() +
  labs(y = "Longueur corporelle (mm)", x = "Sexe") +
  theme_classic()

## test de Fisher

# longueurs moyennes des deux sexes

moy.sexe <- tapply(cistudes$Longueur, INDEX = cistudes$Sexe, FUN = "mean")

# longueur moyenne tous individus confondus

y.moy <- mean(cistudes$Longueur)

# effectif par sexes

n.sexe <- tapply(cistudes$Longueur, INDEX = cistudes$Sexe, FUN = "length")

# somme des écarts intersexes

ecart.femelles <- moy.sexe[1] - y.moy
ecart.males <- moy.sexe[2] - y.moy
somme.ecarts.facteur <-
  sum((n.sexe[1] * (ecart.femelles ^ 2)), (n.sexe[2] * (ecart.males ^ 2)))

# somme des écarts intrasexes

ecart.residu.F <-
  cistudes[which(cistudes$Sexe == "F"), "Longueur"] - moy.sexe[1]
ecart.residu.M <-
  cistudes[which(cistudes$Sexe == "M"), "Longueur"] - moy.sexe[2]
somme.ecarts.residu <- sum((ecart.residu.F ^ 2), (ecart.residu.M ^ 2))

# statistique F

F.stat <-
  (somme.ecarts.facteur / (2 - 1)) / (somme.ecarts.residu / (nrow(cistudes) -
                                                               2))
F.stat

# loi de Fisher

plot(
  density(rf(1000, 2 - 1, 177 - 2)),
  xlab = "Valeurs du F de Fisher à 1 et 175 ddl",
  ylab = "Densité de probabilité",
  main = "",
  bty = "n",
  cex.lab = 1.5
)

# p value

pf(F.stat, 2 - 1, 177 - 2, lower.tail = F)

## Encart 7.8-------------------------------------------------------------------

# Simulation :

set.seed(800)

# on tire une variable x pour de nombreux individus

x <- rnorm(1000, 0, 1)

# une magnitude d'effet forte (pour chaque incrément de 1 x, on multiplie y par 5!!):

alpha <- 2
beta <- 5

# un résidu élevé : la relation y ~ x existe bien, mais elle est environnée de beaucoup de variation

eps <- rnorm(1000, 0, 10)

# simulation des données

y <- alpha + beta * x + eps
df <- data.frame(x, y)

#  maintenant, nous allons essayer de retrouver alpha et beta avec différentes tailles d'échantillons, et regarder la significativité du F à chaque effectif :

n <- c(5, 10, 15, 20, 25, 50, 100)
res.sim <- data.frame()

for (i in 1:length(n)) {
  samp <- sample(1:nrow(df), n[i])
  
  df.i <- df[samp, ]
  
  lm.i <- lm(y ~ x, data = df.i)
  anova.i <- anova(lm.i)
  
  F.val <- anova.i[["F value"]][1]
  F.pval <- anova.i[["Pr(>F)"]][1]
  coef.a <- coef(lm.i)[1]
  coef.b <- coef(lm.i)[2]
  conf.b1 <- confint(lm.i)[2, 1]
  conf.b2 <- confint(lm.i)[2, 2]
  
  res.sim <- rbind(res.sim, c(n[i], F.val, F.pval, coef.a, coef.b, conf.b1, conf.b2))
}

colnames(res.sim) <- c("n",
                       "F.val",
                       "p.value",
                       "alpha.est",
                       "beta.est",
                       "beta.ic.low",
                       "beta.ic.high")

# une colonne qui nous donne la significativité statistique (seuil 5%)

res.sim$sig <- "grey"
res.sim[which(res.sim$p.value < 0.05), "sig"] <- "black"

# représentations graphiques

ggplot(df) +
  aes(x = x, y = y) +
  geom_point() +
  labs(x = "x", y = "y") +
  theme_classic()

ggplot(res.sim) +
  aes(x = n, y = beta.est, col = sig) +
  geom_hline(yintercept = beta, linetype = "dashed") +
  geom_segment(aes(
    x = n,
    y = beta.ic.low,
    xend = n,
    yend = beta.ic.high,
    col = sig
  )) +
  geom_point(size = 2) +
  labs(x = "Effectif de l'échantillon", y = "Pente estimée ± IC à 95%") +
  scale_color_identity() +
  theme_classic()

# 7.4 Relation variable quantitative  - variable  à >2 catégories---------------

## Les données : comparaison biométrique entre trois populations de tortues

# données

cistudes.pop$Population <- factor(cistudes.pop$Population)
summary(cistudes.pop)

# résumé graphique :

p1 <- ggplot(cistudes.pop) +
  aes(x = Longueur) +
  geom_histogram(bins = 10,
                 color = "black",
                 fill = "white") +
  labs(x = "Longueur corporelle (mm)", y = "Effectif") +
  theme_classic()

p2 <- ggplot(cistudes.pop) +
  aes(x = Population, y = Longueur) +
  geom_boxplot() +
  geom_point(col = "steelblue",
             alpha = 0.1,
             size = 1.5) +
  labs(x = "Population" , y = "Longueur corporelle (mm)") +
  theme_classic()

p1 + p2

## le modèle

# Modèle :

mod.cistude.pop <- lm(Longueur ~ Population, data = cistudes.pop)

# on vérifie les résidus:

par(mfrow = c(2, 2))
plot(mod.cistude.pop)

# résumé :

summary(mod.cistude.pop)

# intervalle de confiance à 95% :

confint(mod.cistude.pop)

# ANOVA :

anova(mod.cistude.pop)

# Représentation graphique des effets :

plot(
  ggeffect(mod.cistude.pop),
  add.data = T,
  show.title = F,
  dot.alpha = 0.1
)+
  labs(x = "Population", y = "Longueur (mm)") +
  theme_classic()

# table de résultats

tab_model(
  mod.cistude.pop,
  show.stat = T,
  show.df = T,
  pred.labels = c("Intercept", "Population Brenne","Population  Vigueirat"),
  file = "outputs/C7F32.html"
)

## Le test post-hoc de Tukey sur les contrastes entre catégories

# contraste entre Brenne et Vigueirat :

coef(mod.cistude.pop)["PopulationMarais_du_Vigueirat"] - coef(mod.cistude.pop)["PopulationBrenne"]

# test de Tukey avec la librairie multcomp :

tukey <- glht(mod.cistude.pop , linfct = mcp(Population = "Tukey"))
summary(tukey)

# test post-hoc de Tukey avec la librairie stats (chargée automatiquement à l'ouverture de R). Il y a des différences mineures entre les deux fonctions à cause d'ajustements légèrement différents, qui ne vous impacteront pas :
# on doit d'abord reformuler le modèle avec la fonction aov() (fonction traditionnelle pour les ANOVA)

mod.cistude.popb <- aov(Longueur ~ Population, data = cistudes.pop)

# le test :

posthoc <- TukeyHSD(x = mod.cistude.popb, conf.level = 0.95)
posthoc

## 7.5	Construire un modèle linéaire à plus de deux variables------------------

## Les données : comparer les sorties d’hibernation de plusieurs espèces

# Les données :

amphi$codesp <- factor(amphi$codesp)
summary(amphi)

# Représentation graphique de la variation temporelle de dates de sortie d'hibernation entre espèces

ggplot(amphi) +
  geom_line(aes(x = an, y = date_resc, col = codesp),
            linewidth = 1,
            alpha = 0.7) +
  xlab("Années") +
  ylab("Date de sortie") +
  scale_color_viridis_d(
    name = "Espèce",
    labels = c(
      "Alyte accoucheur (ALYOBS)",
      "Crapaud commun (BUFBUF)",
      "Rainette méridionale (HYLMER)",
      "Grenouille rousse (RANTEM)",
      "Triton marbré (TRIMAR)"
    )
  ) +
  theme_classic()

# distribution des dates de sortie d'hibernation (pas dans l'ouvrage)

ggplot(amphi) +
  aes(x = date_resc) +
  geom_histogram(bins = 10,
                 color = "black",
                 fill = "white") +
  labs(x = "Dates de sortie", y = "Effectif")

# effet espèces

p1 <- ggplot(amphi) +
  geom_boxplot(aes(x = codesp, y = date_resc)) +
  xlab("espèces") +
  ylab("Dates de sortie") +
  theme_classic()

# effet altitude

p2 <- ggplot(amphi) +
  geom_point(aes(x = alt, y = date_resc, color = codesp),
             size = 1.5,
             alpha = 0.7) +
  xlab("altitude (m)") +
  ylab("Dates de sortie") +
  scale_color_viridis_d(
    name = "Espèce",
    labels = c(
      "Alyte accoucheur (ALYOBS)",
      "Crapaud commun (BUFBUF)",
      "Rainette méridionale (HYLMER)",
      "Grenouille rousse (RANTEM)",
      "Triton marbré (TRIMAR)"
    )
  ) +
  theme_classic()

# relation altitude - année :

p3 <- ggplot(amphi) +
  geom_point(aes(x = an, y = alt, color = codesp),
             size = 1.5,
             alpha = 0.7) +
  xlab("années") +
  ylab("altitude (m)") +
  scale_color_viridis_d(
    name = "Espèce",
    labels = c(
      "Alyte accoucheur (ALYOBS)",
      "Crapaud commun (BUFBUF)",
      "Rainette méridionale (HYLMER)",
      "Grenouille rousse (RANTEM)",
      "Triton marbré (TRIMAR)"
    )
  ) +
  theme_classic() +
  theme(legend.position = "none")

(p1 + p2) / (p3 + plot_spacer()) + plot_annotation(tag_levels = 'A')

## Un premier modèle à une seule variable : tendance temporelle toutes espèces confondues

# régression linéaire simple

lm.an <- lm(date_resc ~ an, data = amphi)
summary(lm.an)

plot(
  ggeffect(lm.an),
  add.data = T,
  show.title = F
) +
  labs(x = "Années", y = "Dates de sortie") +
  theme_classic()

## Construire, estimer et interpréter un modèle à plusieurs variables

# régression linéaire multiple

lm.mult <- lm(date_resc ~ an + alt + codesp, data = amphi)

# résidus :

par(mfrow = c(2, 2))
plot(lm.mult)

# résumé :

summary(lm.mult)

## Quelle est la différence concrète entre plusieurs régressions linéaires simples et une régression linéaire multiple ?

# tables des régressions simples :
# sur les années :

tab_model(
  lm.an,
  show.se = T,
  show.stat = T,
  show.df = T,
  file = "outputs/C7F38a.html"
)

# sur l'altitude :

lm.alt <- lm(date_resc ~ alt, data = amphi)
tab_model(
  lm.alt,
  show.se = T,
  show.stat = T,
  show.df = T,
  file = "outputs/C7F38b.html"
)

# sur les espèces :

lm.esp <- lm(date_resc ~ codesp, data = amphi)
tab_model(
  lm.esp,
  show.se = T,
  show.stat = T,
  show.df = T,
  file = "outputs/C7F38c.html"
)

# régression multiple :

tab_model(
  lm.mult,
  show.se = T,
  show.stat = T,
  show.df = T,
  file = "outputs/C7F38d.html"
)

## Représenter les résultats d’une régression linéaire multiple

# graphique incorrect pour les années (pour illustration de ce qu'il ne faut pas faire, à ne pas réutiliser)

p1 <- ggplot(amphi) +
  aes(x = an, y = date_resc) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_abline(intercept = coef(lm.mult)[1], slope = coef(lm.mult)[2]) +
  xlab("Années") + ylab("Dates de sortie") +
  theme_classic()

# graphique incorrect pour les altitudes (pour illustration de ce qu'il ne faut pas faire, à ne pas réutiliser)

p2 <- ggplot(amphi) +
  aes(x = alt, y = date_resc) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_abline(intercept = coef(lm.mult)[1], slope = coef(lm.mult)[3]) +
  xlab("Altitudes (m)") + ylab("Dates de sortie") +
  theme_classic()

p1 + p2

# Droite marginale pour une année moyenne :

ggplot(amphi) +
  aes(x = alt, y = date_resc) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_abline(intercept = coef(lm.mult)[1] - 0.65 * mean(amphi$an),
              slope = coef(lm.mult)[3]) +
  xlab("Altitudes (m)") + ylab("Dates de sortie") +
  theme_classic()

# Graphique correct : représenter les résidus partiels

# Avec ggeffects:

p1 <- ggpredict(lm.mult, terms = "an")
p2 <- ggpredict(lm.mult, terms = "alt")
p3 <- ggpredict(lm.mult, terms = "codesp")

p1b <-
  plot(
    p1,
    show_residuals = T,
    dot_size = 1.5,
    dot_alpha = 0.2
  ) + labs(x = "Années", y = "Dates de sortie marginales", title = "effet marginal des années") +
  theme_classic()

p2b <-
  plot(
    p2,
    show_residuals = T,
    dot_size = 1.5,
    dot_alpha = 0.2
  ) + labs(x = "Altitudes (m)", y = "Dates de sortie marginale", title = "effet marginal de l'altitude") +
  theme_classic()

p3b <-
  plot(
    p3,
    show_residuals = T,
    dot_size = 1.5,
    dot_alpha = 0.2
  ) + labs(x = "Espèces", y = "Dates de sortie marginale", title = "contrastes marginaux entre espèces") +
  theme_classic()

(p1b + p2b) / (p3b + plot_spacer())


# autre solution, avec jtools :

p1 <-
  effect_plot(
    lm.mult,
    pred = an,
    interval = TRUE,
    partial.residuals = TRUE
  ) +
  xlab("Années") + ylab("Dates de sortie marginales")

p2 <-
  effect_plot(
    lm.mult,
    pred = alt,
    interval = TRUE,
    partial.residuals = TRUE
  ) +
  xlab("Altitudes (m)") + ylab("Dates de sortie marginales")

p3 <-
  effect_plot(
    lm.mult,
    pred = codesp,
    interval = TRUE,
    partial.residuals = TRUE
  ) +
  xlab("Espèces") + ylab("Dates de sortie marginales")

p1 + p2 + p3

# Avec sjPlot (essentiellement un clone de ggeffects, plus automatisé)

plot_model(lm.mult, terms = "an", type = "pred")

plot_model(lm.mult,
           terms = "an",
           type = "pred",
           show.data = T) # attention, représente les données brutes et non les résidus partiels (donc incorrect dans la majorité des cas)!

plot_model(lm.mult,
           terms = "alt",
           type = "pred",
           show.data = T) # attention, représente les données brutes et non les résidus partiels (donc incorrect dans la majorité des cas)!

# avec visreg (permet des graphiques en R de base ou en ggplot, nous montrons la version base)
# pratique car très flexible (de nombreux types de modèles implémentés dont modèles mixtes, régressions quantiles, etc), principale alternative à ggeffects

par(mfrow = c(2, 2))
visreg(lm.mult, xvar = "an")
visreg(lm.mult, xvar = "alt")
visreg(lm.mult, xvar = "codesp")

# Avec effects : autre solution en R de base, proche de ggeffects

eff.mod <- effect(term = "an", mod = lm.mult)
plot(eff.mod) # sans les résidus partiels
eff.mod2 <- predictorEffects(lm.mult, residuals = TRUE)
plot(eff.mod2, partial.residuals = list(smooth = F)) # avec les résidus partiels

## Comparer les magnitudes d’effets dans un modèle linéaire

# la fonction scale() permet de centrer réduire (ou de juste centrer, ou juste réduire, selon les arguments choisis)
# A noter qu'on ne centre-réduit pas la variable codesp qui est un facteur. C'est néanmoins possible mais plus complexe et pas forcément intéressant pour notre objectif ici.

lm.mult.sc <- lm(date_resc ~ scale(an) + scale(alt) + codesp, data = amphi)
summary(lm.mult.sc)

# représentation utile pour comparer des effets de variables centrées réduites, avec sjPlot

plot_model(lm.mult.sc,
           terms = c("scale(an)", "scale(alt)"),
           type = "est") +
  theme_classic() +
  labs(x = "Variables", y = "Coefficients ± intervalles de confiance à 95%", title =
         "")
