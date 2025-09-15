##----------------------------------------------------------------------------##

# JY. Barnagaud et O. Gimenez. Analyse de données pour les écologues. Biotope
# Editions, 2025.

#                       -------Annexes-------                               #

## commandes R pour la reproduction des exemples des annexes
## à jour 16 mai 2025

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

pacman::p_load(# génériques
  tidyverse,
  formatR,
  jtools,
  reshape2,
  formattable,
  readxl,
  
  # graphiques
  patchwork,
  viridis,
  ggplot2,
  cowplot,
  factoextra,
  
  # statistiques
  lmtest,
  ade4)

# pour installer manuellement chaque package : install.packages(nomdupackage)
# pour charger manuellement chaque package : library(nomdupackage)

## options ---------------------------------------------------------------------

# enlève fond gris, garde grille et grossit la taille des légendes

theme_set(theme_light(base_size = 16))

# reproductibilité des simulations

set.seed(2020)

## Annexe 1 --------------------------------------------------------------------

mesanges <-
  read.csv2("donnees/mesanges_Caro_et_al_PRSLB_2019.csv", dec = ".")

# on repasse ce jeu de données dans un format "long" plus facile à manipuler

mesanges.long <- pivot_longer(mesanges, cols = c("T0", "T30"))
colnames(mesanges.long) <- c("individu", "traitement", "concentration")

## 2-	Comparer deux échantillons au regard d’une variable quantitative

# histogramme des données :

p1 <- ggplot(mesanges.long) +
  aes(x = concentration, fill = traitement) +
  geom_histogram(bins = 9, color = 'black') +
  labs(x = "Concentration en oestradiol (pg/ml)", y = "Effectif") +
  scale_fill_viridis_d() +
  theme_classic()

p1

ggsave("outputs/A1F2.png", width = 6, height = 6)

# visualisation des liens entre groupes

png(
  filename = "outputs/A1F3.png",
  width = 20,
  height = 20,
  units = "cm",
  res = 300,
  bg = 'white'
)

plot(
  c(rep(0.5, 16), rep(1.5, 16)),
  c(mesanges[, 2], mesanges[, 3]),
  xlim = c(0, 2),
  xaxt = "n",
  bty = "l",
  xlab = "échantillon",
  ylab = "concentration en oestradiol (pg/ml)",
  pch = 21,
  bg = "black",
  cex.lab = 1.5
)

axis(
  side = 1,
  at = c(0.5, 1.5),
  labels = c("T=0 min", "T=30min")
)

segments(
  x0 = rep(0.5, 16),
  x1 = rep(1.5, 16),
  y0 = mesanges[, 2],
  y1 = mesanges[, 3]
)
dev.off()

# La différence entre les deux groupes :

mesanges$delta <- apply(mesanges[, -1], FUN = "diff", MARGIN = 1)
mesanges

# On range les valeurs

mesanges <- mesanges[order(abs(mesanges$delta)), ]
mesanges$rang <- rank(abs(mesanges$delta))
mesanges

# sommes des rangs

sum(mesanges[which(mesanges$delta > 0), "rang"])
sum(mesanges[which(mesanges$delta < 0), "rang"])

# test de Wilcoxon

wilcox.test(mesanges[, 2], mesanges[, 3], paired = T)

## 3-	Comparer des fréquences

proc.migr <- read.csv2("donnees/procellariiformes.csv", dec = ".")

proc.migr$species <- factor(proc.migr$species)
proc.migr$Migratory.mode <- factor(proc.migr$Migratory.mode)
proc.migr$Nest.type <- factor(proc.migr$Nest.type)

head(proc.migr)
dim(proc.migr)
summary(proc.migr)

# table de contingence en effectifs

table(proc.migr[, -1])

# table de contingence en fréquences:

table(proc.migr[, -1]) / 68

# table de contingence théorique sous l'attendu de répartition aléatoire :

chisq.test(proc.migr$Migratory.mode, proc.migr$Nest.type)$expected

# loi de chi² à 1 degré de liberté

png(
  filename = "outputs/A1F5.png",
  width = 20,
  height = 20,
  units = "cm",
  res = 300,
  bg = 'white'
)

chisq1 <- rchisq(100, 1)
curve(
  dchisq(x, df = 1),
  from = 0,
  to = 20,
  main = "",
  bty = "n",
  xlab = "Chi² à 1 degré de liberté",
  ylab = "Densité"
)

abline(
  v = 18.7,
  col = "steelblue",
  lty = "dashed",
  lwd = 2
)
dchisq(18.77, 1)

dev.off()

# test du chi²

chisq.test(proc.migr$Migratory.mode, proc.migr$Nest.type)

## Le cas des échantillonnages déséquilibrés

passereaux <- read.csv2("donnees/traits_passereaux.csv", dec = ".")

passereaux$species <- factor(passereaux$species)
passereaux$Migratory.mode <- factor(passereaux$Migratory.mode)
passereaux$Social.system <- factor(passereaux$Social.system)

head(passereaux)
summary(passereaux)

# table de contingence : on compare les modes migratoires et la socialité

table(passereaux[, c("Migratory.mode", "Social.system")])

# on essaie un chi2

chisq.test(passereaux$Migratory.mode, passereaux$Social.system)

# l'avertissement vient du fait que les hypothèses du chi² ne sont pas respectées : plusieurs effectifs attendus sont très faibles

chisq.test(passereaux$Migratory.mode, passereaux$Social.system)$expected

# on passe par un test exact de Fisher

fisher.test(passereaux$Migratory.mode, passereaux$Social.system)

## Annexe 2 --------------------------------------------------------------------

bufo <- read.csv2("donnees/pheno_bufobufo.csv", dec = ".")

# La régression du chapitre 7

mod.bufo <- lm(date_resc ~ TmaxFev.Mars, data = bufo)
summary(mod.bufo)

# les paramètres à conserver :

pente <- -6
ordonnee_origine <- 140
sd_residuelle <- 10

# L'effectif :

taille_echantillon <- 10

# On simule des températures dans une gamme plausible, en l'occurrence entre le minimum et le maximum observés dans les données sur le crapaud commun :

variable_explicative <- seq(min(bufo$TmaxFev.Mars), max(bufo$TmaxFev.Mars), length = taille_echantillon)
variable_explicative

# On simule la variable de réponse :

moyenne_reponse <- ordonnee_origine + pente * variable_explicative

# On ajoute la variabilité :

variable_reponse <- moyenne_reponse + rnorm(n = taille_echantillon, mean = 0, sd = sd_residuelle)
variable_reponse

# On essaie maintenant de retrouver les paramètres :

reg <- lm(variable_reponse ~ variable_explicative)
summary(reg)

# On extrait la p-value :

p_value <- coef(summary(reg))["variable_explicative", "Pr(>|t|)"]

p_value

## Déterminer la taille d’échantillon nécessaire

# Nombre de simulations :

nombre_simulations <- 400

# On fait varier la taille d'échantillon de 10 à 100, par pas de 5

taille_echantillon <- seq(5, 100, by = 5)

# on prépare un vecteur pour stocker les p-values :

p_value <- rep(NA, nombre_simulations)

# et pour stocker l'indicateur de puissance :

puissance <- rep(NA, length(taille_echantillon))

## La simulation : on boucle sur les tailles d'échantillon qu'on veut tester, et pour chaque taille d'échantillon on fait 400 tirages.

# première boucle : on répète toutes les commandes sur toutes les tailles d'échantillon

for (j in 1:length(taille_echantillon)) {
  variable_explicative <- seq(min(bufo$TmaxFev.Mars), max(bufo$TmaxFev.Mars), length = taille_echantillon[j])
  
  # deuxième boucle, imbriquée dans la précédente : on répète 400 fois la simulation
  
  for (i in 1:nombre_simulations) {
    moyenne_reponse <- ordonnee_origine + pente * variable_explicative
    
    variable_reponse <- rnorm(taille_echantillon[j], mean = moyenne_reponse, sd = sd_residuelle)
    
    reg <- lm(variable_reponse ~ variable_explicative)
    
    p_value[i] <-
      coef(summary(reg))["variable_explicative", "Pr(>|t|)"]
    
  }
  
  puissance[j] <- sum(p_value < 0.05) / nombre_simulations
}

# stocker les tailles d'échantillons et les puissances dans un tableau

df <- data.frame(taille_echantillon, puissance)

# Résultats :

round(df, 2)

## La relation taille d’échantillon – puissance pour un effet faible

# La simulation : on garde tous les paramètres identiques à la simulation précédente, sauf la pente

pente <- -1
ordonnee_origine <- 140
sd_residuelle <- 10

nombre_simulations <- 400

taille_echantillon <- seq(10, 500, by = 10)

p_valeur <- rep(NA, nombre_simulations)
puissance <- rep(NA, length(taille_echantillon))

for (j in 1:length(taille_echantillon)) {
  variable_explicative <- seq(min(bufo$TmaxFev.Mars), max(bufo$TmaxFev.Mars), length = taille_echantillon[j])
  
  for (i in 1:nombre_simulations) {
    moyenne_reponse <- ordonnee_origine + pente * variable_explicative
    
    variable_reponse <- rnorm(taille_echantillon[j], mean = moyenne_reponse, sd = sd_residuelle)
    
    reg <- lm(variable_reponse ~ variable_explicative)
    
    p_valeur[i] <-
      coef(summary(reg))["variable_explicative", "Pr(>|t|)"]
  }
  puissance[j] <- sum(p_valeur < 0.05) / nombre_simulations
}

# Résultat :


data.frame(x = taille_echantillon, y = puissance) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0.8,
             color = "#3b528b",
             lty = 2) +
  labs(x = "taille d'échantillon", y = "puissance") +
  theme_classic()

ggsave("outputs/A2F1.png", width = 5, height = 5)

# Autres paramètres à étudier

pente <- c(-5, -1, -0.5, 0.5, 1, 5)
ordonnee_origine <- 140
sd_residuelle <- 10

nombre_simulations <- 400

taille_echantillon <- seq(10, 500, by = 10)

p_valeur <- rep(NA, nombre_simulations)
puissance0 <- rep(NA, length(taille_echantillon))
puissance <-
  matrix(NA, nrow = length(puissance0), ncol = length(pente))

for (k in 1:length(pente)) {
  for (j in 1:length(taille_echantillon)) {
    variable_explicative <- seq(min(bufo$TmaxFev.Mars), max(bufo$TmaxFev.Mars), length = taille_echantillon[j])
    
    for (i in 1:nombre_simulations) {
      moyenne_reponse <-
        ordonnee_origine + pente[k] * variable_explicative
      
      variable_reponse <- rnorm(taille_echantillon[j], mean = moyenne_reponse, sd = sd_residuelle)
      
      reg <- lm(variable_reponse ~ variable_explicative)
      
      p_valeur[i] <-
        coef(summary(reg))["variable_explicative", "Pr(>|t|)"]
    }
    puissance0[j] <- sum(p_valeur < 0.05) / nombre_simulations
  }
  puissance[, k] <- puissance0
}

df.k <- as.data.frame(puissance)
colnames(df.k) <- paste("pente =", pente, sep = " ")
df.k.w <- as.data.frame(pivot_longer(df.k, cols = 1:length(pente)))
df.k.w$sample.size <- rep(taille_echantillon, each = length(pente))

ggplot(df.k.w) +
  aes(x = sample.size, y = value, col = name) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0.8,
             color = "#3b528b",
             lty = 2) +
  labs(x = "taille d'échantillon", y = "puissance", col = "pente") +
  theme_classic() +
  scale_color_viridis_d()

ggsave("outputs/A2F2.png", width = 8, height = 5)

## Test de puissance pour un GLM

lomolino <- read.csv2("donnees/lomolino.csv", dec = ".")
lomolino$log.surface <- log(lomolino$surface)

# Le GLM du chapitre 8 :

m1 <- glm(
  n_especes ~ log.surface + distance_source + latitude,
  family = "poisson",
  data = lomolino
)

m2 <- glm(n_especes ~ log.surface, family = "poisson", data = lomolino)

lrtest(m1, m2)

# Pour l'analyse de puissance, on reprend les paramètres du m1 :

ordonnee_origine <- 1
pente_surf <- 0.2
pente_dsour <- -0.005
pente_lati <- 0.005

# On tente une dizaine de tailles d'échantillons :

taille_echantillon <- seq(5, 50, by = 5)
nombre_simulations <- 400

# stockage des valeurs d'intérêt

p_value <- rep(NA, nombre_simulations)
puissance <- rep(NA, length(taille_echantillon))

# variables explicatives

for (j in 1:length(taille_echantillon)) {
  var_surf <- runif(
    min = min(lomolino$surface),
    max = max(lomolino$surface),
    n = taille_echantillon[j]
  )
  log.var_surf <- log(var_surf)
  
  var_dsour <- runif(
    min = min(lomolino$distance_source),
    max = max(lomolino$distance_source),
    n = taille_echantillon[j]
  )
  
  var_lati <- runif(
    min = min(lomolino$latitude),
    max = max(lomolino$latitude),
    n = taille_echantillon[j]
  )
  
  # variable de réponse
  
  for (i in 1:nombre_simulations) {
    log_lambda <- ordonnee_origine +
      pente_surf * log.var_surf +
      pente_dsour * var_dsour +
      pente_lati * var_lati
    
    
    variable_reponse <- rpois(n = taille_echantillon[j], lambda = exp(log_lambda))
    
    # hypothèse 1
    
    glm1 <-
      glm(variable_reponse ~ var_surf + var_dsour + var_lati, family = "poisson")
    
    # hypothèse 2
    
    glm2 <- glm(variable_reponse ~ var_surf, family = "poisson")
    
    # ratio de vraisemblances
    
    test <- lrtest(glm1, glm2)
    p_value[i] <- test$`Pr(>Chisq)`[2]
    
  }
  
  puissance[j] <- sum(p_value < 0.05) / nombre_simulations
  
}


# La relation effectif - puissance :

df.glm <- data.frame(taille_echantillon, puissance)

ggplot(df.glm) +
  aes(x = taille_echantillon, y = puissance) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0.8,
             color = "#3b528b",
             lty = 2) +
  labs(x = "taille d'échantillon", y = "puissance") +
  theme_classic() +
  scale_color_viridis_d()

ggsave("outputs/A2F3.png", width = 6, height = 6)

## figures de la couverture du livre -------------------------------------------

#### Figures de la couverture du livre ####

##### première de couverture

###### table

tab.couv <- read_xlsx("donnees/reptiles_languedoc.xlsx")
tab.format <- formattable(tab.couv)
tab.format
tab.format2 <- formattable(tab.couv, align = c("l", rep("r", NCOL(tab.couv) - 1)))
tab.format2

###### densité de probabilité

bufo <- read.csv2("donnees/pheno_bufobufo.csv", dec= ".")

dens <- density(bufo$TmaxFev.Mars)
df <- data.frame(x = dens$x, y = dens$y)
probs <- c(0.1, 0.25, 0.5,  0.9)
quantiles <- quantile(bufo$TmaxFev.Mars, prob = probs)
df$quant <- factor(findInterval(df$x, quantiles))
ggplot(df, aes(x, y))+
  geom_line()+ 
  geom_ribbon(aes(ymin = 0, ymax = y, fill = quant))+
  scale_x_continuous(breaks = quantiles)+
  scale_fill_brewer(guide="none")+
  labs(x = "Températures, °C", y = "Densité de probabilité")

ggsave("outputs/cover_2.png", width = 8, height = 7, bg = "white")

###### capture R

avi <- read.csv2("donnees/avifaune_perche.csv")
summary(avi)

###### ACP
climat <- read.csv2("donnees/climat_LR_ACP.csv", row.names = 1)
pc.clim <- dudi.pca(climat, scannf = F, nf = 2)

fviz_pca_biplot(
  pc.clim,
  title = "ACP : climat d'Occitanie",
  col.ind = "steelblue",
  col.var = "black",
  label = "var",
  pointsize = 3,
  alpha.ind = 0.5,
  labelsize = 4
)+
  xlim(-4,4)+
  ylim(-2,3)
ggsave("outputs/cover_1.png", width = 7, height = 7, bg = "white")

# quatrième de couverture

mesanges <- read.csv2("donnees/mesures_nichoirs.csv", dec = ".")
ggplot(mesanges) +
  aes(x = tarse, y = poids, color = espece) +
  geom_point(size = 3, alpha = 0.7) +
  labs(x = "longueur du tarse (mm)", y = "poids (g)", title = "Mésanges, forêt d'Orléans", col = "espèce :") +
  theme_classic() +
  theme(legend.position = "top")+
  scale_color_manual( values = c("steelblue","goldenrod"), labels = c("Mésange bleue", "Mésange charbonnière"))
# ggsave("outputs/cover-fig.png", width = 5, height = 5)

