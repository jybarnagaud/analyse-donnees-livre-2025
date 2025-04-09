##----------------------------------------------------------------------------##

# JY. Barnagaud et O. Gimenez. Analyse de données pour les écologues. Biotope
# Editions, 2025.

#                       -------Chapitre 5-------                               #

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
  dplyr,
  oce,
  swfscMisc,
  
  # graphiques
  patchwork,
  viridis,
  ggplot2,
  cowplot,
  
  # spatial
  sf,
  mapview
)

# pour installer manuellement chaque package : install.packages(nomdupackage)
# pour charger manuellement chaque package : library(nomdupackage)

## options ---------------------------------------------------------------------

# enlève fond gris, garde grille et grossit la taille des légendes
theme_set(theme_light(base_size = 16))

# reproductibilité des simulations
set.seed(150)

## 5.2 Les conditions clés d’un échantillonnage bien construit -----------------

## La condition de robustesse (simulation)


# nombre de réplicas
N <- 1000

# une très grande population tirée dans une normale
pop10000 <- rnorm(10000, 0, 1)

# 5 sous-populations tirées aléatoirement dans la plus grande
pop50 <- sample(pop10000, 100)
pop200 <- sample(pop10000, 200)
pop500 <- sample(pop10000, 200)
pop1000 <- sample(pop10000, 1000)
pop5000 <- sample(pop10000, 5000)

# on échantillonne dans ces populations, de 1 à 100 individus

param.pop <- data.frame()
for (i in 5:100) {
  for (j in c(50, 200, 500, 1000, 5000, 10000)) {
    pop <- paste("pop", j, sep = "")
    tp <- sample(get(pop), i)
    param.pop <- rbind(param.pop, c(j, i, mean(tp), sd(tp)))
  }
}
colnames(param.pop) <- c("taille.pop", "n.ind", "moyenne", "sd")

# représentation graphique

p1 <- ggplot(param.pop) +
  aes(x = n.ind, y = moyenne) +
  geom_line() +
  geom_hline(
    yintercept = mean(pop10000),
    linetype = "dashed",
    col = "red"
  ) +
  facet_wrap(~ taille.pop) +
  labs(x = "Effectif de l'échantillon", y = "Moyenne") +
  theme_light(base_size = 20)

p1

p2 <- ggplot(param.pop) +
  aes(x = n.ind, y = sd) +
  geom_line() +
  geom_hline(
    yintercept = sd(pop10000),
    linetype = "dashed",
    col = "red"
  ) +
  facet_wrap( ~ taille.pop) +
  labs(x = "Effectif de l'échantillon", y = "Ecart type") +
  theme_light(base_size = 20)

p2

p1 + p2

# sauvegarde
ggsave("outputs/C5F44.png", width = 30, height = 15)

## La condition de randomisation

# Stratégie 3 : le tirage systématique

# données
marees <- read.csv2("donnees/marees_novembre_2021_shom.csv", dec = ".")

# formatage des dates
marees$date2 <- as.POSIXct(marees$date, format = "%d/%m/%Y %H:%M")

# échantillonnage toutes les 12 heures
echt <- seq(min(marees$date2), max(marees$date2), by = '12 hours')

# modèle de marées (pour interpoler les données)

dat.mod <- as.sealevel(elevation = marees$hauteur, time = marees$date2)
maree.mod <- tidem(dat.mod, latitude = 48.676)
maree.fit <- cbind(marees$date2, predict(maree.mod))

# hauteurs d'eau échantillonnées toutes les 12h
obs.tide <- NULL

for (i in 1:length(echt)) {
  maree.echt <- cbind(rep(echt[i], 2), c(0, 120))
  obs.tide <- c(obs.tide, as.numeric(crossing.point(maree.fit, maree.echt)[2]))
}
tide.obs <- data.frame(echt, obs.tide)

# représentation
ggplot(marees) +
  aes(x = date2, y = hauteur) +
  geom_line(aes(y = predict(maree.mod))) +
  labs(x = "Date", y = "Hauteur d'eau") +
  theme_classic() +
  geom_point(data = tide.obs, aes(x = echt, y = obs.tide), col = "red") +
  geom_smooth(
    data = tide.obs,
    aes(x = echt, y = obs.tide),
    col = "red",
    size = 0.5,
    se = F
  ) +
  theme_light(base_size = 16)


ggsave("outputs/C5F53b.png", width = 20, height = 10)


# richesse spécifique sur des points d'écoute

perche <- read.csv("donnees/habitat_Perche.csv",sep=";")
perche_sf <- st_as_sf(perche, coords = c("longitude", "latitude"))
perche_sf <- st_set_crs(perche_sf, 4326)

ggplot(perche) +
  aes(x = longitude, y = latitude, col = rs, size = 1.5) +
  geom_point()+
  scale_colour_viridis_c(option = "viridis") +
  theme_classic() +
    xlim(0.45, 0.62)+
    ylim(48.38, 48.42)+
  guides(size = "none") 

ggsave("outputs/C5F53a.png", width = 10, height = 10)
