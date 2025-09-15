##----------------------------------------------------------------------------##

# JY. Barnagaud et O. Gimenez. Analyse de données pour les écologues. Biotope
# Editions, 2025.

#                       -------Chapitre 4-------                               #

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

# le paquet "ggradar" (dédié à l'un des graphiques) doit être chargé à partir d'un github car il n'est
# pas (encore?) inclus dans un paquet du CRAN. Vous devez donc l'installer manuellement :

install.packages("devtools")
devtools::install_github("ricardo-bion/ggradar", dependencies = TRUE)

# si vous n'effectuez pas cette installation, le script fonctionnera correctement, à l'exception de la
# ligne qui utilise ggradar.

# utiliser cette fonction pour installer et charger automatiquement les packages

install.packages("pacman")
pacman::p_load(
  # génériques
  tidyverse,
  formatR,
  dplyr,
  forcats,
  readr,
  
  # graphiques
  ggplot2,
  patchwork,
  viridis,
  ggmosaic,
  rgl,
  plotly,
  fmsb,
  
  ggradar,
  scales,
  cowplot,
  reticulate,
  mapview,
  webshot,
  
  # autres
  fGarch
  
)

# pour installer manuellement chaque package : install.packages(nomdupackage)
# pour charger manuellement chaque package : library(nomdupackage)

## options ---------------------------------------------------------------------

# enlève fond gris, garde grille et grossit la taille des légendes
theme_set(theme_light(base_size = 16))

# reproductibilité des simulations
set.seed(2020)

## jeux de données -------------------------------------------------------------

mesanges <- read.csv2("donnees/mesures_nichoirs.csv", dec = ".")
oiseaux <- read.csv2("donnees/oiseaux_couesnon.csv")
temp.04 <- read.csv2("donnees/serie_annuelle_temperature.csv")
xy.couesnon <- read.csv2("donnees/xy_couesnon.csv", dec = ".")
climat <- read.csv2("donnees/climat_LR.csv")
reptiles <- read.csv2("donnees/reptiles_LR.csv", dec =
                        ".")
petrel <- read.csv2("donnees/Bonadonna_Mardon_2010.csv", dec =
                      ".")
arbres <- read.csv2("donnees/hauteurs_arbres.csv", dec = ".")
truites <- read.csv2("donnees/truites.csv", dec =
                       ".")
dchen <- read.csv2("donnees/JD_chenilles.csv", dec = ".")
oiseaux2 <- read.csv2("donnees/points_ecoute_couesnon.csv")

## 4.1 Effectifs, fréquences, quantiles ----------------------------------------

## Effectifs et fréquences

mesanges$espece <- factor(mesanges$espece)
mesanges$essence_dominante <- factor(mesanges$essence_dominante)

# résumé

summary(mesanges)

# Nombre de lignes :

nrow(mesanges)

# Nombre de données pour chaque modalité:

summary(mesanges$espece)

# En fréquence:

summary(mesanges$espece) / nrow(mesanges)

# fréquence des petits et grands :

petits <- subset(mesanges, poids < 10)
frequence.petits <- nrow(petits) / nrow(mesanges)
frequence.petits

grands <- subset(mesanges, poids > 18)
frequence.grands <- nrow(grands) / nrow(mesanges)
frequence.grands

## Descripteurs marginaux:

# minimum:

min(mesanges$poids)

# maximum:

max(mesanges$poids)

# étendue:

max(mesanges$poids) - min(mesanges$poids)

# fréquences par classes

classes <- cut(mesanges$poids, 10)
head(classes)
levels(classes)
summary(classes) / length(classes)

# La même chose avec tapply:

tapply(
  classes,
  INDEX = classes,
  FUN = function(x) {
    length(x) / length(classes)
  }
)

## Quantiles

poids.sort <- sort(mesanges$poids)
poids.sort

# à 0.025

0.025 * nrow(mesanges)
poids.sort[10]

# à 0.975

0.975 * nrow(mesanges)
poids.sort[376]
nrow(mesanges) - 376 + 9

# étendue de 95% des poids

poids.sort[375] - poids.sort[10]

# fonction pour les quantiles

quantile(mesanges$poids)

# quelques quantiles des tarses

quantile(mesanges$tarse)

# La fonction quantile ne fonctionnera pas s'il y a des valeurs manquantes. Dans ce cas, on dispose de l'argument 'na.rm' qui dit à R quel comportement adopter en cas de valeur manquante. Cet argument est commun à de très nombreuses fonctions. Le plus souvent, on choisit 'T' (la fonction s'exécute en ne tenant pas compte des NA) ou 'F' (la fonction renvoie un message d'erreur en cas de NA, c'est souvent l'argument par défaut). Il existe parfois d'autres possibilités.

quantile(mesanges$tarse, na.rm = T)

# On peut choisir les quantiles qu'on veut:

quantile(mesanges$poids, p = c(0, 0.025, 0.5, 0.975, 1))
quantile(mesanges$tarse,
         p = c(0, 0.025, 0.5, 0.975, 1),
         na.rm = T)

## 4.3 Représentations graphiques d’une variable -------------------------------

## Encart 4.1

ggplot(mesanges) +
  aes(x = tarse, y = poids) +
  geom_point()

# sauvegarde d'un plot (la taille s'ajuste avec les arguments height et width, en points, cm ou pixels)

ggsave("outputs/C4Encart1_1.png",
       height = 5,
       width = 5)

# nb : la fonction de sauvegarde n'est plus répétée dans les commandes suivantes codant des graphiques,
# il vous suffit de la réutiliser si vous voulez faire un export. Vous pouvez aussi cliquer le bouton
# "Export" dans la fenêtre "Plots" de R Studio.

ggplot(mesanges) +
  aes(x = tarse, y = poids) +
  geom_point() +
  geom_smooth()

# Etiquettes des axes :

ggplot(mesanges) +
  aes(x = tarse, y = poids) +
  geom_point() +
  labs(x = "longueur du tarse (mm)", y = "poids (g)", title = "Mésanges, forêt d'orléans")

# taille des points

ggplot(mesanges) +
  aes(x = tarse, y = poids) +
  geom_point(size = 3, alpha = 0.2) +
  labs(x = "longueur du tarse (mm)", y = "poids (g)", title = "Mésanges, forêt d'orléans")

# apparence générale :

ggplot(mesanges) +
  aes(x = tarse, y = poids) +
  geom_point(size = 3, alpha = 0.2) +
  labs(x = "longueur du tarse (mm)", y = "poids (g)", title = "Mésanges, forêt d'orléans") +
  theme_classic()

# couleur des points (librairie viridis):

ggplot(mesanges) +
  aes(x = tarse, y = poids, color = espece) +
  geom_point(size = 3, alpha = 0.2) +
  labs(x = "longueur du tarse (mm)", y = "poids (g)", title = "Mésanges, forêt d'orléans") +
  theme_classic() +
  scale_color_viridis_d()

# même diagramme, adapté à la couverture du livre

ggplot(mesanges) +
  aes(x = tarse, y = poids, color = espece) +
  geom_point(size = 3, alpha = 0.7) +
  labs(x = "longueur du tarse (mm)", y = "poids (g)", title = "Mésanges, forêt d'Orléans", col = "espèce :") +
  theme_classic() +
  theme(legend.position = "top")+
  scale_color_manual( values = c("steelblue","goldenrod"), labels = c("Mésange bleue", "Mésange charbonnière"))
ggsave("outputs/cover-fig.png", width = 5, height = 5)

# facetting (un graphique par catégorie d'une variable facteur)

ggplot(mesanges) +
  aes(x = tarse, y = poids, color = espece) +
  geom_point(size = 3, alpha = 0.2) +
  labs(x = "longueur du tarse (mm)", y = "poids (g)", title = "Mésanges, forêt d'orléans") +
  theme_classic() +
  scale_color_viridis_d() +
  facet_wrap( ~ essence_dominante)

## Le diagramme en barre : représenter une variable qualitative

# diagramme en barres

summary(oiseaux)

# Les espèces les plus représentées sur la période 1 :

oiseaux.p1 <- subset(oiseaux, periode == 1)

ggplot(oiseaux.p1) +
  aes(x = espece) +
  geom_bar()

# Pas très lisible : commençons par retourner les étiquettes de l'axe des abscisses

ggplot(oiseaux.p1) +
  aes(x = espece) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 90, size = 8))

# c'est toujours difficile à lire. On peut réorganiser ce graphique de l'espèce la plus commune à la plus rare :

ggplot(oiseaux.p1) +
  aes(x = fct_infreq(espece)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 90, size = 8))

# Améliorons un peu ce graphique :

ggplot(oiseaux.p1) +
  aes(x = fct_infreq(espece)) +
  geom_bar(fill = "#3b528b",
           col = "#3b528b",
           alpha = 0.2) +
  theme(axis.text.x = element_text(angle = 90, size = 8)) +
  labs(x = "Espèces", y = "Effectif")

## L’histogramme de fréquences

#  histogramme de fréquences

ggplot(mesanges) +
  aes(x = poids) +
  geom_histogram(fill = "white", col = "black") +
  labs(x = 'Poids (g)', y = 'Effectif')

# modifier le nombre de classes

h1 <- ggplot(mesanges) +
  aes(x = poids) +
  geom_histogram(fill = "white",
                 col = "black",
                 bins = 10) +
  labs(x = 'Poids (g)', y = 'Effectif', title = "10 classes")

h2 <- ggplot(mesanges) +
  aes(x = poids) +
  geom_histogram(fill = "white",
                 col = "black",
                 bins = 20) +
  labs(x = 'Poids (g)', y = 'Effectif', title = "20 classes")

h3 <- ggplot(mesanges) +
  aes(x = poids) +
  geom_histogram(fill = "white",
                 col = "black",
                 bins = 30) +
  labs(x = 'Poids (g)', y = 'Effectif', title = "30 classes")

h4 <- ggplot(mesanges) +
  aes(x = poids) +
  geom_histogram(fill = "white",
                 col = "black",
                 bins = 40) +
  labs(x = 'Poids (g)', y = 'Effectif', title = "40 classes")

(h1 + h2) / (h3 + h4)

## 4.4 Représentations graphiques du croisement entre deux variables -----------

# Deux variables quantitatives : Le nuage de points

ggplot(mesanges) +
  aes(x = tarse, y = poids) +
  geom_point(alpha = 0.2, size = 3) +
  labs(x = "Longueur du tarse (mm)", y = "Poids (g)", title = "")

# Une variable quantitative et une variable qualitative : La boîte à moustaches

ggplot(mesanges) +
  aes(x = espece, y = poids) +
  geom_boxplot() +
  labs(x = "Espèce", y = "Poids (g)", title = "") +
  scale_x_discrete(labels = c("Mésange bleue", "Mésange charbonnière"))

# Deux variables qualitatives : la table de contingences

nichees <- unique(mesanges[, c("nichoir", "espece", "essence_dominante")])
table.nichees <- table(nichees$espece, nichees$essence_dominante)

bp1 <- ggplot(as.data.frame(table.nichees)) +
  aes(Var1, Freq, fill = Var2) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Espèce", y = "Fréquence", fill = "") +
  scale_x_discrete(labels = c("Mésange bleue", "Mésange charbonnière")) +
  scale_fill_viridis_d(direction = -1) +
  theme(legend.position = "none")


bp2 <- ggplot(as.data.frame(table.nichees)) +
  aes(Var1, Freq, fill = Var2) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Espèce", y = "Fréquence", fill = "") +
  scale_x_discrete(labels = c("Mésange bleue", "Mésange charbonnière")) +
  scale_fill_viridis_d(direction = -1)

bp1 + bp2

ggplot(data = as.data.frame(table.nichees)) +
  geom_mosaic(aes(
    weight = Freq,
    x = product(Var1),
    fill = Var2
  )) +
  labs(x = "Espèce", y = "Fréquence", fill = "")

# directement sur les données d'origine (non incluse dans l'ouvrage)

ggplot(data = nichees) +
  geom_mosaic(aes(x = product(essence_dominante, espece), fill = essence_dominante)) +
  labs(x = "Espèce", y = "Fréquence", fill = "")

## 4.5 Représentations graphiques du croisement entre plus de deux variables ----

## Deux variables quantitatives et une variable qualitative

ggplot(mesanges) +
  aes(x = tarse, y = poids, color = espece) +
  geom_point(alpha = 0.4, size = 2) +
  labs(x = "Tarse (mm)", y = "Poids (g)", title = "") +
  scale_color_manual(
    labels = c("Mésange bleue", "Mésange charbonnière"),
    name = "Espèce",
    values = c("#440154", "#fde725")
  ) +
  theme(legend.position = "top")

# histogramme de plusieurs variables

ggplot(mesanges) +
  aes(x = tarse, fill = espece) +
  geom_histogram(col = "black") +
  scale_fill_manual(
    labels = c("Mésange bleue", "Mésange charbonnière"),
    name = "Espèce",
    values = c("#440154", "#fde725")
  ) +
  labs(x = "Tarse (mm)", y = "Effectif") +
  theme(legend.position = "top")

# graphique à facettes

labs <- c("Mésange bleue", "Mésange charbonnière")
names(labs) <- c("PARCAE", "PARMAJ")
ggplot(mesanges) +
  aes(x = tarse, y = poids) +
  geom_point(alpha = 0.4, size = 2) +
  labs(x = "Tarse (mm)", y = "Poids (g)") +
  facet_grid(~ espece, labeller = labeller(espece = labs))

# graphique à facettes (autre manière de coder la même chose, non incluse dans l'ouvrage)

mesanges %>%
  mutate(espece = fct_recode(
    espece,
    "Mésange bleue" = "PARCAE",
    "Mésange charbonnière" = "PARMAJ"
  )) %>% # recode les niveaux du facteur espece
  ggplot() +
  aes(x = tarse, y = poids) +
  geom_point() +
  facet_wrap(~ espece)

# histogramme à facettes

labs <- c("Mésange bleue", "Mésange charbonnière")
names(labs) <- c("PARCAE", "PARMAJ")
ggplot(mesanges) +
  aes(x = tarse) +
  geom_histogram(fill = "white", col = "black") +
  labs(x = "Tarse (mm)", y = "Effectif") +
  facet_grid( ~ espece, labeller = labeller(espece = labs))

# autre possibilité avec un "pipe" (une manière de coder spécifique au tidyverse)

mesanges %>%
  mutate(espece = fct_recode(
    espece,
    "Mésange bleue" = "PARCAE",
    "Mésange charbonnière" = "PARMAJ"
  )) %>% # recode les niveaux du facteur espece
  ggplot() +
  aes(x = tarse) +
  geom_histogram() +
  facet_wrap(~ espece)

## Trois variables quantitatives

nichees2 <- mesanges %>%
  group_by(nichoir, espece) %>% # regroupe par nichoir et espece
  summarise(
    tarse_moyen = mean(tarse),
    # tarse moyen
    poids_moyen = mean(poids),
    # poids moyen
    nb.jeunes = n()
  ) # effectif


ggplot(nichees2) +
  aes(x = tarse_moyen, y = poids_moyen, size = nb.jeunes) +
  geom_point(alpha = 0.4) +
  labs(x = 'Tarse (mm)', y = 'Poids (g)', size = 'Nombre \n de jeunes')

# autre version, avec un gradient de couleur

ggplot(nichees2) +
  aes(x = tarse_moyen, y = poids_moyen, color = nb.jeunes) +
  geom_point(alpha = 0.4, size = 5) +
  labs(x = 'Tarse (mm)', y = 'Poids (g)', color = 'Nombre \n de jeunes')

# plot 3d avec rgl (ne marche pas sur tous les systèmes d'exploitation, non inclus dans l'ovurage)

plot3d(nichees2$tarse, nichees2$poids, nichees2$nb.jeunes)

# plot3d avec plotly

plot_ly(
  data = nichees2,
  x = ~ tarse_moyen,
  y = ~ poids_moyen,
  z = ~ nb.jeunes,
  type = "scatter3d",
  mode = "markers"
) %>%
  layout(scene = list(
    xaxis = list(title = "Tarse (mm)"),
    yaxis = list(title = "Poids (g)"),
    zaxis = list(title = "Nombre de jeunes")
  ))

## Lignes et courbes

# préparation des données climatiques :

temp.04.mpl <- subset(temp.04, POSTE == "Montpellier")

# lignes et courbes, avec ou sans superposition des points de données. Il est vivement recommandé de toujours représenter les points afin de bien mettre en valeur la résolution temporelle des données.

p1 <- ggplot(temp.04.mpl) +
  aes(x = mois, y = temperature_moyenne) +
  geom_point(size = 3) +
  labs(x = "Mois (année 2004)", y = "Température moyenne (°C)") +
  theme_classic()

p2 <- ggplot(temp.04.mpl) +
  aes(x = mois, y = temperature_moyenne) +
  geom_line(linewidth = 1.5, col = "steelblue") +
  labs(x = "Mois (année 2004)", y = "Température moyenne (°C)") +
  scale_x_continuous(breaks = 1:12) +
  theme_classic()

p3 <- ggplot(temp.04.mpl) +
  aes(x = mois, y = temperature_moyenne) +
  geom_line(linewidth = 1.5, col = "steelblue") +
  geom_point(size = 3) +
  labs(x = "Mois (année 2004)", y = "Température moyenne (°C)") +
  scale_x_continuous(breaks = 1:12) +
  theme_classic()

p4 <- ggplot(temp.04.mpl) +
  aes(x = mois, y = temperature_moyenne) +
  geom_smooth(linewidth = 1.5,
              col = "#21918c",
              method = "loess") +
  geom_smooth(linewidth = 1.5,
              col = "#440154",
              method = "gam") +
  geom_smooth(linewidth = 1.5,
              col = "#fde725",
              method = "lm") +
  geom_point(size = 3) +
  labs(x = "Mois (année 2004)", y = "Température moyenne (°C)") +
  scale_x_continuous(breaks = 1:12) +
  theme_classic()

(p1 + p2) / (p3 + p4) + plot_annotation(tag_levels = 'A')

# Une utilisation inappropriée d'une ligne

rs.p1 <-
  aggregate(oiseaux.p1$espece,
            by = list(oiseaux.p1$maille),
            FUN = "length")
colnames(rs.p1) <- c("maille", "rs")

p1 <- ggplot(rs.p1) +
  aes(x = maille, y = rs) +
  geom_line(linewidth = 1.2) +
  labs(x = "Mailles", y = "Richesse spécifique") +
  theme_classic()

p2 <- ggplot(rs.p1) +
  aes(x = maille, y = rs) +
  geom_bar(stat = "identity") +
  labs(x = "Mailles", y = "Richesse spécifique") +
  theme_classic()

(p1 / p2) + plot_annotation(tag_levels = 'A')

# Les cartes

head(xy.couesnon)
xy.oiseaux <- merge(rs.p1, xy.couesnon, by = "maille", all = F)

map1 <- mapview(
  xy.oiseaux,
  xcol = "X" ,
  ycol = "Y",
  zcol = "rs",
  crs = "epsg:2154",
  map.types = "OpenStreetMap"
)

# pour sauvegarder une capture de cette carte :
# mapshot(map1, file = "outputs/C4F20.png")

## 4.6 Les représentations graphiques à éviter ---------------------------------

# Radar (ou diagramme de Kivat, ou diagramme en étoile) :

# Avec fmsb (figure 4.21):

xclim <- aggregate(
  climat$temperature_moyenne,
  by = list(climat$POSTE, climat$annee),
  FUN = "mean"
)

colnames(xclim) <- c("ville", "annee", "temperature")

xclim2 <- reshape(xclim,
                  direction = "wide",
                  idvar = "ville",
                  timevar = "annee")

xclim3 <- xclim2[, -1]

rownames(xclim3) <- xclim2[, 1]
colnames(xclim3) <- 2004:2014

xclim3 <- rbind(rep(17, ncol(xclim3)), rep(13, ncol(xclim3)), xclim3)

radarchart(xclim3, title = "Variations de température moyenne annuelle (°C)")

legend(
  "topright",
  legend = rownames(xclim3)[-c(1, 2)],
  col = 1:4,
  lty = "dashed"
)

radarchart(xclim3, title = "Variations de température moyenne annuelle (°C)")
legend(
  "topright",
  legend = rownames(xclim3)[-c(1, 2)],
  col = 1:4,
  lty = "dashed"
)

# autre radar (nécessite le package ggradar à installer depuis un github, voir en tête de script

xclim <- aggregate(
  climat$temperature_moyenne,
  by = list(climat$POSTE, climat$annee),
  FUN = "mean"
)

colnames(xclim) <- c("ville", "annee", "temperature")

xclim2 <- reshape(xclim,
                  direction = "wide",
                  idvar = "ville",
                  timevar = "annee")

colnames(xclim2) <- c("group", 2004:2014)

xclim2[, -1] <- apply(xclim2[, -1], 2, rescale)

ggradar(xclim2) +
  labs(title = "Variations de température moyenne annuelle (°C)")

# Diagramme circulaire, ou graphique en secteurs :

rept.pie <- tapply(reptiles$codesp, INDEX = reptiles$codesp, FUN = "length")

rept.df <- data.frame(espece = names(rept.pie), effectif =
                        rept.pie)

ggplot(rept.df) +
  aes(x = "",
      y = effectif,
      fill = reorder(espece, -effectif)) +
  geom_bar(width = 1,
           stat = "identity",
           color = "white") +
  coord_polar("y", start = 0) +
  labs(fill = "Espèce", x = "", y = "") +
  theme_void() +
  scale_fill_viridis_d()

# figures 3D faites sous Excel (ce script ne fait que générer les données pour les figures 23 et 24)

petrel2 <- c(as.character(petrel$occupant_S1),
             as.character(petrel$occupant_S2))

resume.petrel <- summary(factor(petrel2))

resume.petrel2 <- data.frame(names(resume.petrel), resume.petrel)

write.table(resume.petrel2,
            file = "outputs/petrels_for_3Dbarplot.txt",
            row.names = F,
            sep = "\t")

xrept <- tapply(reptiles$codesp, INDEX = reptiles$codesp, FUN = "length")

write.table(
  data.frame(espece = names(xrept), effectif = xrept),
  file = "outputs/exemple_pie_3D.txt",
  row.names = F,
  sep = "\t"
)

# Figure de résumé (Fig. 4.25)

x <- rnorm(1000, 0, 1)
x <- sort(x)
y <- 0.5 + 2 * x + rnorm(1000, 0, 2)
z <- y + log(x + 1000)
w <- factor(c(rep("A", 250), rep("B", 500), rep("C", 250)))
v <- factor(c(rep("F", 200), rep("G", 300), rep("H", 100), rep("I", 400)))

df <- data.frame(x, y, z, w, v)

# histogramme

ggplot(df) +
  aes(x = y) +
  geom_histogram(fill = "white", col = "black") +
  theme_classic()

# plot

ggplot(df) +
  aes(x = x, y = y) +
  geom_point(alpha = 0.2, size = 2) +
  theme_classic()

# boxplot

ggplot(df) +
  aes(x = w, y = y) +
  geom_boxplot() +
  theme_classic()

# barplot

ggplot(df) +
  aes(x = w) +
  geom_bar(fill = "black") +
  theme_classic()

# barplot 2

ggplot(df) +
  aes(x = w, fill = v) +
  geom_bar(position = "nudge") +
  scale_fill_viridis_d(direction = -1) +
  theme_classic()

## 4.7 Descripteurs quantitatifs d’une variable --------------------------------

## Descripteurs de tendance centrale

head(arbres)

# histogramme :

ggplot(arbres) +
  aes(x = hmax) +
  geom_histogram(fill = "white", col = "black") +
  labs(x = "Hauteur des arbres (m)", y = "Effectif") +
  theme_classic()

# variabilité de part et d'autre de la moyenne

xmean <- mean(arbres$hmax)

color <- ifelse(arbres$hmax > xmean, '#fde725' , '#3b528b')

ggplot(arbres) +
  aes(x = 1:nrow(arbres), y = hmax) +
  geom_hline(yintercept = xmean, linetype = "dashed") +
  geom_segment(aes(
    x = 1:nrow(arbres),
    y = rep(xmean, nrow(arbres)),
    xend = 1:nrow(arbres),
    yend = hmax
  ), colour = color) +
  geom_point() +
  labs(x = "Arbres", y = "Hauteur (m)") +
  theme(axis.text.x = element_blank(), axis.ticks.x =
          element_blank()) +
  theme_classic()

# moyenne

sum(arbres$hmax) / length(arbres$hmax)
mean(arbres$hmax)

# moyenne pondérée (données simulées)

weighted.mean(c(18, 14, 13, 16.5, 12.3, 8.5), w = c(7, 2, 4, 4, 1, 1))

# truites (données simulées)

truites

# moyenne des longueurs de truites

mean(truites$taille_mm)

# sans les plus grandes

mean(subset(truites, taille_mm < 230)$taille_mm)

# médiane

median(truites$taille_mm)

# table ("tableau croisé")

table(truites$taille_mm)

# Variation de l'écart entre moyenne et médiane selon la forme de la distribution d'une variable

skew.df <-
  data.frame(
    C1 <-
      rnorm(100, 0, 1),
    C2 <-
      rsnorm(100, 0, 1, 5),
    C3 <-
      rlnorm(100, meanlog = 0, sdlog = 0.5),
    C4 <- c(rnorm(80, 0, 1), rnorm(20, 100, 5))
  )

p1 <- ggplot(skew.df) +
  aes(x = C1) +
  geom_histogram(fill = "white", col = "black") +
  labs(x = "x", y = "Effectif", title = "Loi normale") +
  theme_classic() +
  geom_vline(
    xintercept = mean(skew.df$C1),
    lwd = 1.5,
    linetype = 2,
    color = "darkred"
  ) +
  geom_vline(
    xintercept = median(skew.df$C1),
    lwd = 1.5,
    linetype = 2,
    color = "steelblue"
  )

p2 <- ggplot(skew.df) +
  aes(x = C2) +
  geom_histogram(fill = "white", col = "black") +
  labs(x = "x", y = "Effectif", title = "Loi normale asymétrique") +
  theme_classic() +
  geom_vline(
    xintercept = mean(skew.df$C2),
    lwd = 1.5,
    linetype = 2,
    color = "darkred"
  ) +
  geom_vline(
    xintercept = median(skew.df$C2),
    lwd = 1.5,
    linetype = 2,
    color = "steelblue"
  )

p3 <- ggplot(skew.df) +
  aes(x = C3) +
  geom_histogram(fill = "white", col = "black") +
  labs(x = "x", y = "Effectif", title = "Loi lognormale") +
  theme_classic() +
  geom_vline(
    xintercept = mean(skew.df$C3),
    lwd = 1.5,
    linetype = 2,
    color = "darkred"
  ) +
  geom_vline(
    xintercept = median(skew.df$C3),
    lwd = 1.5,
    linetype = 2,
    color = "steelblue"
  )

p4 <- ggplot(skew.df) +
  aes(x = C4) +
  geom_histogram(fill = "white", col = "black") +
  labs(x = "x", y = "Effectif", title = "Loi normale, \n 20 individus extrêmes") +
  theme_classic() +
  geom_vline(
    xintercept = mean(skew.df$C4),
    lwd = 1.5,
    linetype = 2,
    color = "darkred"
  ) +
  geom_vline(
    xintercept = median(skew.df$C4),
    lwd = 1.5,
    linetype = 2,
    color = "steelblue"
  )

(p1 + p2) / (p3 + p4) + plot_annotation(tag_levels = 'A')

## Descripteurs de variance

#Moyenne et médiane

xmean <- mean(arbres$hmax)

xmed <- median(arbres$hmax)

ggplot(arbres) +
  aes(x = hmax) +
  geom_histogram(fill = "white", col = "black") +
  geom_segment(
    aes(
      x = 15,
      y = 25,
      xend = xmed,
      yend = 15
    ),
    arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
    col = "steelblue"
  ) +
  geom_text(aes(x = 15, y = 26, label = "médiane"),
            color = "steelblue",
            size = 5) +
  geom_segment(
    aes(
      x = 32.5,
      y = 25,
      xend = xmean,
      yend = 15
    ),
    arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
    col = "darkred"
  ) +
  geom_text(aes(x = 32.5, y = 26, label = "moyenne"),
            color = "darkred",
            size = 5) +
  labs(x = "Hauteur des arbres (m)", y = "Effectif") +
  theme_classic()

# étendue

range(arbres$hmax)
max(arbres$hmax)  - min(arbres$hmax)

# histogramme avec étendue

xrange <- range(arbres$hmax)

ggplot(arbres) +
  aes(x = hmax) +
  geom_histogram(fill = "white", col = "black") +
  geom_segment(
    aes(
      x = 15,
      y = 25,
      xend = xmed,
      yend = 15
    ),
    arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
    col = "steelblue"
  ) +
  geom_text(aes(x = 15, y = 26, label = "médiane"),
            color = "steelblue",
            size = 5) +
  geom_segment(
    aes(
      x = 32.5,
      y = 25,
      xend = xmean,
      yend = 15
    ),
    arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
    col = "darkred"
  ) +
  geom_text(aes(x = 32.5, y = 26, label = "moyenne"),
            color = "darkred",
            size = 5) +
  labs(x = "Hauteur des arbres (m)", y = "Effectif") +
  theme_classic() +
  geom_segment(
    x = xrange[1],
    xend = xrange[2],
    y = -1,
    yend = -1,
    color = "orange3",
    arrow = arrow(
      length = unit(0.2, "cm"),
      type = "closed",
      ends = "both"
    )
  ) +
  geom_text(aes(x = 15, y = -2, label = "étendue"),
            color = "orange3",
            size = 5) +
  labs(x = "Hauteur des arbres (m)", y = "Effectif")

# écart à la moyenne

xecart <- data.frame(hauteur = arbres$hmax, moyenne = xmean)
xecart$ecart <- xecart[, 1] - xecart[, 2]
xecart[1:10, ]

# écart au carré

xecart$ecart.carre <- xecart$ecart ^ 2
xecart[1:10, ]

# écart moyen au carré

mean(xecart$ecart.carre)

# écart type et erreur standard

sdev <- sd(arbres$hmax)
std.err <- sd(arbres$hmax) / sqrt(nrow(arbres))

# indicateurs de variance sur l'histogramme

ggplot(arbres) +
  aes(x = hmax) +
  geom_histogram(fill = "white", col = "black") +
  geom_segment(
    aes(
      x = 15,
      y = 25,
      xend = xmed,
      yend = 15
    ),
    arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
    col = "steelblue"
  ) +
  geom_text(aes(x = 15, y = 26, label = "médiane"),
            color = "steelblue",
            size = 5) +
  geom_segment(
    aes(
      x = 32.5,
      y = 25,
      xend = xmean,
      yend = 15
    ),
    arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
    col = "darkred"
  ) +
  geom_text(aes(x = 32.5, y = 26, label = "moyenne"),
            color = "darkred",
            size = 5) +
  labs(x = "Hauteur des arbres (m)", y = "Effectif") +
  theme_classic() +
  geom_segment(
    x = xrange[1],
    xend = xrange[2],
    y = -1,
    yend = -1,
    color = "orange3",
    arrow = arrow(
      length = unit(0.2, "cm"),
      type = "closed",
      ends = "both"
    )
  ) +
  geom_text(aes(x = 15, y = -2, label = "étendue"),
            color = "orange3",
            size = 5) +
  labs(x = "Hauteur des arbres (m)", y = "Effectif") +
  geom_segment(
    x = xmean - sdev / 2,
    xend = xmean + sdev / 2,
    y = 11.5,
    yend = 11.5,
    color = "slateblue3",
    arrow = arrow(
      length = unit(0.2, "cm"),
      type = "closed",
      ends = "both"
    )
  ) +
  geom_text(
    aes(x = 18.5, y = 11.5, label = "écart type"),
    color = "slateblue3",
    size = 5
  ) +
  geom_segment(
    x = xmean - std.err / 2,
    xend = xmean + std.err / 2,
    y = 9.5,
    yend = 9.5,
    color = "royalblue1",
    arrow = arrow(
      length = unit(0.2, "cm"),
      type = "closed",
      ends = "both"
    )
  ) +
  geom_text(
    aes(x = 19.5, y = 9.5, label = "erreur standard"),
    color = "royalblue1",
    size = 5
  ) +
  labs(x = "Hauteur des arbres (m)", y = "Effectif")

## Encart 4.1 : moyenne géométrique (données simulées) -------------------------

# taux de croissance constant

N0 <- 10
N3 <- 12
lambda.geom <- (12 / 10) ^ (1 / 3)
lambda.geom

# taux de croissance variable

lambda <- c(0.87, 1.01, 0.98)

N3 <- N0 * lambda[1] * lambda[2] * lambda[3]
N3

# moyenne arithmétique des taux de croissance et projection sur 3 ans

lambda.arith <- mean(lambda)
lambda.arith
N3.arith <- N0 * lambda.arith * lambda.arith * lambda.arith
N3.arith

# moyenne géométrique des taux de croissance

lambda.geom <- prod(lambda) ^ (1 / 3)

lambda.geom
N3.geom <- N0 * lambda.geom * lambda.geom * lambda.geom
N3.geom

# relation entre la moyenne arithmétique et la moyenne géométrique

lambda.geom <- exp(mean(log(lambda.arith)))
lambda.geom

# qu'est - ce que ça change ?

N0 <- c(5, 40, 120, 1000)

# premier vecteur de lambda (3 taux de croissance annuelle)

lambda <- c(1.01, 0.98, 1.07)
lambda.arith <- mean(lambda)
lambda.geom <- exp(mean(log(lambda)))

N3.arith <- NULL
N3.geom <- NULL

for (i in 1:length(N0)) {
  N3.arith[i] <- N0[i] * lambda.arith * lambda.arith * lambda.arith
  N3.geom[i] <- N0[i] * lambda.geom * lambda.geom * lambda.geom
}

N.sim <- data.frame(N0, N3.arith, N3.geom)

# premier vecteur de lambda (taux de croissance annuelle plus faibles)

lambda2 <- c(0.8, 0.3, 0.4)
lambda.arith2 <- mean(lambda2)
lambda.geom2 <- exp(mean(log(lambda2)))

N3.arith2 <- NULL
N3.geom2 <- NULL

for (i in 1:length(N0)) {
  N3.arith2[i] <- N0[i] * lambda.arith2 * lambda.arith2 * lambda.arith2
  N3.geom2[i] <- N0[i] * lambda.geom2 * lambda.geom2 * lambda.geom2
}

N.sim2 <- data.frame(N0, N3.arith2, N3.geom2)

N.sim2

## 4.8 Distributions de variables ----------------------------------------------

## La loi des grands nombres

# simulation

xarbres <- tapply(arbres$hmax, INDEX = arbres$nichoir, FUN = "mean")

moy <- NULL
ect <- NULL
sp <- NULL

for (i in c(3, 5, 7, 10, 15, 30, 45, 60, 96)) {
  tmp <- sample(xarbres, i)
  names(tmp) <- rep(paste0("n = ", i), i)
  moy <- c(moy, rep(mean(tmp), i))
  ect <- c(ect, rep(sd(tmp), i))
  sp <- c(sp, tmp)
}

df <-
  data.frame(
    hauteur = sp,
    sample = names(sp),
    moyenne = moy,
    ecartype = ect
  )

df %>%
  ggplot() +
  aes(x = hauteur) +
  geom_histogram(fill = "white", col = "black") +
  geom_text(
    data = df %>% group_by(sample),
    aes(
      x = 40,
      y = 12,
      label = paste0(round(moyenne, 1), " (", round(ecartype, 1), ")")
    ),
    hjust = 1,
    vjust = 1
  ) +
  labs(x = "Hauteur (m)", y = "Effectif") +
  xlim(10, 40) +
  geom_vline(aes(xintercept = moyenne),
             lty = "dashed",
             color = "red") +
  facet_wrap( ~ fct_relevel(sample, paste0("n = ", sort(
    parse_number(levels(factor(df$sample)))
  )))) +
  theme_classic()

# variation de la moyenne et de l'écart-type

moy.plot <- ggplot() +
  aes(x = c(3, 5, 7, 10, 15, 30, 45, 60, 96), y = unique(moy)) +
  geom_line() +
  geom_point(size = 3, color = "white") +
  geom_point(size = 1, color = "black") +
  labs(x = "Effectif", y = "Moyenne de l'échantillon") +
  theme_classic()

sd.plot <- ggplot() +
  aes(x = c(3, 5, 7, 10, 15, 30, 45, 60, 96), y = unique(ect)) +
  geom_line() +
  geom_point(size = 3, color = "white") +
  geom_point(size = 1, color = "black") +
  labs(x = "Effectif", y = "Ecart - type de l'échantillon") +
  theme_classic()

moy.plot | sd.plot

# histogrammes avec les variations de moyenne - écart type

moy <- NULL
ect <- NULL
ech.45.b = data.frame()

for (i in 1:100) {
  sp = sample(xarbres, 45)
  ech.45.b = rbind(ech.45.b, sp)
}
mn.45 <- apply(ech.45.b, 2, "mean")

moy <- NULL
ect <- NULL
ech.45 <- NULL

for (i in 1:9) {
  tmp <- sample(xarbres, 45)
  names(tmp) <- rep(paste0("n = ", i), 45)
  moy <- c(moy, rep(mean(tmp), 45))
  ect <- c(ect, rep(sd(tmp), 45))
  ech.45 <- c(ech.45, tmp)
}

df <- data.frame(
  hauteur = ech.45,
  sample = names(ech.45),
  moyenne = moy,
  ecartype = ect
)

df %>%
  ggplot() +
  aes(x = hauteur) +
  geom_histogram(fill = "white", col = "black") +
  geom_text(
    data = df %>% group_by(sample),
    # groupe par panneau
    aes(
      x = 40,
      y = 7.5,
      # en haut à droite
      label = paste0(round(moyenne, 1), " (", round(ecartype, 1), ")")
    ),
    # annotation differente pour chaque panneau
    hjust = 1,
    vjust = 1
  ) +
  labs(x = "Hauteur (m)", y = "Effectif") +
  xlim(10, 40) +
  geom_vline(aes(xintercept = moyenne),
             lty = "dashed",
             color = "red") + # ajoute une verticale en la moyenne
  facet_wrap( ~ sample) +
  theme_classic()

## Centrer réduire une variable

# centrage

x1 <- sample(xarbres, 9, replace = F)
x2 <- xarbres - mean(xarbres)
mean(x2)
sd(x2)

# réduction

x3 <- (xarbres - mean(xarbres)) / sd(xarbres)
mean(x3)
sd(x3)

# histogramme :

infmed <- subset(xarbres, xarbres < 23.7)
supmed <- subset(xarbres, xarbres >= 23.7)
df <- data.frame(hauteur = c(infmed, supmed),
                 seuil = c(rep('infmed', length(infmed)), rep('supmed', length(supmed))))

ggplot(df) +
  aes(x = hauteur, fill = seuil) +
  geom_histogram(color = "white") +
  scale_fill_grey()  +
  labs(title = "non centré-réduit", x = "Hauteurs (m)", y = "Effectif") +
  theme_classic() +
  theme(legend.position = "none")

infmeds <- scale(infmed)
supmeds <- scale(supmed)

dfs <- data.frame(hauteur = c(infmeds, supmeds),
                  seuil = c(rep('infmeds', length(infmeds)), rep('supmeds', length(supmeds))))

ggplot(dfs) +
  aes(x = hauteur, fill = seuil) +
  geom_histogram(color = "white") +
  scale_fill_grey()  +
  labs(title = "centré-réduit", x = "Hauteurs (m)", y = "Effectif") +
  theme_classic() +
  theme(legend.position = "none")

## 4.9 Les principales distributions de variables-------------------------------

## La loi normale

ggplot(arbres) +
  aes(x = hmax) +
  geom_histogram(fill = "white", col = "black") +
  labs(x = "Hauteurs d'arbres (m)", y = "Effectif") +
  geom_text(aes(
    x = 30,
    y = 34,
    label = paste0("moyenne = ", round(mean(hmax), 2))
  ), color = "darkred") +
  geom_text(aes(
    x = 30,
    y = 32,
    label = paste0("médiane = ", round(median(hmax), 2))
  ), color = "steelblue") +
  geom_text(aes(
    x = 30,
    y = 30,
    label = paste0("écart - type = ", round(sd(hmax), 2))
  ), color = "slateblue3") +
  theme_classic()

moy <- NULL
ect <- NULL
k <- 0

sarbres <- scale(arbres$hmax)

repeat {
  sp <- sample(arbres$hmax, 40)
  moy <- c(moy, mean(sp))
  ect <- c(ect, sd(sp))
  k <- k + 1
  if (k == 100)
    break
}

ggplot(data.frame(moy = moy)) +
  aes(x = moy) +
  geom_histogram(fill = "white", col = "black") +
  ylim(0, 10) +
  labs(x = "Moyenne sur 40 arbres", y = "Effectif") +
  geom_text(aes(x = 25, y = 7.5, label = "moyenne = 24.2 m"), color = "darkred") +
  geom_text(aes(x = 25, y = 7, label = "médiane = 24.3 m"), color = "steelblue") +
  geom_text(aes(x = 25, y = 6.5, label = "écart type = 0.5 m"), color = "slateblue3")

sd(moy)
sd(arbres$hmax) / sqrt(40)

## La loi binomiale

head(dchen)

# diagramme en barres

dchen %>%
  select(site, nbpins, nbattaq) %>%
  filter(nbattaq > 19) %>% # filtre avec nb attaques > 19
  arrange(desc(nbpins)) %>% # range par ordre decroissant selon le nombre de pins
  head(15) %>% # selectionne les 15 parcelles avec le plus de pins
  pivot_longer(cols = -site,
               names_to = "type",
               values_to = "count") %>% # format long
  
  ggplot() +
  aes(x = site, y = count, fill = type) +
  geom_col(position = "dodge") + # cote a cote
  scale_fill_manual(
    values = c("black", "grey"),
    # couleur de la legende
    labels = c("Nombre de pins infestés", "Nombre de pins"),
    name = ""
  ) +
  labs(x = "Parcelle", y = "Nombre de pins") +
  theme_classic()

# histogramme

pattaq <- data.frame(pattaq = dchen$nbattaq / dchen$nbpins)
head(pattaq)

ggplot(pattaq) +
  aes(x = pattaq) +
  geom_histogram(fill = "white", col = "black") +
  labs(x = "Probabilité d'infestation", y = "Effectif") +
  theme_classic()

# variation de la probabilité d'infestation

prob <- NULL

for (i in c(3:135)) {
  k <- 1
  p <- NULL
  repeat {
    if (k == 10)
      break
    xsamp <- sample(1:135, i)
    tp <- dchen[xsamp, ]
    p <- c(p, tp$nbattaq / tp$nbpins)
    k <- k + 1
  }
  prob <- c(prob, mean(p))
  prob0 <- data.frame(n = 1:length(prob), prob = prob)
}

ggplot(prob0) +
  aes(x = n, y = prob) +
  geom_line() +
  geom_hline(
    yintercept = mean(dchen$nbattaq / dchen$nbpins),
    linetype = "dashed",
    color = "red"
  ) +
  labs(x = "Nombre de parcelles", y = "Probabilité d'infestation") +
  theme_classic()

# histogramme

k <- 1
p <- NULL

repeat {
  if (k == 100)
    break
  xsamp <- sample(1:135, 40)
  tp <- dchen[xsamp, ]
  p <- c(p, mean(tp$nbattaq / tp$nbpins))
  k <- k + 1
}

df <- data.frame(p = p)

ggplot(df) +
  aes(x = p) +
  geom_histogram(fill = "white", col = "black") +
  labs(x = "Probabilité d'infestation moyenne sur 40 placettes", y = "Effectif")

sd(p)
sd(dchen$nbattaq / dchen$nbpins)
sd(dchen$nbattaq / dchen$nbpins) / sqrt(50)

## La loi de Poisson

pouillot <- subset(oiseaux2, espece == "pouillot veloce")
head(pouillot)

ggplot(pouillot) +
  aes(x = total) +
  geom_histogram(fill = "white", col = "black") +
  scale_fill_grey() +
  xlim(0, 7) +
  theme(legend.position = "none") +
  labs(title = "", x = "Comptage de pouillots véloces", y = "Effectif") +
  theme_classic()

## moyenne et variance

mean(pouillot$total)
var(pouillot$total)
range(pouillot$total)

# comptages plus ou moins surdispersés

ggplot(oiseaux2) +
  aes(x = total) +
  geom_histogram(color = "white") +
  scale_fill_grey() +
  theme(legend.position = "none") +
  facet_wrap( ~ espece, scales = "free") +
  theme_light(base_size = 20) +
  labs(x = "Comptage", y = "Effectif (échelle libre)")

# moyennes et variances par espèces

mean.sp <- round(tapply(oiseaux2$total, INDEX = oiseaux2$espece, "mean"), 2)
var.sp <- round(tapply(oiseaux2$total, INDEX = oiseaux2$espece, "var"), 2)

df.pois <- data.frame(names(mean.sp), mean.sp, var.sp)
colnames(df.pois) = c("espece", "moyenne", "variance")
rownames(df.pois) = NULL
df.pois

# simulation de lois de poisson avec différentes espérances

set.seed(200)
n05 <- rpois(100, 0.5)
mean05 = mean(n05)
n1 <- rpois(100, 1)
mean1 = mean(n1)
n5 <- rpois(100, 5)
mean5 = mean(n5)
n10 <- rpois(100, 10)
mean10 = mean(n10)
npois <- data.frame(n05, n1, n5, n10)

p05 <- ggplot(npois) +
  aes(x = n05) +
  geom_histogram() +
  geom_vline(xintercept = mean05,
             lty = "dashed",
             color = "red") +
  labs(title = paste("\u03BB", "= 0,5"),
       x = "Comptages",
       y = "Effectif") +
  theme_light(base_size = 20)

p1 = ggplot(npois) +
  aes(x = n1) +
  geom_histogram() +
  geom_vline(xintercept = mean1,
             lty = "dashed",
             color = "red") +
  labs(title = paste("\u03BB", "= 1"),
       x = "Comptages",
       y = "Effectif") +
  theme_light(base_size = 20)

p5 = ggplot(npois) +
  aes(x = n5) +
  geom_histogram() +
  geom_vline(xintercept = mean5,
             lty = "dashed",
             color = "red") +
  labs(title = paste("\u03BB", "= 5"),
       x = "Comptages",
       y = "Effectif") +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
  theme_light(base_size = 20)

p10 = ggplot(npois) +
  aes(x = n10) +
  geom_histogram() +
  geom_vline(xintercept = mean10,
             lty = "dashed",
             color = "red") +
  labs(title = paste("\u03BB", "= 10"),
       x = "Comptages",
       y = "Effectif") +
  theme_light(base_size = 20)

plot_grid(p05, p1, p5, p10)

# loi des grands nombres

nb <- NULL

for (i in c(3:208)) {
  k <- 1
  p <- NULL
  repeat {
    if (k == 10)
      break
    xsamp <- sample(1:208, i)
    tp <- pouillot[xsamp, ]
    p <- c(p, tp$total)
    k <- k + 1
  }
  nb <- c(nb, mean(p))
  nb0 <- data.frame(n = 1:length(nb), nb = nb)
}

ggplot(nb0) +
  aes(x = n, y = nb) +
  geom_line() +
  geom_hline(
    yintercept = mean(pouillot$total),
    linetype = "dashed",
    color = "red"
  ) +
  labs(x = "Nombre de points d'écoute", y = "Comptage moyen de pouillots véloces")

#  histogramme de la variation des moyennes entre sous échantillons

set.seed(150)
moy <- NULL
ect <- NULL
k <- 0
tp <- scale(pouillot$total)
repeat {
  sp <- sample(tp, 50)
  moy <- c(moy, mean(sp))
  ect <- c(ect, sd(sp))
  k <- k + 1
  if (k == 40)
    break
}

ggplot(data.frame(moy = moy)) +
  aes(x = moy) +
  geom_histogram(fill = "white", col = "black") +
  ylim(0, 6) +
  labs(x = "Moyenne sur 50 points d'écoute", y = "Effectif") +
  geom_vline(aes(xintercept = mean(moy)),
             lty = "dashed",
             color = "darkred")

## 4.10 L’intervalle de confiance ----------------------------------------------

#  Estimation, estimateur et erreur autour de l’estimation

mean(arbres$hmax)

var(arbres$hmax)
sd(arbres$hmax)
sd(arbres$hmax) / sqrt(length(arbres$hmax))

## L’intervalle de confiance  d’une moyenne

# Histogramme des quantiles d'erreur

n <- 1000
y <- rnorm(n, 0, 1)
q.bas <- quantile(y, p = 0.025)
q.haut <- quantile(y, p = 0.975)
ybas <- data.frame(x = y[which(y < q.bas)])
yhaut <- data.frame(x = y[which(y > q.haut)])

ggplot(data.frame(y = y)) +
  aes(x = y) +
  geom_histogram(fill = "white", color = "black") +
  labs(x = "variable centrée réduite", y = "Fréquence") +
  geom_histogram(
    data = ybas,
    aes(x = x),
    fill = 'darkblue',
    color = 'black'
  ) +
  geom_histogram(
    data = yhaut,
    aes(x = x),
    fill = 'darkblue',
    color = 'black'
  ) +
  scale_x_continuous(breaks = c(-4, -1.96, 0, 1.96, 4)) +
  geom_text(
    x = -2.5,
    y = n * 3 / 100,
    label = "2.5%",
    color = 'darkblue'
  ) +
  geom_text(
    x = 2.5,
    y = n * 3 / 100,
    label = "97.5%",
    color = 'darkblue'
  )

# bornes de l'intervalle de confiance

mean(arbres$hmax) - 1.96 * sd(arbres$hmax) / sqrt(length(arbres$hmax))
mean(arbres$hmax) + 1.96 * sd(arbres$hmax) / sqrt(length(arbres$hmax))

# Fonction pour calculer l'intervalle de confiance de la moyenne

IC <- function(x,
               unit = NULL,
               risk = 95,
               type = "Student") {
  n <- length(x)
  ddl <- n - 1
  m <- mean(x)
  se <- sd(x) / sqrt(n)
  q1 <- (100 - risk) / 2
  q2 <- 100 - q1
  if (type == "Student") {
    ic <- qt(c(q1 / 100, q2 / 100), df = ddl) * se
  }
  if (type == "asymptotic") {
    ic <- qnorm(c(q1 / 100, q2 / 100)) * se
  }
  
  res <- list(m = m, se = se, ic = m + ic)
  
  print(paste("moyenne =", round(m, 2), unit, sep = " "))
  
  if (type == "Student") {
    print(
      paste(
        "Intervalle de confiance ?",
        risk,
        "pour",
        ddl,
        "degrés de liberté =",
        "[",
        round(m + ic[1], 2),
        round(m + ic[2], 2),
        "]",
        unit,
        sep = " "
      )
    )
  }
  
  if (type == "asymptotic") {
    print(paste(
      "Intervalle de confiance ?",
      risk,
      "[",
      round(m + ic[1], 2),
      round(m + ic[2], 2),
      "]",
      unit,
      sep = " "
    ))
  }
  
  return(res)
}

# différents quantiles

iconf90 <- IC(arbres$hmax,
              unit = "mètres",
              risk = 90,
              type = "asymptotic")
iconf95 <- IC(arbres$hmax,
              unit = "mètres",
              risk = 95,
              type = "asymptotic")
iconf99 <- IC(arbres$hmax,
              unit = "mètres",
              risk = 99,
              type = "asymptotic")
iconf999 <- IC(arbres$hmax,
               unit = "mètres",
               risk = 99.9,
               type = "asymptotic")

ggplot(arbres) +
  aes(x = hmax) +
  geom_histogram(fill = "white", col = "black") +
  labs(x = "Hauteur des arbres (m)", y = "Effectif") +
  geom_segment(
    x = iconf90$ic[1],
    xend = iconf90$ic[2],
    y = 5,
    yend = 5,
    color = "steelblue",
    size = 1.5
  ) +
  geom_text(
    x = iconf90$ic[1] - 2.5,
    y = 5,
    label = "intervalle à 90%",
    color = "steelblue"
  ) +
  geom_segment(
    x = iconf95$ic[1],
    xend = iconf95$ic[2],
    y = 10,
    yend = 10,
    color = "slateblue3",
    size = 1.5
  ) +
  geom_text(
    x = iconf95$ic[1] - 2.5,
    y = 10,
    label = "intervalle à 95%",
    color = "slateblue3"
  ) +
  geom_segment(
    x = iconf99$ic[1],
    xend = iconf99$ic[2],
    y = 20,
    yend = 20,
    color = "orange3",
    size = 1.5
  ) +
  geom_text(
    x = iconf99$ic[1] - 2.5,
    y = 20,
    label = "intervalle à 99%",
    color = "orange3"
  ) +
  geom_segment(
    x = iconf999$ic[1],
    xend = iconf999$ic[2],
    y = 30,
    yend = 30,
    color = "darkred",
    size = 1.5
  ) +
  geom_text(
    x = iconf999$ic[1] - 2.5,
    y = 30,
    label = "intervalle à 99.9%",
    color = "darkred"
  ) +
  theme_classic()

## 	L'intervalle de confiance de Student

IC(arbres$hmax,
   unit = "mètres",
   risk = 95,
   type = "Student")

# comparer avec l'intervalle de confiance asymptotique

IC(arbres$hmax,
   unit = "mètres",
   risk = 95,
   type = "asymptotic")

# Effet de l’effectif de l’échantillon sur l’intervalle de confiance de Student

h <- sample(arbres$hmax, 3)
h

intconf <- IC(h, unit = "mètres")
intconf

# simulation des variations de l'intervalle de confiance pour 3 à 100 arbres

m <- rep(NA, 100)
se <- rep(NA, 100)
icl <- rep(NA, 100)
icu <- rep(NA, 100)
for (i in 3:100) {
  h <- sample(arbres$hmax, i)
  intconf <- IC(h)
  m[i] <- intconf$m
  icl[i] <- intconf$ic[1]
  icu[i] <- intconf$ic[2]
  se[i] <- intconf$se
}

df <- data.frame(
  moy = m,
  x = 1:100,
  icl = icl,
  se = se,
  icu = icu
)

ggplot(df) +
  aes(x = x, y = moy) +
  geom_ribbon(aes(ymin = icl, ymax = icu), fill = "grey70") +
  geom_point() +
  labs(x = "répétitions", y = "moyenne ± IC(95%)") +
  theme_classic()

# intervalle de confiance avec une loi normale et une loi de student

icn <- apply(
  df,
  MARGIN = 1,
  FUN = function(x) {
    x[1] + qnorm(c(.025, .975)) * x[4]
  }
)

df <- data.frame(
  moy = m,
  x = 1:100,
  icl = icl,
  se = se,
  icu = icu,
  icnl = icn[1, ],
  icnu = icn[2, ]
)

fig <- ggplot(df) +
  aes(x = x, y = moy) +
  geom_ribbon(aes(ymin = icl, ymax = icu), fill = "grey70") +
  geom_line(aes(x = x, y = icnu), lty = 2, color = "red") +
  geom_line(aes(x = x, y = icnl), lty = 2, color = "red") +
  geom_point() +
  labs(x = "répétitions", y = "moyenne ± IC(95%)") +
  theme_classic()

fig

# zoom:

fig +
  coord_cartesian(xlim = c(0, 30))
