##----------------------------------------------------------------------------##

# JY. Barnagaud et O. Gimenez. Analyse de données pour les écologues. Biotope
# Editions, 2025.

#                       -------Chapitre 3-------                               #

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
pacman::p_load(formatR, ggplot2)

# pour installer manuellement chaque package : install.packages(nomdupackage)
# pour charger manuellement chaque package : library(nomdupackage)

## options ---------------------------------------------------------------------

# enlève fond gris, garde grille et grossit la taille des légendes
theme_set(theme_light(base_size = 16))

# reproductibilité des simulations
set.seed(2020)

## 3.3 Une introduction au logiciel R ------------------------------------------

2 + 2

## 3.4 le principe de R --------------------------------------------------------

# Les objets

objet1 <- 2 + 2
objet1

objet2 <- 5 + 5
objet2

objet3 <- objet1 + objet2
objet3

# Les fonctions

exp.objet2 <- exp(objet2)
exp.objet2

base1 <- log(2)
base10 <- log(2, base = 10)

base1
base10

# La librairie ou paquet de fonctions (le “package”)

install.packages("ggplot2", repos = 'http://cran.us.r-project.org') # installer un package

library(ggplot2) # charger un package

# Créer ses objets et ses fonctions

mafonction <- function(x) {
  1 / log(x)
}

test <- mafonction(10)
test

## 3.5 L'aide de R -------------------------------------------------------------

# Accès à l'aide :

help(log)
? log

# L'exemple de l'aide pour la fonction log :

log(exp(3))
log10(1e7) # = 7

x <- 10 ^ -(1 + 2 * 1:9)
cbind(x, log(1 + x), log1p(x), exp(x) - 1, expm1(x))

## 3.6 Le flux de travail ------------------------------------------------------

# l'exemple de script du livre

#--------------------------------#
#### Quelques calculs simples ####
#--------------------------------#

# Ce script permet d'effectuer quelques calculs simples à visée démonstrative
# Auteur: JY Barnagaud (08/11/2021 - mise à jour 19/11/2022)

#------------------------------#
#### Créer quelques objets  ####
#------------------------------#

objet1 <- 8
objet2 <- 6
objet3 <- objet2 - objet1

# objet4=objet3/0 # commande renvoyant un message d'erreur, je la bloque mais la garde pour mémoire

objet4 <- exp(3 / (objet1 + objet2 + objet3)) # attention à la position des parenthèses

objet5 <- log(objet4)

## 3.7 Les données sous R ------------------------------------------------------

## Les formats d'objets

# un objet de mode numeric

x <- 123
mode(x)

# un objet de mode character

z <- ("A")
mode(z)

# un objet de mode logical

y <- (x == 0)
mode(y)

# au fait, quel est le mode de la fonction mode?

mode(mode)

## Créer ses propres objets sous R

# Compiler des données manuellement sous R :

terrier <- c("T1",
             "T2",
             "T3",
             "T4",
             "T5",
             "T6",
             "T7",
             "T8",
             "T9",
             "T10",
             "T11",
             "T12",
             "T13")
terrier

occupant.s1 <- c(
  "Petrel_bleu",
  "Petrel_bleu",
  "Petrel_bleu",
  "Prion_antarctique",
  "Prion_antarctique",
  "Prion_antarctique",
  "Prion_antarctique"
  ,
  "Prion_antarctique",
  "vide",
  "vide",
  "vide",
  "vide",
  "vide"
)

occupant.s2 <- c(
  "Petrel_bleu",
  "Petrel_bleu",
  "vide",
  "Prion_antarctique",
  "Prion_antarctique",
  "Prion_antarctique",
  "Prion_antarctique",
  "vide",
  "Petrel_bleu"
  ,
  "Petrel_bleu",
  "Prion_antarctique",
  "vide",
  "vide"
)

class(terrier)
class(occupant.s1)
class(occupant.s2)

# passer une chaîne de caractère en facteur

occupant.s1 <- factor(occupant.s1)
occupant.s2 <- factor(occupant.s2)

class(occupant.s1)
class(occupant.s2)

occupant.s1
mode(occupant.s1)

# Table de données ("data frame"):

petrels <- data.frame(terrier, occupant.s1, occupant.s2)

# Pour le réécrire sous forme de matrice :

? matrix

mat.petrels <- matrix(c(2, 0, 2, 0, 4, 1, 1, 1, 2), nrow = 3, byrow = T)
mat.petrels <- matrix(c(2, 0, 1, 0, 4, 1, 2, 1, 2), nrow = 3, byrow = F)

## Le système d'indexation

# un grand vecteur:

x <- seq(1, 100, length.out = 500)

# la 100ème valeur:

x[100]

# les valeurs de x, de la première à la vingtième:

x[1:20]

# la cinquième, la dixième et la cinquantième valeurs:

x[c(5, 10, 50)]
x[1:10]

# est équivalent à

x[-c(11:500)]

# l'occupant du septième terrier à la deuxième saison:

petrels[7, 2]

# et les occupants des terriers 3, 9, 12, aux deux saisons:

petrels[c(3, 9, 12), ]

# juste les occupants des terriers à la saison 1:

petrels[, c(1, 2)]

# ou (équivalent):

petrels[, -3]

#  tous les terriers, sauf les 4 derniers:

petrels[-c(10:13), ]

# seulement les identités des occupants lors de la première saision:

petrels[, "occupant.s1"]

# pour un objet de classe data.frame, vous disposez d'un raccourci:

petrels$occupant.s1

# noms de colonnes

rownames(mat.petrels) <- c("s1_Petrel_bleu", "s1_Prion_antarctique", "s1_vide")
colnames(mat.petrels) <- c("s2_Petrel_bleu", "s2_Prion_antarctique", "s2_vide")
mat.petrels

# a:

mat.petrels["s1_Prion_antarctique", ]

#b:

petrels[c(1:4, 10:13), "occupant.s1"]

# c:

mat.petrels["s1_vide", "s2_vide"]

# d:

colnames(petrels) <- c("nom_terrier", "occupant.saison1", "occupant.saison2")
rownames(petrels) <- petrels[, "nom_terrier"]
petrels2 = petrels[, -1] # n'oubliez pas l'assignation si vous ne voulez pas perdre ce jeu de données reformaté!

#e:

petrels2[c("T7", "T8", "T9"), "occupant.saison1"]

## 3.8 Importer des données sous R ---------------------------------------------

# Le répertoire de travail
# NB : attention au sens des "/", inverse par rapport au défaut de Windows. Attention également à encadrer le chemin d'accès de guillemets (à la française : "" ou à l'anglaise : '')

# Chemin d'accès vers le dossier contenant les données du projet (ajustez en fonction de votre propre arborescence):

setwd('D:/analyse_donnees/chapitre_3/donnees')

# Chemin d'accès vers la racine du projet :

setwd('D:/analyse_donnees/chapitre_3')

# Voir le répertoire de travail choisi :

getwd()

# Liste des fichiers contenus dans le répertoire de travail :

list.files()

# Dans le sous-répertoire "donnees":

list.files('donnees')

# Importer un tableau de données au format .txt

# Si on a spécifié le sous-dossier "donnees" comme répertoire de travail :

arbres <- read.table("hauteurs_arbres.txt", header = T, sep = "\t")

# Si on a spécifié la racine du projet comme répertoire de travail (c'est le cas pour les projets R de ce livre) :

arbres <- read.table("donnees/hauteurs_arbres.txt",
                     header = T,
                     sep = "\t")

## Importer un tableau de données au format .csv

arbres <- read.csv2("donnees/hauteurs_arbres.csv", header = T)

## Sauvegarder des données depuis R

# En direction d'un sous-dossier "outputs" qui aura été créé au préalable dans le répertoire de travail:

write.table(
  arbres,
  "outputs/donnees_arbres.txt",
  append = F,
  sep = "\t",
  row.names = F
)

# Pour des fichiers csv:

write.csv(arbres, "outputs/donnees_arbres.csv") # standard anglo-saxon
write.csv2(arbres, "outputs/donnees_arbres.csv") # standard français

# NB : il existe des fonctions pour sauvegarder dans de nombreux formats, y compris
# Excel, Access, fichiers de forme spatiaux, etc.

## Sauvegarder des objets ou ensembles d'objets R

# Tous les objets de la session en cours :

save.image("outputs/masession.RData")

# Juste quelques objets:

save(arbres, petrels, mat.petrels, file = "outputs/donnees_chapitre3.RData")

# En format RDS (ne sauve que le contenu d'un objet, pratique pour les objets R lourds ou complexes):

saveRDS(arbres, file = "outputs/arbres.rds")

# Pour ouvrir un .RData:

load("outputs/masession.RData")

## 3.9 Explorer un tableau de données ------------------------------------------

# Classe et dimensions de l'objet

class(arbres)
dim(arbres)

# pour afficher le nombre de lignes seul

nrow(arbres)

# et les colonnes

ncol(arbres)

# Synthèses d’un objet

View(arbres)
head(arbres)
summary(arbres)
str(arbres)

# synthèse d'un facteur

arbres$cercle <- factor(arbres$cercle)
summary(arbres$cercle)

# réorganisation des niveaux d'un facteur

arbres$cercle <- factor(arbres$cercle, levels = c("N", "E", "S", "W"))

# synthèse d'un objet character

arbres$nichoir <- as.character(arbres$nichoir)
arbres$id <- as.character(arbres$id)
summary(arbres$nichoir)
summary(arbres$id)

## Les filtres de données

subE <- subset(arbres, cercle == "E")

# Réponses à l'exercice :

p1 <- subset(petrels, occupant.s1 == "Petrel_bleu")
p2 <- subset(petrels, occupant.s1 == "vide" &
               occupant.s2 == "Petrel_bleu")
p3 <- subset(petrels, occupant.s2 != "Prion_antarctique")
p4 <- subset(petrels,
             occupant.s1 != "Prion_antarctique" |
               occupant.s2 != "Prion_antarctique")
a1 <- subset(arbres, nichoir %in% c(141, 1441) &
               cercle %in% c("N", "W"))
a2 <- subset(arbres, hmax >= 24)


# Conserver les données correspondant au cercle 'E' (Est) :

arbres[arbres$cercle == "E", ]

# Autre solution :

arbres$cercle == "E"

# Deux conditions :

petrels[petrels$occupant.s1 == "vide" &
          petrels$occupant.s2 == "Petrel_bleu", ]

## 3.10 Synthèses et calculs sur des tableaux de données -----------------------

mesanges <- read.csv2("donnees/mesures_nichoirs.csv", dec = ".")
summary(mesanges)

mesanges$nichoir <- as.character(mesanges$nichoir)

## Moyenne d’une variable ventilée par catégories

# En filtrant des données (ici par espèce):

parmaj <- subset(mesanges, espece == "PARMAJ")

mean(parmaj$poids)

parcae <- subset(mesanges, espece == "PARCAE")

mean(parcae$poids)

# La même chose peut être obtenue automatiquement avec tapply:

tapply(mesanges$poids, INDEX = mesanges$espece, FUN = "mean")
tapply(mesanges$tarse, INDEX = mesanges$nichoir, FUN = "mean")

# Même chose avec aggregate, qui peut appliquer la même fonction sur plusieurs variables d'un coup:

moy.nichoir <- aggregate(mesanges[, c("tarse", "poids")], by = list(mesanges$nichoir), FUN =
                           "mean")
head(moy.nichoir)

## Changer le format d'un tableau de données

# D'un format "long" (une observation = une ligne de données) à un format "large" (observations réparties en lignes et colonnes):

arbres.large <- reshape(
  arbres[, c("nichoir", "cercle", "hmax")],
  v.names = "hmax",
  idvar = "nichoir",
  timevar = "cercle",
  direction = "wide"
)
head(arbres.large)

# D'un format "large" à un format "long":

arbres.long <- reshape(
  arbres.large,
  idvar = "nichoir",
  varying = c("hmax.E", "hmax.N", "hmax.S", "hmax.W"),
  times = c("E", "N", "S", "W"),
  v.names = "hmax",
  direction = "long"
)
head(arbres.long)

## Opérations par lignes et par colonnes

moy.arbres <- tapply(arbres$hmax, INDEX = arbres$nichoir, FUN = "mean")

# La même opération sur toutes les lignes d'une table:

arbres.large$hmax.moy <- apply(arbres.large[, -1], MARGIN = 1, FUN = "mean")
head(arbres.large)

# La même opération sur toutes les colonnes d'une table:

apply(arbres.large[, -1], MARGIN = 2, FUN = "mean")

## Jointures de tableaux de données

arbres.mesanges <- merge(mesanges, arbres.large, by = "nichoir")

## Exercice :

dendro <- read.csv2("donnees/dendrometrie.csv", dec = ".")
summary(dendro)

dendro <- read.csv2("donnees/dendrometrie_corrige.csv", dec = ".")
summary(dendro)

dendro$releve <- factor(dendro$releve)
dendro$ST <- apply(dendro[, c("STF", "STR")], MARGIN = 1, FUN = "sum")

dendro.long <- reshape(
  dendro,
  idvar = "releve",
  varying = c("STCH", "STF", "STR", "PS", "ST"),
  times = c("STCH", "STF", "STR", "PS", "ST"),
  v.names = "surface",
  direction = "long"
)

dendro.arbres <- merge(dendro, arbres.large, by.x = "releve", by.y = "nichoir")

tapply(dendro.arbres$hmax.moy,
       INDEX = dendro.arbres$parcelle,
       FUN = "mean")

resume <- aggregate(dendro.arbres[, c("STCH", "STF", "STR", "PS", "ST")],
                    by = list(dendro.arbres$parcelle),
                    FUN = "mean")
colnames(resume)[1] = "nichoir"
resume

## 3.11 Résolution de problèmes ------------------------------------------------

# Les messages d'erreur

seq(1, a)
arbres[, 10]

# Les messages d'avis
log(-1)

# Exercice

x1 <- matrix(1:100, ncol = 20)
x1[21, ]

x2 <- rep(c("A", "B", "C"), 20)
sum(x2)

ecartype(x1)

sd(z)

x3 <- c(1:15)
x4 <- cbind(x2, x3)

x6 <- matrix(1:20, ncol = 2)
x7 <- matrix(1:9, ncol = 3)
x8 <- rbind(x6, x7)

which(x6 = 7)

x9 <- rnorm(100, 0, 1)
x10 <- runif(100, 1, 2)
mean(x9)
mean(x9, weights = x10)

x12 <- factor("A", "A", "A", "B", "B", "C")

