##----------------------------------------------------------------------------##

# JY. Barnagaud et O. Gimenez. Analyse de données pour les écologues. Biotope
# Editions, 2025.

#                       -------Chapitre 6-------                               #

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
  reshape2,
  ggdendro, 
  corrplot,
  GGally,
  
  # multivarié
  ade4, # la librairie que nous utilisons pour les analyses multivariées
  adegraphics, # améliore les graphiques issus d'ade4 (facultatif)
  factoextra, # graphiques issus d'ade4 en version ggplot (facultatif)
  vegan, # nombreuses fonctions utiles pour l'analyse de communautés écologiques
  FactoMineR, # alternative à ade4
  
  # pour les cartes
  sf,
  mapview
  
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

setwd("chapitre_6")

## jeux de données -------------------------------------------------------------

reptiles <- read.csv2("donnees/reptiles_languedoc.csv", header = T, row.names = 1)
clim.lr <- read.csv2("donnees/climat_LR_ACP.csv", header = T, row.names = 1)
amphi <- read.csv2("donnees/amphibiens.edna.csv", header = T, row.names = 1)
petrel.acou0 <- read.csv2("donnees/petrel_acoustique.csv", header = T)
petrel.morpho <- read.csv2("donnees/petrel_morpho.csv", header = T)

## 6.2	les indices de dissimilarité -------------------------------------------

# Exemple fictif : communautés de libellules

libL <- matrix(c(0,1,1,0,1,0,1,1,0,0,0,1), byrow = F, ncol = 3)
colnames(libL) <- c("Libellula depressa", "Calopteryx splendens", "Calopteryx virgo")
rownames(libL) <- c("Site 1", "Site 2", "Site 3", "Site 4")
libL

# distance de Jaccard avec vegan

jaccard.sim <- vegdist(libL, method = "jaccard", binary = T)
jaccard.sim

## 6.3	Représenter graphiquement les similarités entre échantillons------------

# classification sur distance de Jaccard

arbre <- hclust(jaccard.sim)

arbre %>%
  ggdendrogram() + 
  labs(title = "Ressemblance entre communautés de libellules",
       x = NULL,
       y = "Distance de Jaccard")

# Ce qu'on cherche : une classification à deux dimensions

afc.lib <- dudi.coa(libL, nf = 2, scannf = F)

ggplot(afc.lib$li)+
  aes(x = Axis1, y = Axis2) +
  geom_text(size = 4,col="#21918c",label = rownames(afc.lib$li)) +
  xlim(-2,2) +
  ylim(-2,2) +
  labs(x = "Axe de variation 1", y = "Axe de variation 2") +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_text(data = afc.lib$co, aes(x = Comp1, y = Comp2), col = "#440154", size = 5,label = rownames(afc.lib$co)) +
  theme_classic()

# 6.4	Explorer les compositions de communautés : l'AFC--------------------------

# Données sur les reptiles

dim(reptiles)
head(reptiles, n = c(6,5))

# Analyse factorielle des correspondances (AFC) avec ade4 : 

afc.rept <- dudi.coa(reptiles, scannf = F, nf = 3)

# contenu d'un objet de classe *dudi*: 

afc.rept

# accéder aux différents éléments de l'objet : (code non inclus dans l'ouvrage, correspond à la table 6.7). On retrouvera ces mêmes éléments dans tout objet *dudi*, nommés de la même manière.

# valeurs propres
afc.rept$eig

# coordonnées des sites ( = lignes)
head(afc.rept$li)

# scores normés des sites ( = lignes)
head(afc.rept$l1)

# coordonnées des espèces ( = colonnes)
head(afc.rept$co)

# scores normés des espèces ( = colonnes)
head(afc.rept$c1)

# poids des sites ( = lignes)
head(afc.rept$lw)

# poids des espèces ( = colonnes)
head(afc.rept$cw)

## représentations graphiques :
  
# sites et espèces ensemble (non représenté dans l'ouvrage car peu lisible):

scatter(afc.rept,  posieig = "none")

# alternative avec factoextra : 

fc1 <- fviz_ca(afc.rept,col.row = "#3b528b",col.col = "#fde725",
  title = "AFC - reptiles d'Occitanie")
fc1

# sortie factominer (alternative à ade4, non représentée dans l'ouvrage)

afc.rept2 <- CA(reptiles, ncp = 2)

## A partir d'ici, on poursuit avec les sorties "standard" d'ade4 et adegraphics, mais toutes ces sorties existent aussi avec factoextra. A vous de choisir celles que vous préférez.

# projection des sites (version ade4)

s.label(afc.rept$li)

# la même chose version factoextra

fc2 <- fviz_ca_row(afc.rept,col.row = "#3b528b",
                             title = "AFC - reptiles d'Occitanie : sites")
fc2

# deux sites opposés sur la composante 1

stID953 <- sort(colnames(reptiles)[which(reptiles["ID953",]==1)])
stID261 <- sort(colnames(reptiles)[which(reptiles["ID261",]==1)])
intersect(stID953,stID261)

# trois sites proches 

stID971 <- sort(colnames(reptiles)[which(reptiles["ID971",]==1)])
stID978 <- sort(colnames(reptiles)[which(reptiles["ID978",]==1)])
stID979 <- sort(colnames(reptiles)[which(reptiles["ID979",]==1)])
length(intersect(stID971,stID978))/length(union(stID971,stID978))
length(intersect(stID978,stID979))/length(union(stID978,stID979))
length(intersect(stID971,stID979))/length(union(stID971,stID979))

# deux sites opposés sur la composante 2

stID192 <- sort(colnames(reptiles)[which(reptiles["ID192",]==1)])
intersect(stID953, stID192)
stID417 <- sort(colnames(reptiles)[which(reptiles["ID417",]==1)])
intersect(stID979, stID417)

# éboulis des variances (version ade4)

screeplot(afc.rept)

# éboulis des variances (version factoextra)

fc3 <- fviz_screeplot(afc.rept, addlabels = TRUE,
                      title = "éboulis des variances de l'AFC",
                      xlab = "Dimensions (= composantes)",
                      ylab = "% d'inertie expliquée")
fc3

# plan 2,3 (version ade4)

s.label(afc.rept$li, xax = 2, yax = 3)

# plan 2,3 (version factoextra)

fc4 <- fviz_ca_row(afc.rept,axes = c(2,3),
                   col.row = "#3b528b",
                   title = "AFC - reptiles d'Occitanie : sites - plan [2,3]")
fc4

# explication du site extrême

reptiles["ID953",]
sum(reptiles$HEMTUR)
sum(reptiles$TIMLEP)

# zoom (version ade4)

s.label(afc.rept$li, xax = 2, yax = 3, xlim = c(-3,3), ylim = c(-3,3))

# zoom (version factoextra)

fc5 <- fviz_ca_row(afc.rept,axes = c(2,3), xlim = c(-3,3), ylim = c(-3,3),
                   col.row = "#3b528b",
                   title = "AFC - reptiles d'Occitanie : sites - plan [2,3] : zoom")
fc5

# projections des espèces

s.label(afc.rept$co)
s.label(afc.rept$co, xax = 2, yax = 3)
s.label(afc.rept$co, xax = 2, yax = 3, xlim = c(-3,3), ylim = c(-3,3))

# projections des espèces (factoextra)

fc6 <- fviz_ca_col(afc.rept,
                   xlim = c(-3,3),ylim = c(-4,4),
                   col.col = "darkred",
                   title = "(a) AFC - reptiles d'Occitanie : espèces - plan [1,2]")

fc7 <- fviz_ca_col(afc.rept,axes = c(2,3),
                   xlim = c(-2,4),ylim = c(-11,2),
                   col.col = "darkred",
                   title = "(b) plan [2,3]")

fc8 <- fviz_ca_col(afc.rept,axes = c(2,3), xlim = c(-3,3), ylim = c(-3,3),
                   col.col = "darkred",
                   title = "(c) plan [2,3] : zoom")

(fc6 + fc7) / (fc8 + plot_spacer())

## 6.5	L’analyse en composantes principales (ACP) -----------------------------

summary(clim.lr)
head(clim.lr)

# centrer-réduire les variables

clim.lr.sc <- scale(clim.lr)

# matrice des corrélations

cor.clim <- cor(scale(clim.lr))
cp1 <- corrplot(cor.clim)
recordPlot(cp1)

## L’ACP sous R

# Pour coder une ACP. Par défaut, R affiche une fenêtre et demande de sélectionner le nombre de composantes voulue, puis Entrée:  

pc <- dudi.pca(clim.lr) # sélectionner 2 composantes + entrée
2

# Vous pouvez aussi sélectionner a priori le nombre de composantes à conserver :

pc <- dudi.pca(clim.lr, scannf = F, nf = 2) 

## Choisir le nombre de composantes principales

# éboulis des variances, manuellement : 

barplot(pc$eig, 
        col = "gray70", 
        xlab = "Composantes principales", 
        ylab = "Valeurs propres", 
        names.arg = c(1:9))

# avec ade4 : 

screeplot(pc)

# avec factoextra : 
  
pc.scr <- fviz_screeplot(pc,
                         title = "Eboulis des valeurs propres de l'ACP",
                         ylab = "% Inertie expliquée")
pc.scr

# éboulis des variances - choix en fonction de la cassure (du coude)

bp <- barplot(pc$eig,
              col = "gray70", 
              xlab = "Composantes principales", 
              ylab = "Valeurs propres", 
              names.arg = c(1:9),
              cex.lab = 0.7,cex.axis= 0.7,cex.names = 0.7)
segments(x0 = bp[1],x1 = bp[3],y0 = pc$eig[1],y1 = pc$eig[3],col="red",lwd=2)
segments(x0 = bp[3],x1 = bp[9],y0 = pc$eig[3],y1 = pc$eig[9],col="red",lwd=2)
segments(x0 = bp[3],x1 = bp[3],y0 = pc$eig[3],y1 = 4,col="red",lty= "dotted",lwd = 2)
text(bp[3],4.2,"cassure",col="red",cex = 0.8)

# éboulis des variances - choix d'un % de variance expliquée

in.expl <- 100*pc$eig / sum(pc$eig)
cum.in.expl <- 100*cumsum(pc$eig) / sum(pc$eig)

bp <- barplot(pc$eig,
        col = "gray70", 
        xlab = "Composantes principales", 
        ylab = "Valeurs propres", 
        names.arg = c(1:9),
        cex.lab = 0.7,cex.axis= 0.7,cex.names = 0.7)

text(bp[1],y = 2.5,"% inertie expliquée de chaque valeur propre",cex = 0.5,font = 2,col = "darkred",pos = 4)
text(bp, y = 2, round(in.expl,2),cex = 0.5,col = "darkred")

text(bp[1],y = 1.5,"% d'inertie expliquée cumulé",cex = 0.5,font = 2,col = "darkblue",pos = 4)
text(bp, y = 1, round(cum.in.expl,2), font = 3,cex = 0.5,col = "darkblue")

# quelques exemples fictifs d'éboulis de variance plus ou moins faciles (modifié dans Powerpoint pour ajout des textes d'interprétation)

x3 <- c(0.40,0.35,0.15,0.05,0.025,0.015,0.005,0.004,0.001)
x1 <- c(0.90,0.08,0.009,0.008,0.002,0.0005,0.00025,0.00025)
x6 <- c(0.30,0.25,0.20,0.10,0.09,0.05,0.005,0.004,0.001)
x2 <- c(0.45,0.30,0.15,0.05,0.025,0.0125,0.0124,0.0001)
x4 <- c(0.30,0.25,0.24,0.16,0.05)
x5 <- c(0.30,0.25,0.20,0.15,0.10)

df <- data.frame(y = c(x3, x1, x6, x2, x4, x5),
                 var = factor(c(rep("a: cassure nette", length(x3)),
                         rep("b: une composante domine", length(x1)),
                         rep("c: deux cassures", length(x6)),
                         rep("d: décroissance graduelle (1)", length(x2)),
                         rep("e: cassure tardive", length(x4)),
                         rep("f: décroissance graduelle (2)", length(x5)))),
                 x = c(1:length(x3),
                       1:length(x1),
                       1:length(x6),
                       1:length(x2),
                       1:length(x4),
                       1:length(x5)))

sim.varex <-  ggplot(df) + 
  aes(x = x, y = y) + 
  geom_col() +
  labs(x = "Composantes principales", y = "Valeurs propres") + 
  facet_wrap(vars(var))+
  theme_light(base_size = 16)

sim.varex

## Interpréter les composantes principales : le cercle des corrélations

# Cercle des corrélations (version ade4)

s.corcircle(pc$co)

# Cercle des corrélations (version factoextra)

pc.corcircle <- fviz_pca_var(pc,
                             title = "Variables - ACP sur variables climatiques")
pc.corcircle

# inertie

inertia.dudi(pc, col.inertia = TRUE)

# Projection des individus ( = sites, = lignes), avec ade4

s.label(pc$li)

# Projection des individus ( = sites, = lignes), avec factoextra

pc.ind <- fviz_pca_ind(pc,
                             title = "Sites - ACP sur variables climatiques")
pc.ind

# Coordonnées des individus : 

coords <- pc$li
head(coords)

## Cartographie de la grille (librairies sf et mapview) :

# ajoute les identifiants des cellules de la grille

coords <- pc$li %>%
  mutate(id = as_factor(rownames(.))) %>% 
  mutate(id = str_remove(id, 'ID')) # supprime 'ID' dans le nom des cellules

# couche spatiale vectorielle contenant la grille – jointure avec l’ACP

grille <- st_read("donnees/grid_extent.geojson", crs = 2154) %>% # lit le shapefile la grille
  inner_join(coords, by = c('nms_s3_' = 'id')) # jointure avec les données de l'ACP

mapviewOptions(fgb = FALSE)

carte_grille <-
  grille %>%
  st_boundary() %>%
  mapview(
    color = 'red',
    lwd = 2,
    legend = FALSE,
    map.types = 'OpenStreetMap'
  )
carte_grille # la carte interactive

mapshot(x = carte_grille, file = "outputs/C6F5.png") # sauvegarde une image

# carte du premier axe de l’ACP

carte_pc1 <-
  grille %>%
  mapview(
    zcol = "Axis1",
    legend = TRUE,
    color = 'black',
    lwd = 2,
    map.types = 'OpenStreetMap',
    at = c(-9.7, -6.4, -3.1, -1.6, -0.2, 1.1, 2.3, 3.5),
    layer.name = 'ACP - Composante 1'
  )
carte_pc1
mapshot(x = carte_pc1, 
        file = "outputs/C6F22a.png",
        remove_controls = c("zoomControl", "layersControl", "homeButton"))


# carte du deuxième axe de l’ACP

carte_pc2 <-
  grille %>%
  mapview(
    zcol = "Axis2",
    legend = TRUE,
    color = 'black',
    lwd = 2,
    map.types = 'OpenStreetMap',
    at = c(-2.6, -2, -1.3, -0.6, 0.1, 0.7, 1.2, 2.3),
    layer.name = 'ACP - Composante 2'
  )
carte_pc2 # la carte interactive


## 6.6 L’analyse en coordonnées principales (PCOA)------------------------------

## Les données : analyser un échantillonnage de communautés par ADN environnemental

#données
                                                                    
head(amphi)
summary(amphi)

#graphique variables x variables (on n'en garde que quelques-unes pour la lisibilité)

pairs.amphi <- ggpairs(amphi, columns = c("HYLMER", "LISHEL", "PELCUL", "PELPUN")) +
  theme_light(base_size = 16)

pairs.amphi

# matrice de dissimilarité (avec la librairie vegan, déjà utilisée plus haut)

amphi.dist <- vegdist(amphi, method = "bray")

# ordination (voir l'éboulis des valeurs propres négatives -  à faire tourner dans la console)

amphi.pco <- dudi.pco(amphi.dist)
10

# ordination (on garde toutes les composantes)

amphi.pco <- dudi.pco(amphi.dist, scannf = F)

# transformation en racine carrée des distances afin d'éviter les valeurs propres négatives

amphi.pco2 <- dudi.pco(sqrt(amphi.dist), scannf = F)

# éboulis des valeurs propres par défaut: (à faire tourner dans la console)

amphi.pco2 <- dudi.pco(sqrt(amphi.dist), scannf = F, nf = 4)

# projection des individus

s.label(amphi.pco2$li,
        plabels=list(cex = 0.8,boxes=list(draw = F)),ppoints = list(cex = 0))

# projection des variables (= espèces) : on utilise pour ça l'artifice de les projeter comme des variables supplémentaires, c'est à dire des variables que l'on superpose a posteriori sur le plan sans qu'elles influent sur la construction des composantes principales. Pour cela, ade4 positionne simplement chaque espèce au barycentre des sites où elle a été observée. C'est évidemment une "tricherie" car ces variables ont bien participé à la construction des composantes via la matrice de dissimilarité. 
                                                                    
sp <- supcol(amphi.pco2, amphi)
s.arrow(sp$cosup,
        plabels =
          list(cex = 0.8, boxes = list(draw = F)),
        ppoints = list(cex = 0))
                                                                 
## 6.7	L’analyse multidimensionnelle non métrique (NMDS)-----------------------
                                                                       
## La NMDS sous R

nmds <- metaMDS(amphi, distance = "bray", try = 40)
sppscores(nmds) <- amphi

# projection des sites et des espèces ensemble:

ordiplot(nmds, type = "text")

# qualité de représentation

stressplot(nmds)

## 6.8	L’analyse discriminante-------------------------------------------------
                                                                       
# données

petrel.acou0$Species <- factor(petrel.acou0$Species)
head(petrel.acou0)
                                                                    
# agrégation: moyenne des variables acoustiques par individu
                                                            
petrel.acou <- aggregate(
  petrel.acou0[, -c(1:3)],
  by = list(petrel.acou0$Indiv, petrel.acou0$Species),
  FUN =
    "mean"
)

head(petrel.acou)

colnames(petrel.acou)[1:2] <- c("Indiv", "Species")
rownames(petrel.acou) <- petrel.acou$Indiv

## Première étape : construire une ordination sur notre table de variables
                                                                       
# ACP (on garde toutes les composantes)

pc.petrel <- dudi.pca(petrel.acou[, -c(1:2)], scannf =
                        F, nf = ncol(petrel.acou[, -c(1, 2)]))


# inertie:

pc.petrel.scr <- fviz_screeplot(pc.petrel, title = "Eboulis des variances - variables acoustiques")
pc.petrel.scr

# On détache le package adegraphics pour éviter des conflits de version au cours des analyses qui suivent
                                                        
detach("package:adegraphics", unload =
         TRUE)

# FIGURE 6.30 : cercles de corrélations

s1 <- fviz_pca_var(pc.petrel, axes = c(1, 2), title = "") +
  labs(subtitle =
         "(a)")
s2 <- fviz_pca_var(pc.petrel, axes =
                     c(1, 3), title = "") +
  labs(subtitle =
         "(b)")
s3 <- fviz_pca_var(pc.petrel, axes =
                     c(2, 3), title = "") +
  labs(subtitle =
         "(c)")
(s1 + s2) / (s3 + plot_spacer())

# différenciation des groupes
                                                                
s.class(
  pc.petrel$li,
  fac = factor(petrel.acou$Species),
  label = c("Pétrel bleu", "Prion de la Désolation"),
  col = c("darkblue", "darkred")
)

# analyse discriminante

discr.petrel <- discrimin(pc.petrel,
                          factor(petrel.acou$Species),
                          scannf = F,
                          nf = 2) # 2 groupes donc 1 seul axe discriminant

# contribution des variables

sum(discr.petrel$eig) / sum(pc.petrel$eig)

# test par permutations
                                                                   
test <- randtest(discr.petrel)
plot(test,
     main = "Test par permutations sur l'analyse discriminante",
     ylab = "Fréquence",
     xlab = "sim (ratio variance inter-groupes / variance totale)")

test

## Etudier l’impact du facteur de groupement sur les variables de l’ordination
# ! ne fonctionne que sur les versions d'ade4 postérieures à 2021
# conflit avec adegraphics et/ou factoextra, détacher ces librairies au préalable (voir plus haut)

plot(discr.petrel)

## 6.9	L’analyse de co-inertie-------------------------------------------------

# données

rownames(petrel.morpho) <- petrel.morpho$Indiv
petrel.morpho <- petrel.morpho[, -1]
petrel.morpho <- petrel.morpho[rownames(petrel.acou), ] # réordonner les lignes comme la matrice acoustique
head(petrel.morpho)

## ACP avec la pondération par les lignes de l'ACP sur les variables acoustiques (pc.petrel$lw)
                          
pc.morpho <- dudi.pca(petrel.morpho,
                      row.w = pc.petrel$lw,
                      scannf = F,
                      nf = 2)

# éboulis des variances

scr.pc.morpho <- fviz_screeplot(pc.morpho, title = "Eboulis des variances - variables morphologiques")
scr.pc.morpho

# cercle des corrélations

cor.circ.morpho <- fviz_pca_var(pc.morpho)
cor.circ.morpho

# projection des individus, avec ellipses par espèces (on récupère les espèces dans la matrice de données des variables acoustiques, **attention à bien ordonner les deux matrices exactement de la même manière auparavant, c'est de toute façon nécessaire pour la suite car ade4 ne le fera pas automatiquement**)

s.class(
  pc.morpho$li,
  fac = factor(petrel.acou$Species),
  label = c("Pétrel bleu", "Prion de la Désolation"),
  col = c("darkblue", "darkred")
)

## L’analyse de co-inertie sous R

# analyse de coinertie

coinertie.petrels <- coinertia(pc.petrel, pc.morpho,scannf=F,nf=2)

# test par permutations

rd.coi <- randtest(coinertie.petrels, fixed = 2)

plot(rd.coi,
     main = "Test par permutations sur l'analyse de coinertie",
     xlab = "sim",
     ylab = "Fréquence")

rd.coi

# ensemble de graphiques de synthèse de la coinertie

plot(coinertie.petrels)

# les différentes parties de la sortie graphique

par(mfrow = c(2, 1))
s.corcircle(coinertie.petrels$aX, sub = "ACP acoustique")
s.corcircle(coinertie.petrels$aY, sub = "ACP morphologique")

# projections des lignes

s.arrow(coinertie.petrels$l1, clabel = 0.7)

# projection des colonnes

s.arrow(coinertie.petrels$c1, clabel = 0.7)

# projection des variables

s.class(
  coinertie.petrels$lX,
  fac = factor(petrel.acou$Species),
  label = c("Pétrel bleu", "Prion de la Désolation"),
  col = c("darkblue", "darkred")
)

## Encart 6.1. mesures répétées--------------------------------------------------

# l'ACP sur les données répétées par individus

pc.petrels.rep <- dudi.pca(petrel.acou0[,-c(1:3)],scannf=F,nf=2)  

# la BCA

bca.petrels <- bca(pc.petrels.rep,fac = factor(petrel.acou0$Indiv),scannf = F,nf = 2) 

# l'exploration de la BCA se fait comme pour l'ACP : 

par(mfrow=c(2,1))
s.corcircle(bca.petrels$co,clabel=0.8)
s.label(bca.petrels$li,clabel=0.8) 

# Remarquez que bien qu'il y ait dans la table de données une ligne par enregistrement, la projection représente bien les individus biologiques (oiseaux) - mais leur position tient compte des mesures répétées.
                                                                                                                             
# Ellipses pour les deux espèces:
  
s.class(
  bca.petrels$li,
  fac = factor(petrel.acou$Species),
  label = c("Pétrel bleu", "Prion de la Désolation"),
  col = c("darkblue", "darkred")
)

# l'analyse discriminante sur la BCA

sp.bca <- unique(petrel.acou0[, c("Indiv", "Species")])
discrim.bca <- discrimin(
  bca.petrels,
  fac = factor(sp.bca$Species),
  scannf = F,
  nf = 1
)

# test par permutations sur cette analyse discriminante :

test.discrim.bca <- randtest(discrim.bca)

plot(test.discrim.bca,
     main = "Test par permutations sur l'analyse discriminante avec la BCA",
     ylab = "Fréquence",
     xlab = "sim (ratio variance inter - groupes / variance totale)")

test.discrim.bca <- randtest(discrim.bca)
plot(test.discrim.bca,
     main = "Test par permutations sur l'analyse discriminante avec la BCA",
     ylab = "Fréquence",
     xlab = "sim (ratio variance inter-groupes / variance totale)")

## la coinertie sur la BCA - cette partie n'est pas dans l'ouvrage

# il faut refaire l'ACP sur les données morphométriques pour la pondérer en fonction de la BCA

pc.morpho.wbca <- dudi.pca(petrel.morpho,
                           row.w = bca.petrels$lw,
                           scannf = F,
                           nf = 2)

# la coinertie :

coi.bca <- coinertia(bca.petrels, pc.morpho.wbca, scannf = F, nf = 2)

# la projection des variables a un peu changé par rapport à la coinertie sur données acoustiques moyennées - mais sans beaucoup affecter l'interprétation

s.arrow(coi.bca$l1, clabel = 0.7)
s.arrow(coi.bca$c1, clabel = 0.7)

# de la même manière, la projection des individus reste à peu près la même. C'est assez logique vu la différence entre ces deux espèces, mais l'impact des données répétées serait bien plus fort sur des espèces aux chants plus proches (ou si de manière générale la variance associée aux données répétées était forte par rapport à la variance inter-groupes)

s.class(coi.bca$lX, fac = factor(petrel.acou$Species))

## 6.10	l’analyse canonique des correspondances---------------------------------

## La CCA sous R

# On reprend l'AFC sur les reptiles (afc.rept) du début de ce script, et deux des variables climatiques :

cca.rept <- pcaiv(afc.rept, clim.lr[, c("altitude_mediane", "saison_precip")], scannf = F)

# inertie de la CCA

summary(cca.rept)

# projection

s.label(cca.rept$c1,
        xlim = c(-1, 1),
        boxes = FALSE,
        clabel = 0.7)
s.arrow(cca.rept$cor, add.plot = T, clab = 1)
