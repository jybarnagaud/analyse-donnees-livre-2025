##----------------------------------------------------------------------------##

# JY. Barnagaud et O. Gimenez. Analyse de données pour les écologues. Biotope
# Editions, 2025.

#                       -------Chapitre 2-------                               #

## commandes R pour la reproduction des exemples du chapitre 
## mis à jour 03 juin 2025

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
pacman::p_load(tidyverse, patchwork, reshape2, ggplot2, viridis, formatR)

# pour installer manuellement chaque package : install.packages(nomdupackage)
# pour charger manuellement chaque package : library(nomdupackage)

## options ---------------------------------------------------------------------

# enlève fond gris, garde grille et grossit la taille des légendes
theme_set(theme_light(base_size = 16)) 

# reproductibilité des simulations
set.seed(2020) 

## jeux de données -------------------------------------------------------------

mesanges <- read.csv2("donnees/mesures_nichoirs.csv",dec=".")
climat <- read.csv2("donnees/climat_LR.csv",dec=".")
pheno <- read.csv2("donnees/phenologie_arbres.csv",dec=".")

## 2.1 Pourquoi une analyse statistique des données ? --------------------------

# relation tarse-poids (mésanges)

p1 <- ggplot(mesanges)+
  aes(x=tarse,y=poids)+
  geom_point(size=3,alpha=0.5,col="darkblue")+
  labs(x="Tarse (mm)",y="Poids (g)",title="(a)")+
  theme_classic(base_size=20)

p2 <- ggplot(mesanges)+
  aes(x=espece,y=poids)+
  geom_boxplot(col="darkblue",fill="darkblue",alpha=0.5)+
  labs(x="Espèce",y="Poids (g)",title="(b)")+
  scale_x_discrete(labels=c("Mésange bleue","Mésange charbonnière"))+
  theme_classic(base_size=20)

# variations mensuelles des températures à Perpignan

climat2 <- subset(climat,POSTE=="Perpignan") 

p3 <- ggplot(climat2)+
  aes(x=mois,y=temperature_moyenne)+
  geom_point(size=3,col="darkblue")+
  theme_classic(base_size=20)+
  scale_x_continuous(breaks=seq(1,12,1))+
  labs(x="Mois",y="Température (°C)",title="(c)")

p4 <- ggplot(climat2)+
  aes(x=factor(mois),y=temperature_moyenne)+
  geom_boxplot(alpha=0.5,col="darkblue",fill="darkblue")+
  theme_classic(base_size=20)+
  labs(x="Mois",y="Température (°C)",title="(d)")

# date de débourrement des pins et chênes

p5 <- ggplot(pheno)+
  aes(x=date_jul,y=bmean,col=essence)+
  geom_point(size=3,alpha=0.5)+
  scale_colour_viridis_d(begin=0,end=0.7,name = "Essence", labels = c("Chêne sessile","Pin sylvestre"))+
  theme_classic(base_size=20)+
  labs(x="Date (référence = 1er janvier)",y="Nombre de bourgeons débourrés",title="(e)")

# affichage et sauvegarde de la figure

(p1+p2)/(p3+p4)/(plot_spacer()+p5)
ggsave(filename="outputs/C2F1.png",width = 15, height = 20,units = "in")

## 2.2 La démarche statistique -------------------------------------------------

# exemple simulé : magnitude et incertitude

x <- seq(0,99,1)
y1 <- 0+0.1*x+rnorm(100,0,1)
y2 <- 0+0.1*x+rnorm(100,0,5)
y3 <- 0+1*x+rnorm(100,0,1)
y4 <- 0+1*x+rnorm(100,0,10)
df.mi <- data.frame(x,y1,y2,y3,y4)

p1 <- ggplot(df.mi)+
  aes(x=x,y=y1)+
  geom_point(size=3,alpha=0.5,col="darkblue")+
  theme_classic(base_size = 20)+
  labs(x="Variable 1",y="Variable 2",title="(a)")+ 
  geom_smooth(method='lm',linetype="dashed",col="darkred")

p2 <- ggplot(df.mi)+
  aes(x=x,y=y2)+
  geom_point(size=3,alpha=0.5,col="darkblue")+
  theme_classic(base_size = 20)+
  labs(x="Variable 1",y="Variable 2",title="(b)")+ 
  geom_smooth(method='lm',linetype="dashed",col="darkred")


p3 <- ggplot(df.mi)+
  aes(x=x,y=y3)+
  geom_point(size=3,alpha=0.5,col="darkblue")+
  theme_classic(base_size = 20)+
  labs(x="Variable 1",y="Variable 2",title="(c)")+ 
  geom_smooth(method='lm',linetype="dashed",col="darkred")


p4 <- ggplot(df.mi)+
  aes(x=x,y=y4)+
  geom_point(size=3,alpha=0.5,col="darkblue")+
  theme_classic(base_size = 20)+
  labs(x="Variable 1",y="Variable 2",title="(d)")+ 
  geom_smooth(method='lm',linetype="dashed",col="darkred")


p.12 <- p1+p2+p3+p4 & xlim(0,100) & ylim(0,70)
p.12

ggsave(filename="outputs/C2F2.png",width = 15, height = 15,units = "in")

## 2.3 L’erreur statistique-----------------------------------------------------

# variations annuelles des températures globales

# données réelles

temperature.globale <- read.csv2("donnees/Global_temperature.csv", dec=".")
str(temperature.globale )

ggplot(data = temperature.globale ) +
  aes(x = Year, y = HadCRUT4) + 
  geom_line(linewidth = 1.5 , col = "steelblue") + 
  labs(x = "Années",
       y = "Anomalies de température (°C)")+
  theme_classic(base_size = 20)

ggsave(filename="outputs/C2F6.png",width = 10, height = 10,units = "in")

# données simulées

z <- data.frame()
for(i in 1:nrow(temperature.globale )){
  y1 <- rnorm(40,temperature.globale [i,2],sd=0.5)
  y2 <- rnorm(40,temperature.globale [i,2],sd=0.05)
  x <- cbind(temperature.globale [i,1],y1,y2)
  z <- rbind(z,x)
}
colnames(z) <- c("annee","anomalie_sim_1","anomalie_sim_2")

g1 <- ggplot(data = z) +
  aes(x = annee, y = anomalie_sim_1) + 
  geom_point(color = 'gray80', size = 2) + 
  ylim(-2.5, 2.5) + 
  stat_summary(geom = 'line', fun = 'mean', color = '#fde725',linewidth = 1) + 
  labs(x = "Années",
       y = "Anomalies de température (°C)") + 
  geom_line(data = temperature.globale , aes(x = Year, y = HadCRUT4), color = '#3b528b',linewidth = 1)+
  theme_classic(base_size = 20)

g2 <- ggplot(data = z) +
  aes(x = annee, y = anomalie_sim_2) + 
  geom_point(color = 'gray80', size = 2) + 
  ylim(-2.5, 2.5) + 
  stat_summary(geom = 'line', fun = 'mean', aes(color = 'moyenne simulée')) + 
  #  stat_summary(geom = 'line', fun = 'mean', colour = 'blue') + 
  labs(x = "Années",
       y = "Anomalies de température (°C)") + 
  geom_line(data = temperature.globale , aes(x = Year, y = HadCRUT4, color = 'moyenne réelle'),linewidth = 1) +
  scale_color_manual(name = "",
                     breaks = c("moyenne simulée", "moyenne réelle"),
                     values = c("#fde725", "#3b528b")) +
  theme_classic(base_size = 20)+
  theme(legend.position = c(0.5,0.2))

g1+g2

ggsave(filename="outputs/C2F7.png",width = 10, height = 5,units = "in")

# exemple avec intervalle de confiance

ggplot(data = temperature.globale ) +
  aes(x = Year, y = HadCRUT4) + 
  geom_ribbon(aes(ymin = ICL, ymax = ICU), fill = "#3b528b",alpha = 0.2) + 
  geom_line(linewidth = 1) + 
  labs(x = "Années",
       y = "Anomalies de température (°C)")+
  theme_classic(base_size = 20)

ggsave('outputs/C2F8.png',width = 10, height = 10,units = "in")

## Erreur, biais, précision : exemple simulé (hauteurs d'arbres)

# vraies hauteurs

arbres <- rnorm(n = 10, 
                mean = 15, 
                sd = 2)
arbres.df <- data.frame(id = 1:10,
                        hauteur = arbres)

# mesures des télémètres

telemetre1 <- sapply(arbres, FUN = function(x){rnorm(5, x-2, 0.5)})
telemetre1 <- melt(telemetre1)  # passe format large à long
colnames(telemetre1) <- c('mesure', 'arbre', 'hauteur')
p1 <- ggplot(data = telemetre1) +
  aes(x = arbre, y = hauteur) + 
  geom_point(color = 'gray40', size = 2) + 
  geom_point(data = arbres.df, aes(x = id, y = hauteur), color = 'red', shape = 3, size = 4) + 
  labs(x = "Arbres",
       y = "Hauteur (m)",
       title = 'télémètre 1')+
  scale_x_continuous(breaks=seq(1,10,1))+
  theme_classic()

telemetre2 <- sapply(arbres, FUN = function(x){rnorm(5, x, 0.5)})
telemetre2 <- melt(telemetre2)  # passe format large à long
colnames(telemetre2) <- c('mesure', 'arbre', 'hauteur')
p2 <- ggplot(data = telemetre2) +
  aes(x = arbre, y = hauteur) + 
  geom_point(color = 'gray40', size = 2) + 
  geom_point(data = arbres.df, aes(x = id, y = hauteur), color = 'red', shape = 3, size = 4) + 
  labs(x = "Arbres",
       y = "Hauteur (m)",
       title = 'télémètre 2')+
  scale_x_continuous(breaks=seq(1,10,1))+
  theme_classic()

telemetre3 <- sapply(arbres, FUN = function(x){rnorm(5, x+3, 4)})
telemetre3 <- melt(telemetre3)  # passe format large à long
colnames(telemetre3) <- c('mesure', 'arbre', 'hauteur')
p3 <- ggplot(data = telemetre3) +
  aes(x = arbre, y = hauteur) + 
  geom_point(color = 'gray40', size = 2) + 
  geom_point(data = arbres.df, aes(x = id, y = hauteur), color = 'red', shape = 3, size = 4) + 
  labs(x = "Arbres",
       y = "Hauteur (m)",
       title = 'télémètre 3')+
  scale_x_continuous(breaks=seq(1,10,1))+
  theme_classic()

telemetre4 <- sapply(arbres, FUN = function(x){rnorm(5, x, 6)})
telemetre4 <- melt(telemetre4)  # passe format large à long
colnames(telemetre4) <- c('mesure', 'arbre', 'hauteur')
p4 <- ggplot(data = telemetre4) +
  aes(x = arbre, y = hauteur) + 
  geom_point(color = 'gray40', size = 2) + 
  geom_point(data = arbres.df, aes(x = id, y = hauteur), color = 'red', shape = 3, size = 4) + 
  labs(x = "Arbres",
       y = "Hauteur (m)",
       title = 'télémètre 4')+
  scale_x_continuous(breaks=seq(1,10,1))+
  theme_classic()

p1 + p2 + p3 + p4

ggsave('outputs/C2F9.png',width = 10, height = 10,units = "in")

