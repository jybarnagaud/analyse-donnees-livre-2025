#--------------------------------#
#### Quelques calculs simples ####
#--------------------------------#

# Ce script permet d'effectuer quelques calculs simples à visée démonstrative
# Auteur: JY Barnagaud (08/11/2021)

#------------------------------#
#### Créer quelques valeurs ####
#------------------------------#

objet1<-8
objet2<-6
objet3<-objet2-objet1 # ici, je calcule la différence entre mes deux objets de départ

# objet4=objet3/0 # cette commande renvoie un message d'erreur, je la bloque mais la garde pour mémoire

#------------------------------#
#### Créer d'autres valeurs ####
#------------------------------#

objet4<-exp(3/(objet1+objet2+objet3)) # attention à la position des parenthèses
objet5<-log(objet4)