# analyse-donnees-livre-2025
données et scripts pour la version 2025 du livre pour Biotope sur l'analyse de données

************************ Archive en ligne de l'ouvrage ****************************

Jean-Yves Barnagaud et Olivier Gimenez, 2025. XXX. Biotope Editions.

jean-yves.barnagaud@ephe.psl.eu

Cette archive contient les données et scripts R permettant la réplication des exemples et figures du livre.

Dernière mise à jour : 31/01/2025

************************************************************************************

*** Avertissements ***

- Les scripts fournis ont été vérifiés et fonctionnent correctement à la date de mise en vente de l'ouvrage. 
Le logiciel R et les paquets de fonctions (packages) utilisés étant fréquemment mis à jour, nous ne pouvons
garantir une fonctionnalité pérenne. Dans la très grande majorité des cas, les erreurs pouvant advenir suite
à une mise à jour sont mineures et faciles à résoudre par une recherche rapide sur internet. 

- Les jeux de données sont mis à votre disposition dans le but exclusif de répliquer les analyses décrites 
dans le livre. Ils sont fournis à titre gratuit pour cet unique usage par leurs propriétaires, personnes privées, 
associations, bureaux d'études et agences d'Etat. Les données n'ont pas été modifiées, mais des informations clé
ont été supprimées pour éviter tout usage en dehors du contexte de ce livre. Les jeux de données, leurs
propriétaires et les personnes contact (en date du 31/01/2025) sont listés dans le tableur "index.xlsx".

*** Structure de l'archive ***

l'architecture des dossiers est la suivante : 

/racine
//chapitre_X
///donnees
///scripts
///outputs
///chapitre_X.Rproj

Il y a un dossier par chapitre (sauf le chapitre 1 qui ne nécessite aucun jeu de données). Tous les dossiers 
contiennent les mêmes sous-dossiers + un fichier Rproj qui s'ouvre sous RStudio. Le dossier "outputs" est vide 
initialement, il se remplit lorsqu'un graphique est sauvegardé mais peut être vidé sans dommage. Un fichier 
.Rhistory (une archive temporaire) se crée lorsqu'on lance les scripts et peut être supprimé sans dommage. 

----Usage : 

Ouvrir le fichier Rproj sous RStudio. Accéder au sous-dossier "scripts" depuis l'onglet "Files" de RStudio 
(habituellement en bas à droite). Ouvrir le script depuis cet onglet. Vous pouvez alors lancer les commandes.

En tête de script une liste de paquets de fonctions s'ouvre avec la fonction 'pacman' du paquet {pacman}. Ce 
paquet doit être installé au préalable (commande d'installation dans chaque script, elle peut être neutralisée 
par un # une fois pacman installé). Les paquets de fonctions sont indispensables au bon fonctionnement du script.

Nous n'avons pas cherché à optimiser ces scripts ou les sorties graphiques, mais à utiliser le codage qui nous parait
le plus simple et le plus intuitif. De nombreuses ressources existent en ligne pour construire des scripts plus optimisés.

Nous recommandons d'utiliser la version la plus à jour possible de R et de RStudio.

----Avertissement : 

L'architecture des dossiers, leur nom, leur contenu et le nom des fichiers ne doit pas être modifié, à défaut de 
quoi les scripts ne fonctionneront plus correctement. 

Le script peut être ouvert dans n'importe quel éditeur de texte non interprété (type Bloc Notes) en dehors de RStudio, 
et les commandes envoyées manuellement sous R par des copier-coller. Il vous faudra néanmoins désigner le dossier du 
chapitre voulu comme répertoire de travail en tête de script (fonction setwd(), explications dans le chapitre 3 du livre). 









