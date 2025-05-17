************************ Archive en ligne de l'ouvrage ****************************

Jean-Yves Barnagaud et Olivier Gimenez, 2025. XXX. Biotope Editions.

jean-yves.barnagaud@ephe.psl.eu

Cette archive contient les données et scripts R permettant la réplication des exemples et figures du livre.

Dernière mise à jour : 17/05/2025

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

Pour plus de détails, voir plus bas "Téléchargement de l'archive"

----Avertissement : 

L'architecture des dossiers, leur nom, leur contenu et le nom des fichiers ne doit pas être modifié, à défaut de 
quoi les scripts ne fonctionneront plus correctement. 

Merci de ne pas modifier vous-mêmes l'archive sur GitHub. S'il faut effectuer un changement, prévenez les auteurs par une "issue". 

Le script peut être ouvert dans n'importe quel éditeur de texte non interprété (type Bloc Notes) en dehors de RStudio, 
et les commandes envoyées manuellement sous R par des copier-coller. Il vous faudra néanmoins désigner le dossier du 
chapitre voulu comme répertoire de travail en tête de script (fonction setwd(), explications dans le chapitre 3 du livre). 

** Téléchargement de l'archive **

1 - Installer R -------------------------------------------------------------------------

!! R est un logiciel gratuit. Vous ne devez jamais payer pour le télécharger. 

Télécharger R en suivant ce lien : https://cran.r-project.org/
Dans la section "Download and Install R", choisir votre système d'exploitation (Windows, Mac ou Linux)
Cliquer sur "base" puis sur "Download R..."

Les versions de R changent fréquemment, cela n'est normalement pas très impactant mais il est conseillé de 
faire au moins une mise à jour annuelle. Pour cela, le plus simple est de désinstaller / réinstaller R. 

Sous Windows, vous pouvez faire les mises à jour de R avec la fonction updater() du package installr
Sous Mac, vous pouvez utiliser le package updateR

2 - Installer R Studio ------------------------------------------------------------------

!! R Studio est un logiciel gratuit. Vous ne devez jamais payer pour le télécharger. 

!! Vous devez avoir installé R avant. R Studio n'est qu'un éditeur de script qui exploite le logiciel R.
 
Télécharger RStudio en suivant ce lien : https://posit.co/download/rstudio-desktop/
Le lien "Download Rstudio desktop" démarrera le téléchargement de l'exécutable d'installation.

Des mises à jour R Studio sont fréquemment proposées. Elles sont peu impactantes, mais faites-les afin de rester
sur la version la plus récente. 

3 - Récupérer l'archive -----------------------------------------------------------------

	3a - manuellement (à privilégier pour les utilisateurs débutants ou occasionnels)

		(i) : récupérer l'archive

Allez sur https://github.com/jybarnagaud/analyse-donnees-livre-2025

Cette archive est mise à jour régulièrement par les auteurs à mesure des remarques et corrections qui nous parviennent.
Vous pouvez signaler un problème de fonctionnement des scripts en créant une "issue" : https://github.com/jybarnagaud/analyse-donnees-livre-2025/issues

Pour télécharger l'archive complète : cliquer sur le bouton vert "<> Code" et "Download ZIP"
Décompressez l'archive ZIP.

		(ii) : utiliser l'archive

!! Ne modifiez pas la structure de l'archive. 

Il y a deux manières d'utiliser l'archive : 

- globale : ouvrez le projet "analyse-donnees-livre-2025.Rproj" avec R-Studio. Dans la fenêtre "Files", vous pourrez accéder
à tous les chapitres. Cliquez sur le chapitre voulu, puis sur le dossier "scripts", et ouvrez le script. Vous pouvez alors 
commencer à travailler.

- chapitre par chapitre : chaque chapitre dispose de son propre projet R-Studio nommé "chapitre_XX.Rproj". Ouvrez le projet
du chapitre voulu sous R Studio. Dans l'onglet "Files", choisissez "scripts", et ouvrez le script. Vous pouvez commencer à 
travailler.

	3b - avec un compte GitHub (pour utilisateurs expérimentés)

Cette méthode permet d'importer l'archive du livre depuis R Studio et d'en suivre les mises à jour sans avoir à la 
re-télécharger entièrement. Elle n'est recommandée que pour les utilisateurs habitués de R Studio ayant déjà une 
connaissance minimale de GitHub. Pour cette raison, les deux premières étapes ne sont pas détaillées.

!! Nous utilisons GitHub et non GitLab. 

!! Utilisateurs habitués à Git: merci de ne faire aucun Commit vous-même si vous modifiez les scripts. Créez une "issue" ou 
contactez-nous par courriel s'il faut répercuter une modification. 

		(i) : créer un compte GitHub
		(ii) : télécharger Git
		(iii) configurer R Studio

La configuration de R Studio pour travailler avec votre compte GitHub est détaillée ici : 
https://gist.github.com/Z3tt/3dab3535007acf108391649766409421

		(iv) : récupérer l'archive

Depuis R Studio : File / New Project / Version Control / Git
L'URL du projet est : https://github.com/jybarnagaud/analyse-donnees-livre-2025.git

		(v) : utiliser l'archive

Il y a un projet R pour tout l'ouvrage à la racine (analyse-donnees-livre-2025.Rproj). Les données et scripts sont
répartis dans des dossiers par chapitre. Depuis l'onglet "Files" de RStudio, ouvrez le chapitre de votre choix. 

Pensez à faire tourner la ligne setwd("chapitre_XX") de la section "##options" du script, sinon le chemin d'accès sera incorrect.







