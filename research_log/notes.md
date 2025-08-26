# Research Notes

## Reprise du machin
[DATE] -> 20/06/2025
[CONTENT] -> Portage du repo sur ma nouvel infra linux, analyse des resultats de l'hyperparamétrisation
laissent penser qu'une piste audio longue donne de meilleurs résultats, je lance un test avec des samples
de 10 secondes, on verra bien (par contre la création du dataset prends des plombes, faut que je vois comment
je peux améliorer ça). Et Si je faisais de la concaténation des fichiers audio ? histoire de voir l'intégralité
du génome et de benchmarquer là dessus.
[EOF]

## Rereprise du machin
[DATE] -> 21/08/2025
[CONTENT] -> Ouai alors là je nage un peu, je me souviens juste que le projet avant un ptin de potentiel.
J'ai encore eu une nouvelle infra linux entre temps, je redéploi tout ça avec l'ambition d'avoir une petite
demo sur des faux jeu de données, histoire de faire tourner tout le pipeline et de voir comment ça fonctionne
[EOF]

## Ajout du mode demo
[DATE] -> 22/08/2025
[CONTENT] -> Ajout du mode demo dans le run.py, ça commence à être un peu le bordel dans ce fichier, il y a surement pas mal de code
qui ne sert plus à grand chose. L'idée maintenant ce serait d'intégrer une fonction dans run qui fonctionne sur de vraies donnée,
on pourrait par exemple imaginé que si on detecte le passage d'un fichier de configuration en argument on trigger la fonction avec les
paramètres contenus dans le fichier, pour ça on a besoin d'un petit mais réel jeu de données.
[EOF]

## Ajout du mode demo
[DATE] -> 22/08/2025
[CONTENT] -> Fin de chantier pour cette semaine, j'ai implémenté une fonction pour un cas d'usage réel, pour l'instant focus sur de la classification
binomiale avec des données réelles qu'on peut tirer programatiquement depuis kaggle, ça c'est pas mal, pas besoin de se les trainer sur le repo et
facile à générer / formater à la volée. Le truc c'est que faire un graphe sur tout les gènes j'ai comme dans l'idée que c'est un peu une connerie
et que ça va prendre des plombes (je vais quand même lancer un job avec 3000 genes pour voir si il tourne), pour les besoins du dev je suis partie
sur une selection aléatoire de genes, sauf qu'ils ont peu de chances de former un graphe du coup (d'avoir un truc en commun) à moins de baisser
à mort le threshold de stringDB. Faut que je trouve un système pour gérer les gènes qui se retrouvent pas dans le graphe et un système pour faire 
un découpage du jeu de données que ai du sens en  terme de pathway, je crois que j'avais déjà commencé à bosser dessus, il y a peut etre des choses à 
scvanger de ce coté là.
[EOF]

## update du run
[DATE] -> 26/08/2025
[CONTENT] -> J'ai bossé sur une nouvelle mouture de la fonction run du module run, l'idée est de faire le calcul de l'ordre des gènes
en s'appuyant sur des ressources téléchargées directement depuis stringdb plutot que de passer par l'API et de se faire bouler quand les
requêtes sont trops grosses. Je test un pipeline complet sur les données reelle de Kaggle -> besoin de rien, on peut tout télécharger.
Bcp de gene donc l'algo de calculs de position prends du temps, je sais pas si ça se parallelise, de tte façons on va faire tourner ça sur la nuit.
Je commence à me demander quelles pourraient être les applications au single cell rnaseq, j'ai l'impression que c'est ce qui en vogue en ce moment,
ça vaudrait le coup d'y jeter un coup d'oeil 
[EOF]
