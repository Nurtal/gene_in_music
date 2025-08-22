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
