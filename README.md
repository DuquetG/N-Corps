# Le problème à _n_-corps: simulations de systèmes stellaires

## Aperçu

Il est bien connu que l'univers est en constante évolution au travers du temps. Lorsqu'on étudie un système stellaire pour de grandes échelles de temps, la moindre perturbation peut gravement modifier son allure de manière chaotique. Ce programme a comme but d'étudier comment sont modifiés les trajectoires des objets célestes d'un système à l'équilibre ou non de _n_-corps en raison de différents jeux de paramètres, configurations et pertubations

Ce projet a été réalisé dans le cadre du cours de __Physique Numérique - Projet (MU4PY108 - S1-23)__ du programme M1 de physique de la Sorbonne Université.

Collaborateurs:

- __ARSENAULT__, Marek
- __BRY__, Charles-Éric
- __DUQUET__, Gabriel

## Initialisation du code

### Créer un environnement virtuel
```
python -m venv venv
```
** À noter ici que si la commande `python` ne fonctionne pas, il se peut que l'utilisation de la commande `python3` corrige le problème. Cela est souvent du au fait que deux versions de python sont installées simultanément sur la console.

### Activer l'environnement virtuel
__Pour Linux/macOS :__
```
source venv/bin/activate
``` 
__Pour Windows :__
```
venv\Scripts\activate
```

### Installer les exigences dans l'environnement virtuel
```
pip install -r requirements.txt
```

## Fonctionnement général du code à deux dimensions

Le code commence par lancer une boucle qui itère autant de fois que le nombre de simulations. Dans la boucle, le code apelle d'abord le sous-programme `Simulation2D` pour créer, dans des fichiers `.csv` ou `.dat` les listes des valeurs selectionnées pour chaque pas de temps. Pour sélectionner les listes valeurs à produire, voir la section __Paramètres - Booléens__.  
  
Si vous choississez de produire les trajectoires et/ou les vitesses à tout moment, les valeurs de chaque simulation seront placées bout à bout dans un fichier `.csv` rassemblant toutes les trajectoires ou toutes les vitesses.

Le code traitera ensuite ces fichiers avec la fonction `heal` de `csv_healer.py` pour concatener les _n_ lignes de chaque simulation sur la même ligne et ainsi de suite. Cela permettra de superposer les simulations pour faciliter la visualisation des différences de trajectoire grâce à `Animation2D.py`. 

Avant tout, cela facilite le traîtement des données pour l'utilisation du sous-programme `Lyapunov`. Ce dernier effectue le calcul de l'exposant de Lyapunov pour chaque simulation.

Finalement, le code appelle `graphiques.py` pour produire les différents graphiques utiles à l'étude du système.

Nous recommandons d'indenter les requêtes d'execution de commandes inutilisées dans le `main.f90` pour éviter d'alourdir le code.

Pour comprendre le fonctionnement des programmes en particulier, nous recommandons d'aller lire la documentation propre à chaque programme.

__Pour lancer le code, insérez le code suivant dans votre terminal sans le '+':__

```diff
+ gfortran main2D.f90 -o main2D
```

## Paramètres

Les paramètres sont libres d'être modifiés pour étudier les intéractions du système pour différentes combinaisons.

### Entiers (integer):
- Nombre de corps : _nbCorps_
- Nombre d'itérations : _Nstep_
- Masse dont l'exposant de Lyapunov est étudié : _chosen_mass_
- Nombre de simulations : _nbSimulations_
- Nombre de pas sautés à chaque calcul de l'exposant de Lyapunov : _pasCalcul_

### Nombres réels (real(8)):
- Temps de l'indentation en secondes : _dt_

### Listes (real(8)):
- Liste des masses : _M_ (dimension(_nbCorps_))
- Liste des conditions initiales (positions et vitesses): _X_ (dimension(nbCorps, 2*__le nombre de dimensions__))

### Booléens (logical):
- Calcul des trajectoires des corps : _wtraj_
- Calcul des fluctuations d'énergie du système : _wenergy_
- Calcul de l'énergie moyenne à chaque pas de temps : _wviriel_
- Calcul des vitesses des corps : _wvelocity_

## À savoir

- Les fichiers `.csv`, `.dat`, `.exe` et `.mod` sont ignorés.
