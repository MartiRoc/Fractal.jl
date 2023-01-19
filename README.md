# Fractal.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://MartiRoc.github.io/Fractal.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://MartiRoc.github.io/Fractal.jl/dev/)
[![Build Status](https://github.com/MartiRoc/Fractal.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/MartiRoc/Fractal.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/MartiRoc/Fractal.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/MartiRoc/Fractal.jl)

## Description / Documentation:  

Ce package a pour fonction principale **JFractalMR()** qui permet à partir d'un complexe c donné en paramètre de générer une approximation de l'ensemble de julia de paramètre c. \
De nombreux paramètres graphiques sont disponibles pour rogner / déformer / changer les couleurs de l'image générée. Ce package a été créé sous la version 1.8 de Julia.

Les fonctions du package ainsi que leurs paramètres sont exhaustivement décrits dans les commentaires du script disponible à l'adresse :

https://github.com/MartiRoc/Fractal.jl/blob/master/src/Julia_fonctions.jl

**PointDeDepart** ligne 17, \
**Julia1** ligne 44, \
**Julia2** ligne 75, \
**Julia3** ligne 112, \
**Julia4R** ligne 170, \
**Julia4Rv2** ligne 243, \
**VecToMat** ligne 353, \
**MatToImage** ligne 420, \
**JFractalMR** ligne 505.

ATTENTION : Seule la fonction **JFractalMR()** peut être appelée sans utiliser le préfixe "**Fractal.**". Toutes les autres nécéssitent cet identifiant. 

## Installation:

### 1-ère Méthode par le RPEL uniquement

Taper "]" puis *Entrée* dans le RPEL de Julia (julia>) pour accéder au package manager mode du RPEL (>pkg), taper ensuite :

`add "https://github.com/MartiRoc/Fractal.jl.git"`

Cela va ajouter le package à l'environnement de travail et le précompiler (cela peut prendre quelques minutes). Attention ce n'est pas tout à fait l'adresse de ce dépôt, ne pas oublier le ".git" à la fin. Sortir ensuite du package manager mode (*Ctrl + C* sur windows, ou *Backspace* devant >pkg). Taper ensuite dans le RPEL de Julia, 

`using Fractal`

puis *Entrée*. Le package et les fonctions qu'il contient sont maintenant disponibles.

### 2-ième Méthode en utilisant Pkg.jl

Dans un script ou dans le RPEL de Julia, resp. insérer ou enchainer les instructions suivantes : 

`using Pkg` \
`Pkg.add(url = "https://github.com/MartiRoc/Fractal.jl.git")` \
`using Fractal`

Comme pour la 1-ère méthode, les deux premières instructions servent à ajouter à l'environnement de travail le package et le précompiler, cela peut prendre quelques minutes. 

### Remarque

Une fois que le package a été ajouté à l'environnement de travail, seul l'instruction `using Fractal` est nécéssaire pour commencer à l'utiliser dans cet environnement. 
