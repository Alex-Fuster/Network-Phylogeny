---
title: "Projet Macro-évolution"
geometry: margin = 2cm
latex_engine: xelatex
header-includes:
 - \usepackage{caption}
 - \usepackage{graphicx}
 - \usepackage{float}
 - \usepackage{setspace}
---

# Introduction

1. Communautés et théories de l'évolution

2. Richesse spécifique et phylogénie des communautés comme indicateurs de l'évolution des communautés

3. Structures de traits et des réseaux d'intéractions, lien avec la phylogénie

4. Objectif : Détection de la signature du type d'interactions dans l'évolution des communautés

\bigskip

# Methods

\begin{figure}[H]
\centering
  \includegraphics[width = 150mm]{Figures/conceptual_figure.png}
\caption{Conceptual illustration of the model used for the simulations. (a) Important information about the model. The niche space is defined as a linear and continuous space between 0 and 1 (a1). Each species is defined by two numerical values : the niche optimum which is a number between 0 and 1 and the interaction range, centred on the niche value (a2). There is an interactions if the niche value of a species is located into the interaction range of another one (a3). Once that is set, the model can run following always the same frame (b). First, each species is tested for origination (speciation and establishment of new species)(steps 1 to 3). Second they are tested for extinction (step 4). From the time n to the time n+1, there are 5 different possible situations : 1. the species stay, 2. it goes extinct, 3. it speciates and the two news species establish, 4. and 5. it speciates but only one new species establishes.}
\end{figure}

\bigskip

# Results

\begin{figure}[H]
\centering
  \includegraphics[width = 100mm]{Figures/community_dynamic100.png}
\caption{Community dynamic. Green represents indices calculated for positive interaction communities, Orange is used to illustrate indices calculated for negative interaction communities. A. illustrates species richness over time. Both interaction type community show the same type of structure (logistic shape) in the evolution of their species richness. B. illustrates the evolution of species turnover depending on time. Here again, the two interaction types show the same structure.}
\end{figure}

Figure 3 : Phylogenie et traits

Figure 4 : Phylogenie et réseaux

\bigskip

# Discussion

1. Pour obtenir des communautés qui étaient capables de se maintenir dans le temps, nous avons été obligé de définir des paramètres très différents pour la probabilité d'établissement des nouvelles espèces et celle d'extinction. Turnover similaire alors qu'on pourrait justement s'attendre à des différences, dû aux différences de paramètres.

2. Malgré le fait que l'on observe de la pression de sélection sur les traits et les interactions entre espèces, pas de lien visible entre ces indicateurs pourtant lien mis en évidence dans la littérature (conservatisme phylogénétique, etc.)

3. Modèle qui sépare distinctement les types d'interaction que l'on trouve dans les communautés, mais ce n'est pas qqch que l'on retrouve dans les environnements naturels ni en laboratoire.

4. En conclusion, avec ce modèle, pas possible de détecter une signature du type d'interaction dans l'évolution de communautés.
