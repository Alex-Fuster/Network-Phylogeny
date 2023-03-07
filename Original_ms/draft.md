---
title: "Macro-evolution doesn't really care if you're friendly or not"
geometry: margin = 2cm
latex_engine: xelatex
header-includes:
 - \usepackage{caption}
 - \usepackage{graphicx}
 - \usepackage{float}
 - \usepackage{setspace}
---
<!--\onehalfspacing-->
\doublespacing

# Introduction

Since the early stage of evolutionary ecology, species interactions, defined here as any possible interaction between to or more species such as predation, competition or mutualism, have been integrated into theory.
Erlich and Raven (1964) (*RRR*) already included interactions inside their escape and radiate model. They did highlight the importance of trophic interactions on diversification rate plant/herbivore relationships. Few years later Van Valen (1973) (*RRR*) developed his *Red Queen hypothesis* in which interactions, such as competition or predation, were integral part of the equations. More recently, Thompson (2005) (*RRR*) based his *Geographic Mosaic Theory of Coevolution* on the spatial variation of interactions (what he cold hot and cold spots) to explain differential evolution between populations of the same species.
Despite the ubiquity of interactions on evolutionary theory, their imprint on communities evolution (i.e. macroevolution) is not always detectable and straightforward to decrypt. Most of species interactions are weak and rapidly fluctuate over time. Based on these facts, two conclusions can be formulated. First, because interactions vary largely and that their strength is most of the time very weak, their imprint get lost in few generations. This is particularly true when working with fossils or phylogenies (*RRR*). The opposite thinking can be done. Even if interactions are weak, their impact at one time in the community can be so strong that it will mark it an will strongly influence the evolution of this community (*RRR*). Even if there is no one unique voice proclaiming if we can observe or not signature of interactions, their implication on macro-evolution direction can not be contested.

\bigskip

Interactions, depending on their type (i.e. positive or negative), generate different effect and out-put on communities evolution, which can be difficult to distinguish.The amount of studies on both types of interactions is unbalanced. Most of studies have investigate the impact of negative interactions on community structure and evolution,an especially competition. It's only quite recently that the interest for the influence of positive interactions, always focused on mutualism, grown up. Most of studies on negatives interactions converge towards the same results. Competition has a strong selective strength, with a  large impact on species fitness, causing traits displacement (Silvestro et al. 2015) and an increase of evolutionary rate (*RRR*). Conversely, studies about mutualistic interactions haven't reach a consensus (*RRR* <!-- review - Weber ?-->). Depending on the model used or the species you look at, it's possible to have really contrasting results.  
* Peralta 2016 : effect of phylogenies on species interaction and network structure also dependent on interaction type ; networks evolve  
* network struture (nested for mutualism ; modular for competition - Thébault & Fontaine 2010)  
* Weber et al. 2017  
Differences community structure and evolution are visible when looking at it through different interaction type.

\bigskip

To study macro-evolution, tree main approaches were developed. Traits, phylogeny and more recently network approaches are usually used to track evolution and a combination of at least two of these approaches is becoming more and more common. Traits and phylogeny are the most used approaches talking of evolution. The evolution of traits leaves a strong and easily detectable signature on community evolution, such as trait convergent, divergence or nothing. Aside of that, phylogeny have been widely used for evolutionary studies (communities related or not), with the possibility of covering long period of time. In the 2000's, based on the fact that traits and phylogenetic trees structure could be very similar, models and experiments start to couple phylogenetic tree structure and the evolution of species traits (*RRR*). One way on doing that is to look at evolution through the eyes of phylogenetic conservatism. This approach is also used as a proxy to predict (or guess) species interaction (*RRR*). Indeed, observable interaction networks at a specific time are a reflection, among other, of this conservatism. This makes phylogenetic approach a staple tool to follow community evolution, and a possible link between traits and interaction network approaches. Apart of that, networks have been used first to understand organization of communities, but it's only more recently that they start to be used in a context of evolutionary studies (*RRR*), with for instance Peralta (2016) (*RRR*) showing that ecological, instead of being statics, evolves over time.
Even if the relationship between the three approaches (traits, phylogeny and networks) is implicit for most of evolutionary studies, only few have tried to explicitly combine them and highlight the potential influence they have on on the others (*RRR*). *Ending sentence*


<!--![]
legende figure conceptuelle
In ecological communities, evolution can be followed using different approaches. The most used are phylogenetic and traits approaches and since more recently, network approaches. We know that the relationship between these approaches are still unclear, even less when we decide to explicitly include interaction in the equation.
-->

\bigskip

In this paper we want to investigate the influence of interaction type on communities evolution. We will take into account only positive and negative interactions, such as facilitation and competition respectively for instance. More specifically, we want to highlight the potential differences in the evolution of phylogenetic, traits and network structures that interaction types could induce at the community level. To do so, we will use a stochastic model based on speciation-extinction dynamic, which rebuild a community from one species and allowing the establishment of new species based on the number of there interactions. **Main result of the paper**


# Material and Methods

## Model

### General information
We build a origination-extinction model based on a niche model to track the signature of interactions type on community evolution, on unipartite communities (competition and facilitation). This model allows us to follow community evolution on a long-term perspective, where one time step represent a certain amount of time that we can consider as a generation.

In this model, we differentiate positive and negative interactions, respectively competition and facilitation, through diversification rate. This key rate is defined by two probabilities: origination and extinction probabilities. The shape of these probabilities differs depending on the interaction type, as we will explain later.
For each simulation, we have the choice to have only positive or negative interactions in the community. The community is defined as a linear continuous space limited between zero and one (i.e. niche space). During each simulation we record information about each child species entering the community: its position into the linear space and its interactions with other species, the time it appears and disappears. These information allows us to rebuild phylogenetic trees of a specific community, the possible interactions between each species and the evolution of their traits.

### Model functioning
Each species of the community is characterized by two values that we call traits: the niche optimum and the interaction range. The niche optimum indicates the position of each species in the niche space. It's defined by a numerical value between 0 and 1. The interaction range provides an area inside which each present species will interact with the focused species. Its centre is the value of the niche optimum.  
Each community start with one unique species for which the position on the niche space is picked randomly between 0 and 1. From this unique species, child species will emerge, inheriting its parent traits.

To be able to see a child species establish into the community (i.e. origination probability, *fffig1*), there is two embedded conditions that have to be fulfilled. We take separately every species present in the community, and test (i) the probability of  these species to speciate and (ii) the probability of new formed species (i.e. child species) to establish in the community.
The probability of speciation is a fixed number defined on the parameters. If speciation is possible, two potential species will be tested to enter into the community. The probability of establishment is then calculated, depending on the number of interactions the potential child species would have. This number is obtained by counting the number of species having their niche optimum trait value in the interaction range of the focal species.

The probabilities of establishment within communities with only positive or negative interactions are defined by the following equations:
$$
P(establishment_{pos}) = u_{0pos} + u_{1pos} . e^{-a_u . I_m},
$$

$$
P(establishment_{neg}) = u_{0neg} + u_{1neg} . e^{-a_{uneg} . I_{m}},
$$

where $I_m$ is the sum of interactions a species have at a specific time step, $u_{0pos}$ and $u_{1pos}$ defined how fast and high the establishment probability will increase, $u_{0neg}$ and $u_{1neg}$ defined how fast and high the establishment probability will decrease.

Every new species receive two trait values (niche optimum and interaction range), inherited from its parent species with some variation. The child niche optimum is randomly picked into a normal distribution which have the parent niche optimum value as mean and a fixed variation (see table 1, $sd$). The child interaction range value is determined as the parent interaction range value more or less a number randomly picked in the normal distribution around 0 and with a variance of 0.2.

Once all the species have been tested for origination, we test all species, even newly established species, for extinction. For each species, the extinction probability is calculated based on the number of interaction a species have (*fffig2*). The probability differs with the interaction type. If the extinction probability of a specific species is higher than a randomly picked number on a uniform distribution, than the species is removed from the community, it goes extinct, and the time step of removal is recorded.

The extinction probabilitis within communities with only positive or negative interactions are defined by the following equation:
$$
P(extinction_pos) = e_{0pos} + e_{1pos} . e^{-a_{epos} . I_{m}},
$$

$$
P(extinction_neg) = e_{0neg} . (1 - e^{-a_{eneg} . I_{m}}).
$$

**Définition des termes des équations**

\begin{figure}[H]
\centering
  \includegraphics[width = 100mm]{figures/facil_extinct_estab_prob.png}
\caption{Positive interactions communities probabilities. Black line represents origination probability, green line represents extinction probability.}
\end{figure}

\begin{figure}[H]
\centering
  \includegraphics[width = 100mm]{figures/compet_extinction-estab_prob.png}
\caption{Negative interactions communities probabilities. Black line represents origination probability, green line represents extinction probability.}
\end{figure}


For each time step, the same framework is followed: test for origination, then test for extinction. A simulation stop when it has ran for a certain amount of time steps, or when a maximal number of species in the community has been reached ($S_{max}$).

## Simulations
We ran our model and kept 100 simulations per interaction type. This amount was sufficient to catch relevant information from the output. Only simulations with communities composed of more than 20 species (when the community had reach a stable species richness) were kept for analysis. We already observe a difference between positive and negative interactions communities. Where in negative interaction communities, *%* of simulations were successful, while only *%* were accepted for positive interactions communities.
At the end of each simulation, four output were extracted: (i) species present at each time step, (ii) traits values of each species (niche optimum and interaction range), (iii) establishment time of each parent and child species, and (iv) extinction time of each species.

\bigskip

This model has 21 parameters. Some are common for the two types of interactions, some are specific to one or the other (*table 1*). Shared parameters were, for instance, the maximal amount of interactions per species (*I_max*). It allows us to avoid logistic explosion of species richness in the communities, which saturates the model and makes the simulations endless. The total number of created species during each simulation (*S_max*) is also a way of avoiding endless simulations. It has been setted at 1000 species per simulation, which is largely enough to have the emergence of stable communities in term of species richness.
The speciation probability (*u_max*) is also common to all simulations/interaction types. It's the same for every species and fixed in time at 0.15, which is around 1 speciation each 7 time step.
Concerning the extinction and establishment probabilities, parameters allow us to shape curves of these probability in a way that we can obtain stable communities. Details of these parameters are in *table 1*. All the values appeared to me in a dream, because only values matter (Perceval, 1134). Parameters related to these probabilities were determinant for the success of simulations.

\bigskip

Table 1 : Parameters, signification and values.  
:warning: *$u_1-u_0$ and $e_1-e_0$ are the establishment and extinction rates in absence of interactions.*

| Parameter | Signification                                                                                                              |Value|
|-----------|----------------------------------------------------------------------------------------------------------------------------|-----|
|$int$      |Interaction type, 0 for competition, 1 for facilitation                                                                     |     |
|$S_{max}$  |Maximal number of species in the system                                                                                     |1000 |
|$av_r$     |Half range of the niche of the first species                                                                                |0.1  |
|$sd$       |Standard deviation of the normal distribution used to calculate the niche optimum trait of a child species      |2*$av_r$ + 0.0001|
|$u_{max}$  |Speciation probability                                                                                                      |0.1  |
|$a_{eneg}$ |Shape of the exponential decay of the negative extinction - interaction relationship                                        |0.025|
|$e_{0neg}$ |Asymptotic extinction probability with infinite negative interactions                                                       |0.5  |
|$e_{1neg}$ |Extinction probability with absence of competition interaction                                                              |1    |
|$a_{epos}$ |Shape of the exponential decay of the positive extinction - interaction relationship                                        |1.2  |
|$e_{0pos}$ |Asymptotic extinction probability with infinite positive interactions                                                       |0.075|
|$e_{1pos}$ |Extinction probability with absence of facilitation interaction                                                             |5.19 |
|$a_{uneg}$ |Shape of the exponential decay of the colonization - interaction relationship                                               |0.075|
|$u_{0neg}$ |Asymptotic establishment probability with infinite competition interactions                                                 |0.075|
|$u_{1neg}$ |Establishment probability with absence of competition interaction                                                           |2    |
|$a_u$      |Shape of the exponential decay of the colonization - interaction relationship                                               |0.45 |
|$u_{0pos}$ |Asymptotic establishment probability with infinite facilitation interactions                                                |1    |
|$u_{1pos}$ |Establishment probability with absence of facilitation interactions                                                         |-1   |
|$d$        |Decrease speed of the establishment probability                                                                             |0.5  |
|$B_{spe}$  |Constant minimal number of interaction per species (establishment prob)                                                     |4    |
|$B_{ext}$  |Constant minimal number of interaction per species (extinction prob)                                                        |4    |
|$I_{max}$  |Maximal number of interacting species                                                                                       |40   |


## Analysis

With this model, we wanted to see if we were able to identify a signal of interaction type on the evolution of communities. First, we need to be sure our model was working, meaning we were able to have stable communities for both, positive and negative, interactions. To do so, we followed the species richness dynamic of every simulation and calculated the mean species richness dynamic and its 95% confidence interval. This calculation has been done for almost all indices showed in this paper (*fig. 3*).

### Phylogenetic analysis
Because we tracked which species is the parent of any specific species of a community, we can build the real phylogenetic tree of each community. Thanks to this tree, we calculated the $\alpha$-value at each time step, for every community. The $\alpha$-value is a measure of diversification rate, based on the $\gamma$ statistic (*RRR*), which allows us to compare different phylogenetic trees. It compares the branches length at the top and the bottom of a tree. If top branches are relatively longer than bottom ones, $\alpha$ is higher than 0. This means that the diversification process increased/occurred recently in the community history. The opposite, $\alpha$ lower than 0, means that diversification was faster at the beginning of the community history.
For the analysis, we built bifurcation trees with the information on parentality and extinction time. It also allows us to extract branch length of every species in communities, which indicates species *age*.

### Network analysis
We calculated the connectance of species interaction networks at each time step, for each community. The connectance gives us the relative number of realized interaction (i.e. link) in the community in regards to all the links we could have between the species of a community. As our communities form unipartite and directed networks (i.e. an interaction goes from a specific species to another, the reciprocal is not always true, giving non-symmetrical adjacency matrices), we used the following equation to calculate connectance:

$Co = L/S^{2}$,

where $L$ is the number of links and $S$ the number of species of the community.


We also calculated the degree, $k$, of each species at each time step (i.e. the number of interaction a species have). For example, to find the  degree of a species $i$, we take the amount of species inside the interaction range of this focal species at a specific time step. This is possible, by gathering together information of species presence and trait values of species.  
For both, connectance and degree, we calculated the mean and the 95% confidence interval, for each simulation.

### Traits analysis
We looked at the evolution of the niche occupancy and the distance between species inside the niche space.
*As to be done*.

# Results

## Model evaluation - Vital functions check up

### Species richness dynamics

Le type de courbe obtenue pour la dynamique de la richesse spécifique, que ce soit pour les communautés d'interactions positives ou négatives, sont normal considérant le modèle utilisé. Le fait que les phases de croissances pour atteindre le plateau au niveau du nombre d'espèces est étroitement lié à la probabilité d'établissement des nouvelles espèces. Au début des simulations, cela met bien en valeur 3 processus ayant des directions opposées.
- Chez les com. pos. l'établissement de nouvelles espèces est plus probable quand le nombre d'interaction de l'espèce arrivant est plus élevé. Il sera donc plus difficile pour les espèces arrivant au début de l'histoire de la communauté de l'établir, étant donné qu'il y a moins d'espèces disponibles avec lesquelles interagir. (P(est.))
- Plus le nombre d'espèces augement dans la comm. plus la possibilité de spéciation augmente également, de façon exponentielle
- ?

C'est à mettre en relation avec la différence de succès lors des simulations (les débuts difficiles des com. pos.).

com. pos. => au début, on va avoir un mécanisme poussant à garder les espèces avec des traits similaires
com. neg. => au début, garde les espèces avec des traits différents

On peut aussi constater que la variation dans les simulations des com. pos. est beaucoup plus forte que dans celles des simulations des com. neg.. Quand les comm. ont atteint une stabilité au niveau du nombre d'espèces, cette variation devient similaire pour les deux types de communautés.

On peut donc voir grâce à cela qu'il y a un déterminisme plus important pur les com. neg. que pour les positives.

### Alpha-statisitque

Au début de l'histoire de la communauté, les branche sont assez longue car il y a peu d'espace de niche disponible (relativement parlant, la probabilité qu'une espèce fille apparaisse proche de son ancêtre est plus forte que l'apparition de traits complètement différents, donc la niche est déjà occupée localement). Par la suite, on observe une *perte d'héritage* des traits dans la communauté, une diversitfication des traits et par conséquent, de la communauté.

**À voir** stat. K ou lambda de *Blumbers*  
Description des traits dans la phylogenie est-elle conservée ou aléatoire ?

**Lien phylogénie - trait**  
Au début, observe t-on un héritage évolutif ? et avec le temps, une perte de cet héritage, caractérisé par une structure brownienne de la distribution des traits dans les clades ?

**Lien phylogénie - réseaux**  
Observe t-on moins de hasard dans les réseaux ? les espèces proches physogénétiquement ont une niche similaire et interagissent entre elles. Est-ce qu'on a cette perte de relation là dans le temps ? (conservatisme phylogénétique)

idée basée sur les analyse du projet de Madelaine:

- PCoA sur les distance phylo
- PCoA sur les distance dans les réseaux
- Procrust pour relier les deux

voir si c'est faisable et que ça donne qqch en regardant (i) à différent moments dans le temps ; (ii) en prenant toutes les espèces, même celles éteintes ; (iii) en ne prenant que les espèce qui ne sont pas éteintes.
(ii) et (iii) permettent de voir l'état de l'occupation de l'espace dans les réseaux.
=> Au début, les espèces sont proches les unes des autres (réseau relatif), avec le temps il y a plus d'occupation de l'espace réseau.

**voir le procrust à différents pas de temps + cumulatif**

### Connectance et degré

**Vérifié la calcul de la connectance (pas supposé être plus grand que 1) voir ce qui ne fonctionne pas**

**Redimensionner l'axe y (c(o:0.5))**

*À placer après le pannel de traits et de degré*
Pour les com. neg., la largeur moyenne du range d'interaction diminue avec le temps, cela impacte donc directement le nombre d'interactions que les espèces pourront avoir. On observe donc une diminution du degré moyen pour ces communauté. De ce fait, en plus d'augmenter le nommbre d'espèces dans la communauté et donc la probabilité de réaliser des liens, cela restreint la connectance, qui diminue rapidement et de façon très claire.

Pour les com. pos., la diminutions n'est pas aussi nette, et n'est pas la conséquence des même mécanismes. Ici, la diminution de la connectance est dûe à l'augmentation du nombre d'espèces dans la communauté, et donc l'impossibilité croissante d'interagir avec toutes les espèces présentes. La diminution moins lisse de la connectance montre une diversité de configuration de réseaux plus grande que pour les com. neg. .

**Ajouter les écarts-type pour les degré moyens**
- diversification ou conservatisme des structures au cours du temps ?
- on pourrait regarder le diamètre ou le *rayon spectral*

**À voir et lire** *Stanichenko** => l'effet de la compétition sur le rayon spectral, étude qui n'a pas été fait sur des simulations, donc ça serait intéressant de voir si on obtient le même genre de résultats ou non.
En gros, ils font une mesure de la nestedness sur des réseaux unipartites.

\begin{figure}[H]
\centering
  \includegraphics[width = 150mm]{figures/basic_figures_100.png}
\caption{Communities evolution considering species richness and the three approaches. Basic indices,such as the species richness dynamic, the $\alpha$-value and the connectance, don't show any real difference in the shape of their evolution. Positive interactions communities seem to reach the carrying capacity slightly slower than the negative ones. A little difference is also visible with the $\alpha$-value, where the minimum value is in average higher for positive interaction communities than for negative ones. However, there is a selective pressure on the number of interactions and the interaction range related to the interaction type of the community. The mean interaction range trait value increase over time in positive interaction communities, while the opposite is observed for negative interaction communities. The mean number of interactions per species also differs from positive to negative interactions communities, it increase over time for communities with positive interactions and decrease for negative ones.}
\end{figure}

## Age-Degree

Normalement, on s'attend à observer une augmentation du degré avec l'age dans les com. pos.. (plus on a de degré, plus on a de chance de persister dans la communauté). Mais ici, ça n'est pas ce que l'on observe.

Avec notre modèle, on a une espèce d'héritabilité dans le nombre de degré, cela nous donne des espèces jeunes avec un degré élevé.

Ce qui est important ici, c'est le haut de l'enveloppe, la forme, pas ce qu'il y a dedans.

Pour le vioplot, ça explique le alpha qui tend vers 0.

**À voir**

- Distribution des branches qui ses sont éteintes sans avoir fait de spéciation. Quelle est la différence avec que l'on a là ?

**À faire**

Vioplot tous les 10 pas de temps, pour voir l'âge moyen (min et max). Ça nous donne une information sur le turn-over des espèces dans les communautés.
Normalement, on devrait avoir plus de turn-over dans les communautés négatives que dans les positives. Ça peut aussi être un proxy ou une équivalence du taux de diversification (ou au mins, ça nous permet de faire un lien)

Ce que ça pourrait nous permettre de faire c'est de montre que la structure phylogénétique n'est pas un processus évolutif stable dans le temps, ce qui est souvent considérer dans les études phylogénétiques (**à vérifier**).  
Avec ça, on pourrait contrait un agrument souvent utilisé : Si on a un changement dans la structure de l'arbre, on a un ou des changements dans les processus de sélection. Car nous observerions un changement dans la structure phylo, mais pas dans les processus de sélection.

\begin{figure}[H]
\centering
  \includegraphics[width = 100mm]{figures/age-degree_figures_100.png}
\caption{Relationship between species degree and branch length/species age. The distribution of species degree is related to the interaction type of the community. Species from negative interaction communities tend to have in average less interactions than species from positive interactions communities. This is can be a consequence of the interaction range trait value evolution difference between positive and negative interactions (\textit{fig.3}).}
\end{figure}

## Niche occupancy

**Voir les statistiques spatiales sur les patrons de point** (à quelle point est ce que la répartition des traits est aléatoire ?)

Avec notre modèle, on a des mécanismes qui favorisent l'uniformisaté.

Donc avec ces statitiques, on pourrait avoir une description de la niche sans prendre en compte l'effet combiné de la position et du range d'interaction (les deux traits attitré à chaque espèce). DOnc indépendant dans la largeur du range d'interaction, est-ce qu'on a bien une distribution uniforme ?

\begin{figure}[H]
\centering
  \includegraphics[width = 150mm]{figures/figure_niche_occ.png}
\caption{Evolution of niche occupancy and distance between species over time. The evolution of niche occupancy over time indicates that the niche is homogeneously filled from the early time of community history. Once again, there is no major difference in the shape of this evolution depending on the interaction type. The same lack of shape difference is observed for the distance between species in the niche space. This means that, even if it creates selective pressure on interaction range trait value, our model can not differentiate trait structure depending on interaction type.}
\end{figure}

<!--# Discussion

# Conclusion
--->
