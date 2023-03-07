---
title: "Macro-evolution model"
geometry: margin = 2cm
latex_engine: xelatex
header-includes:
 - \usepackage{caption}
 - \usepackage{graphicx}
 - \usepackage{float}
---

<!---# Diversification dynamics model--->

*The script is available at : https://github.com/DominiqueGravel/ms_diversification (code folder) and in the dopbox Macroevo_Interaction (same folder)*

The aim of the model is to recreate establishment/extinction events in a community of a particular type of interactions (*i.e.* facilitation, competition or predator-prey). It allows us to access the structure of the species interactions network over time, and to build the phylogenetic tree of the community.

:warning: Here we are talking about **unipartite** networks.

<!---## Species and interactions

In this model, each species is defined by a set of three traits : $n$, $r$ and $o$.
$n_i$ represents the niche position of the species $i$, $r_i$ the range of the niche (which is fixed for all species, given by the parameter $av_r$), and $o_i$ the optimum of the niche of the species $i$.
In case of competition and facilitation interactions, a species $i$ interacts with other species $j$ whose niche position falls within the range $[n_i - r_i, n_i + r_i]$.
In case of predator-prey interactions, a predator with niche position $n_i$ feeds on preys whose niche position $n_j$ falls within the range $[o_i - r_i, o_i + r_i]$.--->

# Model functioning

Before running the model, we have to indicate :\newline
- the number of time steps ($nsteps$),\newline
- the seed, \newline
- a set of parameters (defined at the end of the document)

## Before all : Prepare tables and the <!--- basal and the ---> first species

<!---The basal species can be considerate as the network resource. To sustain, the first species has to interact with at least one of the basal species.
The number of basal species is given in the parameters ($pars\$Sbasal$). Each of those species is characterized by a random number between 0 and 0,2 (*it can be changed*). These basal species are not allowed to mutate, they cannot lead to speciation or extinction.
To start the simulations, a first species is generated. Its traits ($n$ and $o$) are randomly selected into a uniform distribution between 0 and 1, and $r$ is fixed at 0.2.

Once the basal and the first species have been setted, the 2 next steps are repeated until the maximum amount of time ($nstep$) has been reached.--->

### Setting tables $traits_{mat}$, $anc$, $extinction$ and $pres$.

The aim here is to build empty tables and matrices, which will be filled all along the simulation.

$\bullet$ $traits_{mat}$ table compiles all the traits combinations of established species. The traits of interest here (for competition and facilitation) are the traits $n$ (niche optimum) and $r$ (half of niche range).

| Trait n | Trait r | Trait o |
|---------|---------|---------|
|         |         |         |
|         |         |         |

$\bullet$ $anc$ compiles the time step at which a species establish into the community and its ancestor.

| time step | ancestor | new species |
|-----------|----------|-------------|
|           |          |             |
|           |          |             |

$\bullet$ $extinction$ compiles the species extinctions and the time step at with the extinction occurs.

| time step | extinct species |
|-----------|-----------------|
|           |                 |
|           |                 |

$\bullet$ $pres$ compiles all the present species in the community at a time $t$. Each column is a species and each row a time step. $pres$ is pre-build with $S_{max}$ columns and $nsteps$ row. It is a binary matrix where 1 means the species is present and 0 absent.

|        | Species1 | Species2 | Species3 | ... |
|--------|----------|----------|----------|-----|
| time 1 |          |          |          |     |
| time 2 |          |          |          |     |
| ...    |          |          |          |     |


Each species is represented by a number going from 1 to $S_{max}$, attributed by arrival order.

### Setting the first species

To start the simulation, a first species is generated using the function `rand_traits_anc`. Its traits ($n$ and $o$) are randomly selected into a uniform distribution between 0 and 1, and $r$ is fixed by the parameter $av_{r}$.

Then the number of species, $S$, in the community is setted to 1 (as we have only one species for now).

**Now that the model is ready to really start, the two later parts of the model (establishment and extinction) will be executed at each time step, as part as a loop, as long as $S_{max}$ or $nsteps$ has been reached.**

## Inside the loop, First step : See if new species can estabish

The first thing is to check the total number of species in the community (not only the actual number of species at the time $t$). If the limit $S_{max}$ has been reached, the model stops.

For each present species, a random number is picked on a uniform distribution between 0 and 1. Mutation is allowed if this number is below the threshold of the mutation probability ($u_{max}$), defined by

$$
P(mutation) = \frac{u_{max}}{(1 + e^(d * (S_{t} - I_{max})))},
$$

where $u_{max}$ is , $d$ a constant, $S$ species number at $t$ time/step and $I_{max}$ the maximal number of species in the community.

### Mutation allowed? Yes!

<!-- When mutation in allowed on a species $i$, new traits will be calculated from the traits of this species $i$, using the function `rand_traits_mut`. <!---To avoid the possibility of super-generalist, because there is no explicit cost for generality in this model,---> A *mutant* trait $n_{m}$ is picked on a normal distribution (Fig. 1) where the ancestor trait $n$ is the mean, and the standard deviation is defined by a number higher than the total interaction range.

$$
P(n_{m}) = \frac{1}{\sqrt{2\pi\sigma^2} } e^{ -\frac{(n_{m}-n)^2}{2\sigma^2} }
$$,

$$ \sigma = 2 * av_{r} + \epsilon $$.

\begin{figure}[]
\centering
  \includegraphics[width = 100mm]{figures/normal_distrib.png}
\caption{Examples of normal distribution used to generate the trait n of a new species.}
\end{figure}

A new value for the trait $o$ ($o_m$) is recorded as the value of $n$.
The trait $r$ is also subject to mutation. The niche range of the new species , $r_m$ is defined based on $r$ plus a number picked on a normal distribution of mean 0 and standard deviation of 0,2.

Once the *mutant* traits have been setted, we calculate the number of possible interaction this *mutant* species can have. First, the function `get_L_vec`, calculate the number of species with a trait $n$ within the range of $[n_m - r_m, n_m + r_m]$ of the *mutant* species.

Then, we test for establishment the *mutant* species. The probability of establishment is different depending on the interaction type. For facilitation interactions (Fig. 2), the probability follows:

$$
P(establishment_f) = u_0 + u_1 . e^{-a_u . I_m},
$$

where $u_0$ represents the asymptotic establishment probability with infinite interactions, $u_1$ ($u_1 = -u_0$ for facilitation interactions) is the establishment probability with absence of interactions, $a_u$ is the shape of the exponential decay of the colonization-interaction relationship and $I_m$ is the sum of the possible interactions of mutant species with the rest of the community.

$I_{m}$ is obtained by:

$$
I_{m} = \sum (I . P) + basal_{estab},
$$

where $I$ is a vector obtained with the function `get_L_vec`, $P$ is the species present at the time $t$ (from $pres$) and $basal_{estab}$ is the minimal number of basal species needed for the community not to collapse.

When simulating competition interactions (Fig. 3), the probability follows:

$$
P(establishment_c) = u_{0neg} + u_{1neg} . e^{-a_{uneg} . I_{m}},
$$

Once the probability of establishment of the *mutant* species has been calculated, a random number is picked on a uniform distribution between 0 and 1. If this random number is lower than the probability of establishment, the *mutant* species can establish into the community.

### Establishment allowed? Yes!

If establishment is allowed, the traits of new species are recorded into the $traits_{mat}$, the time of arrival and related ancestor are recorded into $anc$, the presence of the species is recorded into $pres$ at the row $t$, and the number of species $S$ is increased by 1.


## Inside the loop, Second step : See if the existing species goes extinct

Once every species of the community have been checked for speciation/establishment, we test yet for extinctions.

First, we want to know which species co-occurs and interacts with which. To do so, we start by creating a co-occurrences matrix ($cooc$) between present species at the time $t$. Then, thanks to the function `get_L_mat`, the potential interactions matrix ($L_p$) between every species in the community is calculated. Species $i$ interacts with species $j$ if $n-j$ is within $[n_i - r_i, n_i + r_i]$. And because no self-interaction is allowed, the diagonal will always be 0.
The realized interaction matrix is obtained by multiplying $L_p$ and $cooc$.

### Extinction probability

$\bullet$ In case of facilitation interaction

First, the sum of the realized interaction per species is calculated by

$$
I_{m} = \sum I + basal_{ext},
$$

where $I$ is the realized interaction per species, and $basal_{ext}$ is the minimal number of basal species needed for the community not to collapse to fast.
$I_{m}$ is the vector containing the sum of species interactions per species.

Then the extinction probability (Fig. 2) is defined as :

$$
P(extinction) = e_{0pos} + e_{1pos} . e^{-a_{epos} . I_{m}},
$$

where $e_{0pos}$ is the asymptotic extinction probability with infinite positive interactions, $e_{1pos}$ is ($1 - e_{0pos}$), $a_{epos}$ is the shape of the exponential decay of the positive extinction-interaction relationship and $I$ is the total number of interaction per species.

$P(extinction)$ is a vector containing all the probability of extinction per species.

$\bullet$ In case of competition interaction

Same thing than for facilitation interactions, the sum of the realized interaction per species $I_{m}$ is calculated. Then the extinction probability (Fig. 3) is defined as :

$$
P(extinction) = e_{0neg} . (1 - e^{-a_{eneg} . I_{m}}),
$$

where $e_{0neg}$ is the asymptotic extinction probability with infinite negative interactions, $a_{epos}$ is the shape of the exponential decay of the positive extinction-interaction relationship and $I_{m}$ is a vector with the total number of interaction per species.

### Extinction allowed?

For an extinction to occur, a random number is picked, into a uniform distribution between 0 and 1, for every species. The extinction probability is compared to each random number. If a number is smaller than the extinction probability, the related species goes extinct.

For every species $i$ that goes extinct, the extinction is recorded into $ext$, with the species number $i$ and the time step $t$ ; and $pres$ is updated, replacing $pres[t, i]$ by 0.

\begin{figure}[H]
\centering
  \includegraphics[width = 100mm]{figures/facil_extinct_estab_prob.png}
\caption{Example of establishment and extinction probabilities for facilitation interactions. establishement parameters : u0 = 0.9, u1 = 0.9, a = 0.45, d = 0.5, Imax = 12 ; extinction parameters : e0pos = 0.05, e1pos = 5.19, aepos = 0.9, Imax = 12.}
\end{figure}

\begin{figure}[H]
\centering
  \includegraphics[width = 100mm]{figures/ext-est_compet.png}
\caption{Example of establishment and extinction probabilities for competition interactions. establishement parameters : u0 = 0.075, u1 = 2, a = 0.5 ; extinction parameters : e0neg = 0.15, e1neg = -0.15, aeneg = 0.025.}
\end{figure}

### Network properties

To avoid recording every interaction matrix, we calculate the network connectance and each species degree (number of interaction per species).

Because we have a undirected unipartite network where self-interactions are not allowed, the connectance is given by

$$
Co = L / S^{(S-1)/2},
$$

where $L$ is the number of links and $S$ is the number of species in the network.

# Model output

The model output gives you access to : \newline
- $traits_{mat}$ table,  \newline
- $anc$ table,  \newline
- $extinction$ table,  \newline
- $pres$ matrix,  \newline
- $L$, the realized interactions matrix from the last time step.


# Parameters

$u_{max}$ : Mutation probability (asymptotic establishment rate)\newline
$u_0$ : Asymptotic establishment probability with infinite facilitation interactions\newline
$u_1$ : Establishment probability with absence of facilitation interactions\newline
$a_u$ : Shape of the exponential increase of the colonization - interaction relationship (facilitation)

$u_0neg$ : Asymptotic establishment probability with infinite competition interactions\newline
$u_1neg$ : Establishment probability with absence of competition interactions\newline
$a_uneg$ : Shape of the exponential decay of the colonization - interaction relationship (competition)


$e_{0neg}$ : Asymptotic extinction probability with infinite negative interactions\newline
$e_{1neg}$ : Extinction probability with absence of interactions ($-pars\$e_0neg$)\newline
$e_{0pos}$ : Asymptotic extinction probability with infinite positive interactions\newline
$e_{1pos}$ : ($1 - pars\$e_0pos$)

:warning: *$u_1-u_0$ and $e_1-e_0$ are the establishment and extinction rates in absence of interactions.*

$a_{eneg}$ : Shape of the exponential decay of the negative extinction - interaction relationship\newline
$a_{epos}$ : Shape of the exponential decay of the positive extinction - interaction relationship

$av_r$ : Range of the niche\newline
$sd$ : standard deviation of the normal distribution used to calculate the trait n of a new species

$int$ : Interaction type (0 for competition, 1 for mutualism, 2 for predation)
$Bspe$ : Number of basal species needed for establishment
$Bext$ : Number of basal species needed for extinction
$Smax$ : Maximal number of species allowed during the simulation

$d$ : Shape of the exponential decay of the probability of establishment/speciation
$I_{max}$ : Maximal number of species in the community/system


 <!---
================================================================================
The rate of change in species richness $R$ is a dynamic balance between speciation events and extinction events, represent as follows:

$\frac{dR}{dt}\frac{1}{R} = S(R) - E(R)$

We consider a successful speciation event to be the combination of a mutation leading to speciation and the acquisition of traits that are ecological suitable (i.e. they provide preys, mutualists or minimize competition). In what follows we use exponential equations for these functions, but other functional forms could be considered as well, depending on the assumptions considered. Thus, we define:

$S(R) = u_{max}(u_0 + u_1 e^{-\alpha I}}$

and

$E(R) = e_{max}(e_0 + e_1 e^{-\alpha I})$

where $u_{max}$ and $e_{max}$ are the asymptotic speciation and extinction rates respectively, and $u_1-u_0$ and $e_1-e_0$ are the speciation and extinction rates in absence of interactions. There are multiple ways to parameterize those functions, for different types of interactions. These are summarized at Table 1 and the functions are illustrated at Fig. 1-3. In short, interactions are modifiers of the $u_{max}$ and $e_{max} and the shape of the function depends on the type of interactions. Maximal speciation rate could happen either at null diversity (e.g. in absence of interactions) or at infinite interactions (e.g. with mutualism).
--->
