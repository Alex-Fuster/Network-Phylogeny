# Methods

## Model structure

dynamic niche model
community level, where each population or species are represented by two traits. The first one is numerical trait value, that one can consider as niche optimum value. The second trait represent the interaction range of the species/population, the niche trait value is considered as the center of this range. For two species to interaction, the niche trait value of a species _i_ should be into the range of interaction of the species _j_. (**figure avec representation des valeurs de traits et des range d'interaction**)
start with one population, which is allow or not to diverse, depending on origination and extinction probabilities.


## Simulations

Table 1 : Parameters, signification and values. :warning: *$u_1-u_0$ and $e_1-e_0$ are the establishment and extinction rates in absence of interactions.*

| Parameter | Signification                                                                                                              |Value|
|-----------|----------------------------------------------------------------------------------------------------------------------------|-----|
|$int$      |Interaction type, 0 for competition, 1 for facilitation, 2 for predation                                                    |     |
|$S_{max$   |Maximal number of species in the system                                                                                     |1000 |
|$av_r$     |Half range of the niche of the first species                                                                                |0.1  |
|$sd$       |Standard deviation of the normal distribution used to calculate the niche optimum trait of a new species        |2*$av_r$ + 0.0001|
|$u_{max}$  |Speciation probability                                                                                                      |0.1  |
|$a_{eneg}$ |Shape of the exponential decay of the negative extinction - interaction relationship                                        |0.025|
|$e_{0neg}$ |Asymptotic extinction probability with infinite negative interactions                                                       |0.5  |
|$e_{1neg}$ |Asymptotic extinction probability with infinite negative interactions                                                       |1    |
|$a_{epos}$ |Shape of the exponential decay of the positive extinction - interaction relationship                                        |1.2  |
|$e_{0pos}$ |Asymptotic extinction probability with infinite positive interactions                                                       |0.075|
|$e_{1pos}$ |                                                                                                                            |5.19 |
|$a_{uneg}$ |Shape of the exponential decay of the colonization - interaction relationship                                               |0.075|
|$u_{0neg}$ |Asymptotic establishment probability with infinite competition interactions                                                 |0.075|
|$u_{1neg}$ |Establishment probability with absence of competition interaction                                                           |2    |
|$a_u$      |Shape of the exponential decay of the colonization - interaction relationship                                               |0.45 |
|$u_0$      |Asymptotic establishment probability with infinite facilitation interactions                                                |1    |
|$u_1$      |Establishment probability with absence of facilitation interactions                                                         |-1   |
|$d$        |Decrease speed of the establishment probability                                                                             |0.5  |
|$B_{spe}$  |Constant minimal number of interaction per species (establishment prob)                                                     |4    |
|$B_{ext}$  |Constant minimal number of interaction per species (extinction prob)                                                        |4    |
|$I_{max}$  |Maximal number of interacting species                                                                                       |40   |

## Analysis

### Phylogenetic struture
- alpha value
- sackin index
- branch length
- distance

### Network structure
- degree

### Traits evolution
- ?
