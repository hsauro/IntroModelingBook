#!/usr/bin/env python
# coding: utf-8

# # Kinetics in a Nutshell
# 
# 
# Reaction kineticsis the study of how fast chemical reactions take place, what factors influence the rate of reaction, and what mechanisms are responsible.
# 
# A generalized reaction scheme is shown below:
# %
# $$ a \text{A} + b \text{B} + \ldots \rightarrow + p \text{P} + q \text{Q} + \ldots $$
# 
# where $a, b, \ldots, p, q, \ldots$ are called the stoichiometric amounts.
# 
# ## Defintions
# 
# **Rate of Change**
# 
# The rate of change is defined as the rate of change in concentration or amount of a designated molecular species.
# 
# $$ \text{Rate of Change} = \frac{dA}{dt} $$
# 
# **Stoichiometric Coefficients**
# 
# The {\bfseries stoichiometric coefficient}, $c_i$, for a molecular species A$_j$, is the difference between the molar amount of species, $i$ (or {\bfseries stoichiometric amount}) on the product side and the molar amount of the same species on the reactant side.
# 
# $$ c_i = \text{Molar Amount of Product} - \text{Molar Amount of Reactant} $$
# 
# In the reaction, 2 A $\longrightarrow$ B, the molar amount of A on the product side is zero while on the reactant size it is two. Therefore the stoichiometric coefficient of A is given by $0-2 = -2$. In many cases a particular species will only occur on the reactant or product side. As a result, reactant stoichiometric coefficients tend to be {\em negative} while product stoichiometric coefficients tend to be {\em positive}.
# 
# **Reaction Rates**
# 
# The reaction rate, often denoted by the symbol $v$, is measured with respect to a given molecular species normalized by the species stoichiometric coefficient. This definition ensures that no matter which molecular species in a reaction is measured, the reaction rate is uniquely defined for that reaction. More formally, the reaction rate for the reaction:
# 
# $$ a \text{A} + b \text{B} + \ldots \rightarrow p \text{P} + q \text{Q} + \ldots $$
# 
# is:
# 
# $$ v  = \frac{1}{c_a} \frac{d\!A}{dt} = \frac{1}{c_b}
# \frac{d\!B}{dt} \ldots= \frac{1}{c_p} \frac{d\!P}{dt} =
# \frac{1}{c_q} \frac{d\!Q}{dt} \ldots $$
# 
# where $c_x$ are the stoichiometric coefficients. Alternatively, we can express the rate of change in terms of the reaction rate as:
# 
# $$ \frac{d\!A}{dt} = c_a v $$
# 
# 
# ## Elementary Mass-Action Kinetics
# 
# An elementary reaction is one that cannot be broken down into simpler reactions. Such reactions will often display simple kinetics called mass-action kinetics. For a reaction of the form:
# 
# $$ a \text{A} + b \text{B} + \ldots \rightleftharpoons p \text{P} + q \text{Q} + \ldots $$
# 
# the mass-action kinetic rate law is given by:
# 
# $$ v = k_1 A^a B^b \ldots - k_2 P^p Q^q \ldots $$
# 
# $k_1$ and $k_2$ are the forward and reverse rate constants, respectively.
# 
# # Enzyme Kinetics
# 
# Enzymes are protein molecules that can accelerate a chemical reaction without changing the reaction equilibrium constant.
# 
# ## Michaelis-Menten Kinetics
# 
# The standard model for enzyme action\index{enzyme action} describes the binding of free enzyme to the reactant forming an enzyme-reactant complex. This complex undergoes a transformation, releasing product and free enzyme. The free enzyme is then available for another round of binding to new reactant.
# 
# $$ \text{E} + \text{S} \rightleftharpoons \text{ES} \stackrel{k_2}{\longrightarrow} \text{E} + \text{P} $$
# 
# 
# where $k_1, k_{-1}$ and $k_2$ are rate constants, $S$ is substrate, $P$ is product, $E$ is the free enzyme, and $\mathit{ES}$ the enzyme-substrate complex.
# 
# By assuming a steady state condition on the enzyme substrate complex, we can derive the Briggs-Haldane equation\index{Briggs-Haldane} relation (sometimes mistakenly called the Michaelis-Menten equation):
# 
# $$ v = \frac{Vm\ S}{K_m + S} $$
# 
# where $V_m$ is the maximal velocity, and $K_m$ the substrate concentration that yields half the maximum velocity.

# In[1]:


a = 2
print('my 1st line')
print(f'my {a}nd line')


# ## Reversible Rate laws
# 
# An alternative and more realistic model is the reversible form:
# 
# $$ \text{E} + \text{S} \rightleftharpoons \text{ES} \rightleftharpoons \text{E} + \text{P} $$
# 
# The aggregate rate law for the reversible form of the mechanism can also be derived and is given by:
# 
# $$ v = \frac{V_f\ S/K_S - V_r\ P/K_P}{1 + S/K_S + P/K_P} $$
# 
# ## Haldane Relationship
# 
# For the reversible enzyme kinetic law there is an important relationship:
# 
# $$ K_{eq} = \frac{P_{eq}}{S_{eq}} = \frac{V_f\ K_P}{V_r\ K_S} $$
# 
# This equation sh\-ows that the four kinetic constants, $V_f, V_r, K_P$ and $K_S$ are not independent. Haldane relationships can be used to eliminate one of the kinetic constants by substituting the equilibrium constant in its place. This is useful because equilibrium constants tend to be known compared to kinetic constants which may be unknown. By incorporating the Haldane relationship, we can eliminate the reverse maximal velocity ($V_r$) from~\ref{eqn:revMM} to yield the equation:
# 
# $$ v = \frac{V_f/K_S (S - P/K_{eq})}{1 + S/K_S + P/K_P} $$
# 
# Separating out the terms makes it easier to see that the above equation can be partitioned into a number of distinct parts:
# 
# $$ v = V_f\ \cdot\ (1 - \Gamma/K_{eq})\ \cdot\ \frac{S/K_s}{1 + S/K_S + P/K_P} $$
# 
# where $\Gamma = P/S$. The first term, $V_f$, is the maximal velocity; the second term, $(1 - \Gamma/K_{eq})$, indicates the direction of the reaction according to thermodynamic considerations. The last term refers to the fractional saturation with respect to substrate. Thus we have a maximal velocity, a thermodynamic and a saturation term.
# 
# ## Competitive Inhibition
# 
# There are many molecules capable of slowing down or speeding up the rate of enzyme catalyzed reactions. Such molecules are called enzyme inhibitors and activators. One common type of inhibition, called competitive inhibition, occurs when the inhibitor is structurally similar to the substrate so that it competes for the active site by forming a dead-end complex.
# 
# ![InhibitorMechanisms50.png](attachment:InhibitorMechanisms50.png)
# 
# 
# Competitive and Uncompetitive Inhibition. $P$ is the concentration of product, $\mathit{E}$ is the free enzyme, $\mathit{ES}$ the enzyme-substrate complex, and $\mathit{ESI}$ the enzyme-substrate-inhibitor complex. }
# 
# 
# The kinetic mechanism for a pure competitive inhibitor is shown in Figure above (a) where $\mathit{I}$ is the inhibitor and $\mathit{EI}$ the enzyme inhibitor complex. If the substrate concentration is increased, it is possible for the substrate to eventually out compete the inhibitor. For this reason the inhibitor alters the enzyme's apparent $K_m$, but not the $V_{\text{m}}$.
# 
# $$
#  v = \frac{V_m\ S}{S + K_m\left(1 + \displaystyle\frac{I}{K_i}\right)} \\[8pt]
#    = \frac{V_m\ S/K_m}{1 + S/K_m + I/K_i}
# $$
# 
# At $I=0$, the competitive inhibition equation reduces to the normal irreversible Michaelis-Menten equation. Note that the term $K_m (1 + I/K_i)$ in the first equation more clearly shows the impact of the inhibitor, $I$, on the $K_m$. The inhibitor has no effect on the $V_m$.
# 
# A reversible form of the competitive rate law can also be derived:
# 
# $$ v = \frac{\displaystyle\frac{V_m}{K_{s}} \left( \displaystyle S - \displaystyle\frac{P}{K_{\text{eq}}}\right)}{1 + \displaystyle\frac{S}{K_{s}} + \displaystyle\frac{P}{K_{p}} + \displaystyle\frac{I}{K_i}}
# $$
# 
# where $V_m$ is the forward maximal velocity, and $K_{s}$ and $K_{p}$ are the substrate and product half saturation constants.
# 
# Sometimes reactions appear irreversible, where no discernable reverse rate is detected, and yet the forward reaction is influenced by the accumulation of product. This effect is caused by the product competing with substrate for binding to the active site and is often called product inhibition. Given that product inhibition is a type of competitive inhibition, we will briefly discuss it. An important industrial example of this is the conversion of lactose to galactose by the enzyme $\beta-$galactosidase where galactose competes with lactose, slowing the forward rate.
# 
# To describe simple product inhibition with an irreversible reaction, we can set the $P/K_{eq}$ term in the reversible Michaelis-Menten rate law~\eqref{eqn:revMM} to zero. This yields:
# 
# $$ v = \frac{V_{m} S}{S + K_m \left( 1 + \displaystyle\frac{P}{K_p}\right)} $$
# 
# It is not surprising to discover that equation has exactly the same form as the equation for competitive inhibition. As the product increases, it out competes the substrate and therefore slows down the reaction rate.
# 
# We can also derive the equation by using the following mechanism and the rapid-equilibrium assumption:
# 
# $$ \text{E} + S \rightleftharpoons \text{ES} \longrightarrow \text{EP} \rightleftharpoons  \text{E} + \text{P} $$
# 
# where the reaction rate $v$ is assumed to be proportional to $\mathit{ES}$.
# 
# 
# ![InhibitorMechanisms50.png](attachment:competitiveInhibitionGraph60.png)
# 
# 
# Competitive inhibition: Effect of substrate concentration on the reaction rate as a function of different inhibitor concentrations: $V_m = 1, K_m = 1, K_i = 1$. Solid circles mark the apparent $K_m$ values and the light circle the $K_m$ in the absence of inhibitor. As the inhibitor concentration increases the apparent $K_m$ increases in value.
# 
# CODE EXPERIMENTS WITH COMPETITIVE INHIBITION
