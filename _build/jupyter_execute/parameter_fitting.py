#!/usr/bin/env python
# coding: utf-8

# # Parameter Fitting for Mortals
# 
# Parameters are constants used in models, such as kinetics contants for reactions.
# Although there is often knowledge of likely ranges of parameter values, the appropriate value is often unknown
# for a specific model.
# 
# So, how do we proceed without knowing parameter values? The answer is that we estimate parameter values from
# observed data.
# Here is a simple example.
# Consider the irreversible reaction ``A + B`` $\rightarrow$ ``C`` with mass action kinetics.
# That is, the kinetics expression is $kAB$.
# We do not know $k$.
# At first glance, it seems that we cannot go further.
# However, suppose that we have observed data for ``C``.
# So, we can *guess* an initial value of $k$, run a simulation, and compare
# the simulated ``C`` with the observed ``C``.
# We then repeatedly change $k$ in the direction so that we get a small difference
# between observed ``C`` and simulated ``C``.
# 
# The following sections elaborate on this procedure using examples and providing code.

# In[1]:


# Imports 
import SBstoat


# ## Fitting Basics

# Some key concepts are:
# - A **parameter** is a constant whose value is unknown.
# - **Observed data** are data that are used for fitting parameter values.
# - **Simulated data** are data produced by a simulation for an assignment of values to parameters.
# - **Residuals** are the difference between observed data and simulated data.

# The basic algorithm for fitting is described below.
# Below, we outline the approach for a single parameter.
# 
# 1. Pick an initial value for $k$ and simulate the reaction.
# 2. Calculate the **residuals**, the difference between simulated ``C`` and observed ``C``. Increase $k$.
# 3. Obtain a new simulated ``C`` and calculate new residuals.
# 4. If the new residuals are very small (e.g., very small variane), return $k$.
# 5. If the new residuals are smaller than the last residuals (e.g., smaller variance), change $k$ in the same direction.
# 6. If the new residuals are larger than the last residuals, change $k$ in the opposite direction (with a smaller amount of change).
# 7. Go to Step 3.

# ## Fitting With SBstoat

# ### Running Example

# ### Doing the fit

# ### Evaluating Quality of Fit

# ## A More Realistic Example

# ## Estimating Parameter Variances With Bootstrapping

# ## Tips for Fitting

# There is considerable art in fitting. The quality of parameter fits
# can be improved by exploring the following:
# - **What you fit**. This means, which species concentrations you include in the residuals, what part of the time course
# you try to fit, and the ranges you choose for parameter values.
# - **How you fit**. There are several different algorithms for doing non-linear minimizations...
# - **How much compute**. Fitting is computationally intensive. More compute time is often needed for a better result. ``SBstoat`` helps
# in this regard by making use of multiple cores.
