#!/usr/bin/env python
# coding: utf-8

# # Chapter X: Running deterministic and stochastic simulations with Tellurium and libRoadRunner 

# <br>
# <div align='center'>
#     <figure>
#     <img src="https://raw.githubusercontent.com/vporubsky/tellurium-libroadrunner-tutorial/master/tellurium-utilities.png" width="75%" style="padding: 20px">
#         <figcaption>Tellurium is an integrated modeling environment which makes available libraries to support biochemical model building and simulation. Reproduced from [cite Tellurium paper/docs].</figcaption>
#     </figure>
# </div>

# ## How does numerical simulation help us model network kinetics?
# 
# <ul>
#   <li>Provides a method to approximate analytical solutions for complex (often non-linear) systems</li>
#   <li>Kinetic laws describe the rates of change of species in the system, which can be modeled mathematically</li>
#   <li>For a sufficiently large network, simulators must be efficient and perform rapid numerical integration</li>
# </ul>
# 
# ## How is numerical integration performed?
# 
# <ul>
#   <li>Requires differential equations describing the model, constant values, and initial conditions for all variables which are evolving in time</li>
#   <li>Numerical integration algorithms assume local linearity, and use tangent lines over small step sizes to gradually approximate the solution, starting from some initial value</li>
#   <li>Euler's method can be examined as a simple example</li>
# </ul>

# In[1]:


# First, install Tellurium, which comes with libRoadRunner 
get_ipython().system('pip install tellurium')


# In[2]:


import tellurium as te # Python-based modeling environment for kinetic models
import roadrunner as rr # High-performance simulation and analysis library
import numpy as np # Scientific computing package
import random # Generate random numbers 
import matplotlib.pylab as plt # Additional Python plotting utilities


# ## Simulating the simple Antimony model

# In[3]:


Ant_str = """
model test               # name the model
    compartment C1;      # specify compartments
    C1 = 1.0;            # assign compartment volume
    species S1, S2;      # specify species
    S1 in C1; S2 in C1;  # allocate species to appropriate compartment

 
    J1: S1 -> S2; k1*S1; # reaction; reaction rate law;
    
    S1 = 10.0;           # assign species initial conditions
    S2 = 0.0;

    k1 = 1.0;            # assign constant values to global parameters
end
"""

r = te.loada(Ant_str)  # create an executable model by loading the string to a RoadRunner object instance
result = r.simulate(0, 10, 25) # simulate(time_start, time_end, number_of_points)
print(result)


# ## Plotting the simple Antimony model

# In[4]:


r.plot(title = 'Uni-uni mass-action model', xtitle = 'Time', ytitle = 'Concentration', figsize = (8, 6)) 


# ## Writing a simulation with events

# In[5]:


Ant_str = """
model test               # name the model
    J1: S1 -> S2; k1*S1; # reaction; reaction rate law;
    
    S1 = 10.0;           # assign species initial conditions
    S2 = 0.0;

    k1 = 1.0;            # assign constant values to global parameters
    E1: at (time > 5): S1 = 10; # add an event - spike in S1
end
"""
r = te.loada(Ant_str)  # create an executable model by loading the string to a RoadRunner object instance
r.simulate(0, 10, 100) # simulate(time_start, time_end, number_of_points)
r.plot(title = 'Uni-uni mass-action model with event', xtitle = 'Time', ytitle = 'Concentration', figsize = (8, 6)) 


# ## Simulating a complex model with interesting dynamics

# In[6]:


repressilator_str = """
# Compartments and Species:
species M1, P3, P1, M2, P2, M3;

# Reactions:
J0:  -> M1; a_m1*(Kr_P3^n1/(Kr_P3^n1 + P3^n1)) + leak1;
J1: M1 -> ; d_m1*M1;
J2:  -> P1; a_p1*M1;
J3: P1 -> ; d_p1*P1;
J4:  -> M2; a_m2*(Kr_P1^n2/(Kr_P1^n2 + P1^n2)) + leak2;
J5: M2 -> ; d_m2*M2;
J6:  -> P2; a_p2*M2;
J7: P2 -> ; d_p2*P2;
J8:  -> M3; a_m3*(Kr_P2^n3/(Kr_P2^n3 + P2^n3)) + leak3;
J9: M3 -> ; d_m3*M3;
J10:  -> P3; a_p3*M3;
J11: P3 -> ; d_p3*P3;

# Species initializations:
M1 = 0.604016261711246;
P3 = 1.10433330559171;
P1 = 7.94746428021418;
M2 = 2.16464969760648;
P2 = 3.55413750091507;
M3 = 2.20471854765531;

# Variable initializations:
a_m1 = 1.13504504342841;
Kr_P3 = 0.537411795656332;
n1 = 7.75907326833983;
leak1 = 2.59839004225795e-07;
d_m1 = 0.360168301619141;
a_p1 = 5.91755684808254;
d_p1 = 1.11075218613419;
a_m2 = 2.57306185467814;
Kr_P1 = 0.190085253528206;
n2 = 6.89140262856765;
leak2 = 1.51282707494481e-06;
d_m2 = 1.05773721506759;
a_p2 = 8.35628834784826;
d_p2 = 0.520562081730298;
a_m3 = 0.417889543691157;
Kr_P2 = 2.71031378955001;
n3 = 0.44365980532785;
leak3 = 3.63586125130783e-11;
d_m3 = 0.805873530762994;
a_p3 = 4.61276807677109;
d_p3 = 1.54954108126666;

"""

repressilator = te.loada(repressilator_str)
repressilator.simulate(0, 100, 500)
repressilator.plot(figsize = (10, 8), linewidth = 3)


# ## Simulating a complex model of a true biological system
# 
# Using the model describe in Chapter X, we will demonstrate both deterministic and stochastic simulation studies.

# In[7]:


antimony_str = """
model sars_cov2_infection

// Define model equations
production_genomic_ssRNA:               -> gen_ssRNA;   k1*rep_ssRNA
transcription_genomic_ssRNA:  gen_ssRNA -> rep_ssRNA;   k2*gen_ssRNA
translation_of_viral_envelope:          -> envelope;    k3*rep_ssRNA
degradation_of_ssRNA:         rep_ssRNA -> ;            k4*rep_ssRNA
degradation_of_envelope:       envelope -> ;            k5*envelope
formation_of_free_virus_complex: gen_ssRNA + envelope -> virus_free;  k6*gen_ssRNA*envelope
binding_of_virus_complex_to_host:  virus_free + $ACE2 -> virus_bound; k7*virus_free*ACE2

// Annotate reactions
transcription_genomic_ssRNA.sboTerm = SBO:0000183
translation_of_viral_envelope.sboTerm = SBO:0000184

// Set global constants
k1 = 0.968; k2 = 0.025;
k3 = 1000; k4 = 0.222
k5 = 2; k6 = 7.5E-6
k7 = 0.112 

// Annotate kinetic constant distribution information
# k1.confidenceInterval = {0.883, 0.892}
# k4.confidenceInterval = {0.184, 0.186}

// Set initial conditions
rep_ssRNA = 1; gen_ssRNA = 0;
envelope = 0; ACE2 = 100;
virus_free = 0; virus_bound = 0;

end
"""


# In[8]:


import tellurium as te
import matplotlib.pyplot as plt
plt.style.use('ggplot')

# Load the model string into a RoadRunner object instance
model = te.loada(antimony_str)

# Generate deterministic results using CVODE
det_results = model.simulate (0, 200, 100, ['time', 'gen_ssRNA'])

# Set a seed value for reproducible stochastic output
model.seed = 124

plt.figure(2)
for i in range(10):
    # Reset variable concentrations to initial conditions
    model.reset()
    # Generate stochastic results using Gillespie's algorithm
    stoch_results = model.gillespie (0, 200, 100, ['time', 'gen_ssRNA'])
    # Plot stochastic simulation trajectory
    plt.plot(stoch_results['time'], stoch_results['gen_ssRNA'], linewidth=4, alpha=0.4)

# Plot deterministic results
plt.plot(det_results['time'], det_results['gen_ssRNA'], color='black', linewidth=2,)
plt.xlabel('Time')
plt.ylabel('[gen_ssRNA]')

# Save figure for Docker implementation and show
# plt.savefig('curated_k7_sars_cov2_infection_simulation.jpg', dpi = 300)
plt.show()

print(f'\nTimeseries of gen_ssRNA using fitted values for k1 and k4, and curated data for k7:')
print(det_results)


# <br>
# <html>
#    <head>
#       <title>Bibliography</title>
#    </head>
#    <body>
#       <h1>Bibliography</h1>
#       <ol>
#          <li>
#             <p>K. Choi et al., <cite>Tellurium: An extensible python-based modeling environment for systems and synthetic biology</cite>, Biosystems, vol. 171, pp. 74–79, Sep. 2018.</p>
#          </li>
#          <li>
#             <p>E. T. Somogyi et al., <cite>libRoadRunner: a high performance SBML simulation and analysis library.,</cite>, Bioinformatics, vol. 31, no. 20, pp. 3315–21, Oct. 2015.</p>         
#           <li>
#             <p>L. P. Smith, F. T. Bergmann, D. Chandran, and H. M. Sauro, <cite>Antimony: a modular model definition language</cite>, Bioinformatics, vol. 25, no. 18, pp. 2452–2454, Sep. 2009.</p>
#          </li>
#          <li>
#             <p>K. Choi, L. P. Smith, J. K. Medley, and H. M. Sauro, <cite>phraSED-ML: a paraphrased, human-readable adaptation of SED-ML</cite>, J. Bioinform. Comput. Biol., vol. 14, no. 06, Dec. 2016.</p>
#          </li>
#          <li>
#             <p>J. Wolf, H. Y. Sohn, R. Heinrich, and H. Kuriyama, <cite>Mathematical analysis of a mechanism for autonomous metabolic oscillations in continuous culture of Saccharomyces cerevisiae</cite>, FEBS Lett., vol. 499, no. 3, pp. 230–234, Jun. 2001.</p>
#          </li>
#       </ol>
#    </body>
# </html>
