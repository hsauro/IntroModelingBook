#!/usr/bin/env python
# coding: utf-8

# # Chapter X: Building a mechanistic model with the Antimony language 
# 
# <br>
# <div align='center'>
#     <figure>
#         <img src="https://raw.githubusercontent.com/vporubsky/tellurium-libroadrunner-tutorial/master/data_aggregation_logo.png" width = "50%" style="padding: 0px">
#         <figcaption>Model construction requires several stages including data aggregation, annotation and creation of a mathematical representation of the system which can be represented in a computational environment.</figcaption>
#     </figure>
# </div>
# <br>
# <br>
# 
# <span style="color:red">Include data aggregation, writing reactions, adding species and setting values</span>.
# 

# ## Data aggregation
# 
# Data aggregation is the collection of data from multiple experiments, scientific papers, and online data sources. Typically you will need to curate your aggregated data to ensure the quality of measurements are suffiencient to include in your model.
# 
# <ul>
#   <li>SABIO-RK: biochemical reaction kinetics database</li>
#      <ul class="square">
#       <li>Describes chemical reactions and kinetics</li>
#       <li>Contains information about participants and modifiers in reactions</li>
#       <li>Metabolic and signaling network reactions</li>
#      </ul>
#   <li>BRENDA: the comprehensive enzyme information system</li>
#      <ul class="square">
#       <li>Enzyme information classified by the biochemical reaction it catalyzes</li>
#       <li>Kinetic information about substrates and products is available</li> 
#      </ul>
#   <li>ChEBI: dictionary of "small" chemical compounds</li>
#   <li>KEGG: collection of pathway/genome/diesease/drug databases</li>
#   <li>BioCYC: collection of pathway/genome databases</li>
#        <ul class="square">
#       <li>Search for genes, proteins, metabolites or pathways, and the occurence of your term will be located in multiple databases</li> 
#      </ul>
#   <li>BioModel: repository of mathematical models of biological systems</li>
#       <ul class="square">
#       <li> *Will be covered in more detail later in the course</li> 
#      </ul>
# </ul>
# 
# <a href="https://www.sciencedirect.com/science/article/abs/pii/S0958166917301428?via%3Dihub">Appendix A of Goldberg et al. (2018)</a> provides a useful and more comprehensive list of data sources containing intracellular biochemical data. 
# 
# 
# <ul>
#   <li>Metadata: data that describes biochemical data</li>
#   <li>Collect information about:</li>
#      <ul class="square">
#       <li>Units</li>
#       <li>Estimates of measurement accuracy</li>
#       <li>Annotations</li>
#       <li>Ontology terms defining the annotations</li>
#       <li>etc.</li>
#      </ul>
#   <li>Collect provenance data:</li>
#      <ul class="square">
#       <li>Lab which generated the data</li>
#       <li>Experimental conditions</li>
#       <li>Protocol used to generate the data</li>
#       <li>Paper which reported the measurement</li>
#       <li>etc.</li>

# In[1]:


import requests
import io
import pandas as pd

QUERY_URL = 'http://sabiork.h-its.org/sabioRestWebServices/kineticlawsExportTsv'

# Specify search fields and search terms
query_dict = {"Organism":'"SARS coronavirus"', "Product":'"Spike glycoprotein RBD-SD1-ACE2 complex"'}
query_string = ' AND '.join(['%s:%s' % (k,v) for k,v in query_dict.items()])

# Specify output fields and send request
query = {'fields[]':['EntryID', 'Organism', 'UniprotID','ECNumber', 'Parameter'], 'q':query_string}
request = requests.post(QUERY_URL, params = query)
request.raise_for_status()

# Print query results
print('\nThe results of SABIO-RK query for the formation of spike glycoprotein RBD-SD1-ACE2 complex in SARS coronavirus:\n')
print(request.text)
data_string = request.text
data = io.StringIO(data_string)
df = pd.read_csv(data, sep="\t")

parameter_types = list(df['parameter.type'])
if 'rate const.' in parameter_types:
    k7 = list(df['parameter.startValue'])[parameter_types.index('rate const.')]
    k7_units = list(df['parameter.unit'])[parameter_types.index('rate const.')]

print(f"\nThe formation of the spike glycoprotein RBD-SD1-ACE2 complex occurs at a rate of {k7:.3f} {k7_units}.")


# ## Modeling biochemical reaction kinetics
# 
# <ul class="roman">
#  <li>Speed of the biochemical reaction involved in the network determines evolution of the system</li>
#   <li>Factors that influence the rates of reactions must be considered and modeled appropriately</li>
#   <ul class="square">
#   <li>e.g. the presence of an activator or repressor would change the baseline speed of that reaction</li>
#  </ul>
#  <li>Consider the mechanisms involved and determine appropriate simplifications</li>
# </ul>
# 
# 
# <ul class="roman">
#  <li>Models which describe how variables in a system evolve over time</li>
#  <ul class="square">
#   <li>e.g. floating species concentrations</li>
#  </ul>
#  <li>Quantities are derived from the variables</li>
#   <ul class="square">
#   <li>e.g. pathway flux</li>
#  </ul>
#  <li>Some parameters of the model are fixed by the modeler</li>
#    <ul class="square">
#   <li>e.g. rate constants</li>
#   <li>e.g. enzyme concentrations</li>
#   <li>e.g. boundary species concentrations</li>
#  </ul>
#  <li>Deterministic differential equations models are useful when we can assume there are a large number of participants in the chemical reactions</li>
#  <li>Stochastic models are useful for dilute systems in which reactions may not occur at every timepoint</li>
# </ul>
# 
# 
# <br>
# <br>
# <br>
# <div align='center'><img src="https://raw.githubusercontent.com/vporubsky/tellurium-libroadrunner-tutorial/master/network_to_diff_eq.png" width = "60%" style="padding: 0px"></div>

# ## Writing a simple Antimony string 
# <br>
# Consider the system depicted below, in which species S1 is converted into species S2 with some chemical reaction rate, k1. Assuming this model follows mass-action kinetics as described in Chapter X, the rate law is shown: k1*[S1]. 
# <br>
# <div align='center'><img src="https://raw.githubusercontent.com/vporubsky/tellurium-libroadrunner-tutorial/master/simple-antimony.png" width="75%" style="padding: 20px"></div>

# In[2]:


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


# ## Editing a complex model of a true biological system
# <br>
# It is common to build upon existing models to enable new analyses...
# <br>
# <div align='center'><img src="" width="75%" style="padding: 20px"></div>

# In[3]:


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
