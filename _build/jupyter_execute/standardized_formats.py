#!/usr/bin/env python
# coding: utf-8

# # Chapter X: SBML, SED-ML, and the COMBINE archive: standard formats supported in Tellurium 

# ## Import and export capabilities with Tellurium
# <br>
# Models can be imported from the BioModels Database, given the appropriate BioModel ID. However, due to a firewall on nanoHUB, we are unable to perform this import. Instead, upload BIOMD0000000090.xml to your Tellurium tool on nanoHUB.
# 
# This is a model of respiratory oscillations in Saccharomyces cerevisae by <a href="https://www.ebi.ac.uk/biomodels/BIOMD0000000090">Jana Wolf et al. (2001):</a> </div>
# <br>
# 
# <div align='center'><img src="https://raw.githubusercontent.com/vporubsky/tellurium-libroadrunner-tutorial/master/wolf_publication.PNG" width="65%" style="padding: 20px"></div>
# <br>
# <div align='center'><img src="https://raw.githubusercontent.com/vporubsky/tellurium-libroadrunner-tutorial/master/wolf_network.PNG" width="65%" style="padding: 20px"></div>

# In[1]:


#%% You can load an SBML model directly from the BioModels Database, given the BioModel ID
# wolf = te.loadSBMLModel("http://www.ebi.ac.uk/biomodels-main/download?mid=BIOMD0000000090")
wolf = te.loadSBMLModel("BIOMD0000000090.xml")
wolf.simulate(0, 100, 1000, ['time', 'oxy']) # note that specific species can be selected for recording concentrations over the timecourse
wolf.plot(figsize = (10, 6), xtitle = 'Time', ytitle = 'Concentration')
 


# In[ ]:


# Export the model you just accessed from BioModels to the current directory as an Antimony string
wolf.reset()
print(wolf.getFloatingSpeciesConcentrationIds())
wolf.exportToAntimony('wolf.txt', current = True) 


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
