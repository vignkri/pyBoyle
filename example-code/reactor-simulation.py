
# coding: utf-8

# In[1]:


import boyle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# # Reactor 1 Simulation

# In[2]:


path_data = "./path/to/client/project/data/folder/"
ph_settings = {"method": "fixed", "value": 7}


# In[3]:


_meta_ = {"name": "Reactor-01 Mesophilic Value",
          "data": path_data,
          "description": "Testing mesophilic constants",
          "tags": "Testing, Reactor01-Agrana"}
_settn_ = {"process": "debug",
           "constants": "standard",
           "step_size": 1,
           "ph": ph_settings,
           "solver": {"method": "bdf", "order": 1,
                      "nsteps": 100000, "relative": 1e-4,
                      "absolute": 1e-8}}


# ## 1. Simulate Results

# In[4]:


_config = {"metadata": _meta_, "settings": _settn_}
simulator = boyle.Manager(_config)


# In[5]:


dataset = simulator.start()


# ## Column Names

# In[6]:


HEADER_START = ["run_no", "time"]
HEADER_DEBUG = ["mu_1", "mu_2", "mu_3", "mu_4", "mu_5", "mu_6",
                "mu_7", "mu_8", "pH", "raw_flow", "T_flow_in"]
HEADER_VOL = ["volume"]
HEADER_CORE = ["carb_ins", "carb_ine", "carb_sol", "prot_ins",
               "prot_ine", "amino", "lipids", "ac_lcfa", "ac_prop",
               "ac_buty", "ac_val", "ac_ace", "dg_nh4",
               "dg_ch4", "dg_co2", "dg_h2s", "io_z", "io_p", "io_a"]
HEADER_DEGRADERS = ["dead_cell", "degr_carb", "degr_amino", "degr_lipid",
                    "degr_lcfa", "degr_hprop", "degr_butyr", "degr_valer",
                    "degr_acet", "gf_nh3", "gf_ch4", "gf_co2", "gf_h2s"]
HEADER_END = ["gasrate"]
# -- Debug Headers
new_debug_headers = ["run_no", "time",
                     "mu_1", "mu_2", "mu_3", "mu_4", "mu_5",
                     "mu_6", "mu_7", "mu_8",
                     "pH", "T_flow_in", "flow_in",
                     "carb_is", "carb_in", "carb",
                     "prot_is", "prot_in", "amino",
                     "lipids", "lcfa", "hpr", "hbut",
                     "hval", "hac", "ch4", "co2", "h2s", "z+"]


# In[7]:


debug = pd.DataFrame(dataset.debug[1:], columns=new_debug_headers)
debug.to_csv("./Agrana/Exports/Data/debug_dataset_r1_bmp.csv")


# In[7]:


y_hat = pd.DataFrame(np.asarray(dataset.y_hat))
y_hat.columns = HEADER_START + HEADER_VOL + HEADER_CORE + HEADER_DEGRADERS


# In[9]:


interpolate = pd.DataFrame(
    boyle.tools.interpolateData(np.asarray(dataset.y_hat), time_step=24))
interpolate.columns = HEADER_START + HEADER_VOL + \
    HEADER_CORE + HEADER_DEGRADERS


# In[16]:


plt.figure(figsize=(18, 3))
plt.plot(((y_hat.gf_ch4.shift(-1) - y_hat.gf_ch4) * 24) / 1000)
plt.show()
