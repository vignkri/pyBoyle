
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
# --
import boyle.preprocessing as p


# In[2]:


# AMINO is missing
COLUMN_ORDER = ["time", "temp", "Qin", "Qout",
                "carb_ins", "carb_ine", "carb_sol", "prot_ins",
                "prot_ine", "prot_sol", "lipi_sol", "ac_lcfa",
                "ac_prop", "ac_buty", "ac_vale", "ac_ace",
                "dg_nh4", "dg_ch4", "dg_co2", "dg_h2s",
                "io_z", "io_p", "io_a"]
ADDITIONAL_COLUMNS = ["DeadCells", "degrcarb", "degrAmino",
                      "degrLipids", "degrLCFA", "degrProp",
                      "degrButy", "degrValer", "degrAcet"]


# In[80]:


temperature = pd.read_csv("./load/data/from/inputs/reactor_temperatures.csv",
                          parse_dates=[0],
                          dtype={"T_r1": np.float, "T_r2": np.float,
                                 "T_r3": np.float})
temperature.median(axis=0)


# In[3]:


substrate = pd.read_csv("./load/data/from/inputs/substrate_database.csv",
                        index_col=[0, 1], header=[0, 1])
substrate.index = substrate.index.droplevel(1)
# substrates are in g/l
substrate.head(5)


# In[4]:


substrate_cols = [col[0] for col in substrate]
new_index = [substrate_cols,
             ["AVG"] * 20 + ["STDEV"] * 20]
substrate.columns = new_index


# In[5]:


feed = pd.read_csv("./load/data/from/client/inputs/reactor_feed_data.csv",
                   header=[0, 1, 2], index_col=[0], parse_dates=[0]) \
    .sort_index()
feed.sample(5)


# In[6]:


# -- Identify Common elements
feed_names = set(feed.columns.levels[2])
substrate_names = set(substrate.columns.levels[0])
intersection_columns = substrate_names.intersection(feed_names)


# In[7]:


# -- Create substrate dataset from the intersected column names
valid_columns = []
for name in intersection_columns:
    valid_columns.append((name, p.createSampledComposition(substrate[name])))
# --
data = []
for sub in valid_columns:
    substrate_name = sub[0]
    records = sub[1]
    record_data = {}
    for row in records:
        record_data.update(
            {"{}".format(row[0]): np.random.choice(row[1], 1)[0]}
        )
    # --
    record_data.update({"param": substrate_name})
    data.append(record_data)
# -- set up as pandas
sampledSubstrate = pd.DataFrame.from_dict(data).set_index("param").sort_index()


# In[8]:


# -- fix the feed units of the dataset
_feed = feed.fillna(0)
_feed_columns = _feed.columns.droplevel([0, 1])
_feed.columns = _feed_columns
idx = pd.IndexSlice
_feed = _feed[idx[list(intersection_columns)]].copy().sort_index(axis=1)
# -- correct solid columns
solids = [col for col in _feed.columns if col.startswith("SF")]
liquids = [col for col in _feed.columns if col.startswith("LF")]
# -- convert units
_feed.loc[:, solids] = _feed.loc[:, solids] * 1000
_feed.loc[:, liquids] = _feed.loc[:, liquids] * 1000


# In[15]:


# -- create composition matrix
composition = (_feed @ sampledSubstrate)
flow_sum = _feed.sum(axis=1)
composition = composition.truediv(flow_sum, axis=0)


# In[81]:


# ---
composition["Qin"] = flow_sum
composition["Qout"] = composition.Qin.values
# --
composition["temp"] = 38.115
composition["time"] = np.hstack(
    [np.array([720]), 720 + np.arange(1, composition.shape[0]) * 24]
)


# In[99]:


# -- reindex to contain all valid columns in the dataset
composition = composition.reindex(COLUMN_ORDER, axis=1)
for col in ADDITIONAL_COLUMNS:
    composition[col] = 0


# In[100]:


# -- set unusable values to zero including infinity values
composition.fillna(0, inplace=True)
composition[composition == np.inf] = 0
composition[composition == (-np.inf)] = 0


# In[101]:


# composition["BMP"] = 0
composition.columns


# In[102]:


np.save("./output/to/simulation/folder/feed.npy", composition)
