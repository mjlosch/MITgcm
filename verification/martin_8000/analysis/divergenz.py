#Analysis using XMITGCM

from xmitgcm import open_mdsdataset
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

data_dir = "../run_2"
ds = open_mdsdataset(data_dir)

print("Es wurden ", ds.time.size, "Zeitschritte gemacht von ", pd.to_timedelta(ds.time[0].values), "zu ", pd.to_timedelta(ds.time[-1].values))
index = int(input("Plots nach wie vielen Zeitschritten?"))
#print(ds)
print(ds.SIuice[index].values)