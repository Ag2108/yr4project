import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ab_phase=0.0
chartname='occ.png'

data = pd.read_csv('loctest.dat',delim_whitespace=True,header=None)

data_0 = data[np.isclose(data[0],ab_phase,atol=1e-7)]

vals=data_0.iloc[0,1:].values

sites=range(1,len(vals)+1)

plt.bar(sites,vals,color='skyblue')

plt.xticks(sites, [f'{i}' for i in sites])
plt.xlabel('States')
plt.ylabel('State Occupancy at Lowest Band')
#plt.title()

plt.savefig(f'{chartname}')
print(f'Plot successfully generated and saved as{chartname}')