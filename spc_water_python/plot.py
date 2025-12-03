import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv('stats.dump')
print(df)

df['potential_energy'].plot()
plt.show()
