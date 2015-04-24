import matplotlib.pyplot as plt
import numpy as np
from get_data import df

df = df[np.isfinite(df.F1)]
df = df[np.isfinite(df.W10)]

ax = plt.subplot(111)
logF1 = np.log10(np.abs(df.F1))
logW10 = np.log10(df.W10)
ax.plot(logF1, logW10, ".")

plt.show()

