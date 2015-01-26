import matplotlib.pyplot as plt 
import GetLyneData
import TNtools as TN

TN.PlotDefaults()

x, y = GetLyneData.GetData(D=None,
file_name="/home/greg/timing-noise/ExtractDataFromLyne/B1828_W10_01.txt")


fig, (ax) = plt.subplots(ncols=1, figsize=(6, 3))
ax.plot(x, y, "o-", markersize=3)
ax.set_xlabel("Modifed Julian Date [days]", size=14)
ax.set_ylabel("$W_{10}$ [ms]", size=14)
ax.set_xlim(53000, 55000)
ax.set_yticks(ax.get_yticks()[:-1])
plt.tight_layout()

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

plt.savefig("/home/greg/Neutron_star_modelling/Posters/OxfordPoster/img/raw_data.pdf")
plt.show()
