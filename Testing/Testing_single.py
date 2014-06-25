from nsmod.one_component_model import main
from nsmod import Plot
import matplotlib.pyplot as plt

file_name = main()

Plot.Spherical_Plot(file_name)

plt.show()
