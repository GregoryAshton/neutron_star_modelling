import urllib2
from BeautifulSoup import BeautifulSoup
import matplotlib.pyplot as plt
import numpy as np
from nsmod import Plot

url = "http://www.jb.man.ac.uk/pulsar/glitches/original.html"
soup = BeautifulSoup(urllib2.urlopen(url).read())

table = soup.find("table")

names = []
dFs = []
dFs_err = []
for row in table.findAll("tr")[5:-3]:
    cols = row.findAll('td')

    n = cols[1].getText()
    dF = cols[6].getText()
    dF_err = cols[7].getText()
    if dF_err == "*":
        dF_err = np.nan
    if float(dF) != 0:
        names.append(n)
        dFs.append(float(dF) * 1e-6)
        dFs_err.append(float(dF_err) * 1e-6)

dFs = np.array(dFs)
dFs_err = np.array(dFs_err)
names = np.array(names)

print "Using {} data points".format(len(names))

log10dF = np.log10(dFs)
fig, ax = plt.subplots()
out = ax.hist(log10dF, bins=13, histtype="step", color="k")
ax.set_xlabel("$\log_{10}{\Delta f}$", size=22)
ax.set_ylabel("Count", size=22)
plt.savefig("EspinozaGlitchCat_dF.pdf")

print "Maximum glitch size = {}".format(dFs.max())
