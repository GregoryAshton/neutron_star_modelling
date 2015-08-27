import urllib2
from BeautifulSoup import BeautifulSoup
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from nsmod import Plot
import os


def GetData(download=False):
    pickled_data = "Espinoza_saved_data.csv"

    if download is False and os.path.isfile(pickled_data):
        print("Using saved data")
        df = pd.read_csv(pickled_data)
        return df
    else:
        print("Downloading the data")

    url = "http://www.jb.man.ac.uk/pulsar/glitches/original.html"
    soup = BeautifulSoup(urllib2.urlopen(url).read())

    table = soup.find("table")

    ddict = {'names': [],
             'dF': [],
             'dF_err': [],
             'dF_F': [],
             'dF_F_err': [],
             'dF1': [],
             'dF1_err': [],
             'dF1_F1': [],
             'dF1_F1_err': [],
             }
    for row in table.findAll("tr")[5:-3]:
        cols = row.findAll('td')

        idxs = [1] + range(6, 14)
        data = np.array([cols[i].getText() for i in idxs])
        data[data == "*"] = np.nan

        n, dF, dF_err, dF_F, dF_F_err, dF1, dF1_err, dF1_F1, dF1_F1_err = data
        ddict['names'].append(n)
        ddict['dF'].append(float(dF) * 1e-6)
        ddict['dF_err'].append(float(dF_err) * 1e-6)
        ddict['dF_F'].append(float(dF_F) * 1e-9)
        ddict['dF_F_err'].append(float(dF_F_err) * 1e-9)
        ddict['dF1'].append(float(dF1) * 1e-15)
        ddict['dF1_err'].append(float(dF1_err) * 1e-15)
        ddict['dF1_F1'].append(float(dF1_F1) * 1e-3)
        ddict['dF1_F1_err'].append(float(dF1_F1_err) * 1e-3)

    df = pd.DataFrame(ddict)
    df = df[df.dF != 0]

    df.to_csv(pickled_data)

    return df


if __name__ == "__main__":
    df = GetData()
    names = df.names.values
    dF = df.dF.values
    print "Using {} data points".format(len(names))

    log10dF = np.log10(dF)
    fig, ax = plt.subplots()
    out = ax.hist(log10dF, bins=13, histtype="step", color="k")
    ax.set_xlabel("$\log_{10}{\Delta f}$", size=22)
    ax.set_ylabel("Count", size=22)
    plt.savefig("EspinozaGlitchCat_dF.pdf")

    print "Maximum glitch size = {}".format(dF.max())
