import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as m


def GetData():

    names = ["NAME", "F0", "F0_err", "F1", "F1_err", "F2", "F2_err",
             "P0", "P0_err", "P1", "P1_err", "BINARY", "TYPE", "AGE"]
    df_full = pd.read_csv("ATNF_data_file.txt", sep=" ", skiprows=[0, 1, 2, 3],
                          skipinitialspace=True, skipfooter=1, engine='python',
                          names=names, index_col=False,
                          na_values="*")

    df_full = df_full.query("P0 > 20e-3 & P1 > 1e-17").copy().reindex()
    df_full['F2_rel_err'] = np.abs(df_full['F2_err'] / df_full['F2'])
    df_full = df_full.query("F2_rel_err < 0.75").copy().reindex()
    df_full = df_full[df_full.BINARY.isnull()].copy().reindex()
    df = df_full[~np.isnan(df_full['F2'])].copy().reindex()

    df['nobs'] = df.F0 * df.F2/df.F1**2
    df['tch'] = -df.F0 / (2*df.F1)
    df['tch_years'] = df['tch'] / (86400 * 365)

    print df.info()
    return df


def NobsPlot(df):
    fig, ax = plt.subplots()
    ax.semilogx(df['tch_years'], df['nobs'], "ok")
    ax.set_yscale("symlog")
    ax.set_xlabel(r"$\tau_{\mathrm{ch}}$ [years]")
    ax.set_ylabel(r"Observed braking index $n_{\mathrm{obs}}$")
    return ax

if __name__ == "__main__":
    df = GetData()
    ax = NobsPlot(df)
    plt.savefig("BirukovData.pdf")
