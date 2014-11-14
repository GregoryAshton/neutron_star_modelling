# Python script to turn cleaned OCR read txt into data base

import pandas as pd
import glob


def txtTodf(file_name):
    f = open(file_name, "r")

    names = []
    peaks_s = []
    peaks_phase = []

    lists = [names, peaks_s, peaks_phase]
    i = 0
    j = 0
    for row in f:
        items = row.split(" ")
        if len(items) > 1 and items[0][0].isalnum():
            for item in items:
                lists[i].append(item.rstrip("\n"))
            i+=1
        else:
            i = 0

        if i > 3:
            print "Houston we have a problem.. with the file at line {}".format(j)
            break
        j +=1
    df = pd.DataFrame({'names' : names,
                       'peaks_s' : peaks_s,
                       'peaks_phase' : peaks_phase,
                       'source_file' : [file_name for i in range(len(names))]
                       })
    return df

if __name__ == "__main__":
    file_names = glob.glob("stripped-cleaned-*txt")
    df = None
    for file_name in file_names:
        if df:
            df = df.append(txtTodf(file_name), ignore_index=True)
        else:
            df = txtTodf(file_name)
    print df
    df.to_csv("HobbsLyneKramer2010Figure3Data.dat", sep=" ")

