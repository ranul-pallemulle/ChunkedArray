#!/usr/bin/python

import sys
import matplotlib.pyplot as plt


def load_data(filename):
    vTime = list()
    vV = list()
    with open(filename, 'r') as datafile:
        for line in datafile:
            time, V, _ = line.split('   ')
            time = time.split(':')[1].strip()
            V = V.split(':')[1].strip()
            vTime.append(float(time))
            vV.append(float(V))
    return vTime, vV


if __name__ == '__main__':
    _, ax = plt.subplots()
    ax.set(xlabel='time', ylabel='V', title='Action Potential')
    ax.grid()
    filenames = sys.argv
    for name in filenames[1:]:
        time, V = load_data(name)
        ax.plot(time, V, label=name)
    plt.legend()
    plt.show()
