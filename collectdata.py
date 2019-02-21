import ugradio
import numpy as np
import matplotlib.pyplot as plt


def collect(volt_range, blocks, numrepeat, name):
    """volt_range: str """

    for i in range (0,numrepeat):
        data = ugradio.pico.capture_data(volt_range, divisor=10, dual_mode=True, nblocks=blocks)
        np.savetxt(name+str(i),data)

def insert_upper(volt_range):
    data = ugradio.pico.capture_data(volt_range, divisor=10, dual_mode=True)
    np.savetxt('uppertestsig', data)

    return np.loadtxt('uppertestsig')


def insert_lower(volt_range):
    data = ugradio.pico.capture_data(volt_range, divisor=10, dual_mode=True)
    np.savetxt('lowertestsig', data)

    return np.loadtxt('lowertestsig')



