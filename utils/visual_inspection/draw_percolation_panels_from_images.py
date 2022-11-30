import numpy as np
import cv2 as cv
import argparse
import matplotlib.pyplot as plt

#temperature = [0.16,0.14,0.12,0.10,0.09,0.07,0.05,0.04,0.01]
#phi = [0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5]

if __name__ == '__main__':

    # read data either through files system via glob or via db
    parser = argparse.ArgumentParser()
    parser.add_argument('-ptype', type=str, choices=['double_manta_asymm_1',
                        'double_mouse_symm_1', 'double_mouse_symm_2', 'double_mouse_asymm_1'])
    parser.add_argument('-delta', type=float, choices=[0.2, 0.3, 0.4])

    args = parser.parse_args()

    temperature = [0.16, 0.14, 0.12, 0.10, 0.09, 0.07, 0.05, 0.04, 0.01]
    phi = [0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]

    N_temperature = len(temperature)
    N_phi = len(phi)

    fig, ax = plt.subplots(N_temperature, N_phi)

    print("start drawing")
    for row, pi in zip(ax, phi):
        for col, ti in zip(row, temperature):
            fname = "{}_phi_{}_delta_{}_temp_{}_rhombi.png".format(
                args.ptype, pi, args.delta, ti)
            print(fname)
            print(col)
            col.spines['top'].set_visible(False)
            col.spines['right'].set_visible(False)
            col.spines['left'].set_visible(False)
            col.spines['bottom'].set_visible(False)

            col.get_xaxis().set_visible(False)
            col.get_yaxis().set_visible(False)

            img = cv.imread(fname)
            col.set_title("$\phi=${},$T=${}".format(pi, ti))
            col.imshow(img)

    plt.savefig("{}_{}_snapshot.png".format(args.ptype, args.delta))
