import matplotlib.pyplot as plt
import numpy as np
import os.path as p
import json

binary = "./SiWiR2_EX1"
levels = [3, 4, 5, 6, 7, 8]
number_of_iterations = 15

json_path = "times.json"
time_json = {str(level): {} for level in levels}


def get_array_of_times():
    out = []
    for lvl in levels:
        out.append(time_json[str(lvl)]["time"])
    return out


def get_array_of_err():
    out = []
    for lvl in levels:
        out.append(time_json[str(lvl)]["error"])
    return out


def plot_all():
    times = get_array_of_times()
    plt.plot(levels, np.array(times))

    plt.title("Runtime for different levels")
    plt.xlabel("Level")
    plt.ylabel("Runtime in s")
    plt.show()


if __name__ == '__main__':
    if p.exists(json_path):
        with open(json_path, 'r') as file:
            time_json = json.load(file)
        plot_all()
    else:
        print("no json you idiot run timing.py first")
