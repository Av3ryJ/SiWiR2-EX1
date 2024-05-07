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
        out.append(time_json[str(lvl)]["err"])
    return out


def plot_all():
    times = get_array_of_err()
    plt.loglog(pow(2, np.array(levels)), np.array(times))

    plt.title("Loglog graph of error vs mesh-width")
    plt.xlabel("mesh width")
    plt.ylabel("Final Error")
    plt.show()


if __name__ == '__main__':
    if p.exists(json_path):
        with open(json_path, 'r') as file:
            time_json = json.load(file)
        plot_all()
        for i in range(3,9):
            print(pow(2,i*-1))
    else:
        print("no json you idiot run timing.py first")
