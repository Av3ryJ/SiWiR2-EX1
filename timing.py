import subprocess
import os.path as p
import json
import re

binary = "./SiWiR2_EX1"
levels = [3, 4, 5, 6, 7, 8]
number_of_iterations = 15

json_path = "times.json"
time_json = {level: {} for level in levels}


def run_bin(level):
    result = subprocess.run([binary, str(level), str(number_of_iterations)], capture_output=True, text=True)
    return result.stdout


def parse_output(output):
    output_format = {"iterations": 0, "time": 0.0, "err": 0.0, "norms": []}
    regex_string_lvl_iter = r"Level: (\d+).* Iterations: (\d*)"
    regex_string_time_err = r"runtime: ([\d|\.|e|-]*); total error: ([\d|\.|e|-]*)"

    temp = re.search(regex_string_lvl_iter, output)
    # output_format["level"] = int(temp.group(1))
    output_format["iterations"] = int(temp.group(2))
    temp = re.search(regex_string_time_err, output)
    output_format["time"] = float(temp.group(1))
    output_format["err"] = float(temp.group(2))

    for i in range(1, output_format["iterations"]+1):
        search_str = rf"iteration {i}: discrete residual norm = (.*);"
        temp = re.search(search_str, output)
        output_format["norms"].append(float(temp.group(1)))

    print(output_format)
    return output_format


def time_all():
    # iterate over sizes and threads
    for level in levels:
        print(f"running: {level}")
        returned = run_bin(str(level))
        parsed = parse_output(returned)
        time_json[level] = parsed
    with open(json_path, 'w') as f:
        json.dump(time_json, f)


if __name__ == '__main__':

    if p.exists(json_path):
        print("times.json already exists!")
        with open(json_path, 'r') as file:
            time_json = json.load(file)
    else:
        print("Timing...")
        time_all()
        print("---Done---")
