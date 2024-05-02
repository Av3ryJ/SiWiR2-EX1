import subprocess
import os.path as p
import json
import regex

binary = "./mgsolve"
levels = [3, 4, 5, 6, 7, 8]
number_of_iterations = 10

json_path = "times.json"
time_json = {level: 0.0 for level in levels}


def run_bin(level):
    result = subprocess.run([binary, str(level), str(number_of_iterations)], capture_output=True, text=True)
    return result.stdout

def parse_output(output):
    return output

def time_all():
    # iterate over sizes and threads
    for thread_num in thread_numbers:
        print(f"running: {thread_num}")
        returned = run_bin(str(thread_num))
        parsed = parse_output(returned)
        print(f"time was: {time}")
        time_json[thread_num] = time
    with open(json_path, 'w') as f:
        json.dump(time_json, f)


if __name__ == '__main__':
    if p.exists(json_path):
        with open(json_path, 'r') as file:
            time_json = json.load(file)
    else:
        time_all()