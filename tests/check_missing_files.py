import os
import os.path
from pathlib import Path

#from ..run.test/UA/inputs import defaults.json

file_call_double = []
file_call_single = []
missing_files = []

root_folder = Path(__file__).parents[1]
input_files = root_folder / "src/inputs.cpp"
json_files = root_folder / "share/run/UA/inputs/defaults.json"
aether_files = root_folder / "share/run/aether.json"

with open(input_files, "r") as f:
    for line in f:
        start = line.find('settings["')
        if(start != -1):
            middle = line.find('"]["', start)
            end = line.find('"]', middle + 1)
            if(middle == -1):
                file_call_single += [line[(start + 10): (end)]]
            else:
                file_call_double += [line[(start + 10): (middle)]]
                file_call_double += [line[(middle + 4): end]]



with open(json_files, "r") as j:
    json_string = j.read()
    for item in file_call_single:
        if not item in json_string:
           missing_files += [item]
    
    setting_size = len(file_call_double) - 1
    for i in [x for x in range(setting_size) if x % 2 == 0]:
        if not file_call_double[i] in json_string:
           missing_files += ["[" + file_call_double[i] + ", " + file_call_double[i + 1] + "]"]
        else:
           start = json_string.find(file_call_double[i])
           end = json_string.find("}", start)
           if not file_call_double[i + 1] in json_string[start: end]:
               missing_files += ["[" + file_call_double[i] + ", " + file_call_double[i + 1] + "]"]

print("Missing files are listed below: ")
print(missing_files)