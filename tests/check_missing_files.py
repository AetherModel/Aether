import os.path
from pathlib import Path

file_call_double = []
file_call_single = []

root_folder = Path(__file__).parents[1]
input_files = root_folder / "src/inputs.cpp"
default_file = root_folder / "share/run/UA/inputs/defaults.json"
aether_file = root_folder / "share/run/aether.json"

#read inputs.cpp and record all the settings it calls
#reads single item settings (eg ["Seed"]) in a separate list
#because it's easier to deal with that way
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

#function to check each json file for the settings called in inputs.cpp
def check_json_file(file, input_arr):
    #first checks to see if the single setting names are in the json file
    for item in file_call_single:
        if not item in file:
           input_arr += [item]
    
    #then checks for 2-part setting names
    setting_size = len(file_call_double) - 1
    for i in [x for x in range(setting_size) if x % 2 == 0]:
        if not file_call_double[i] in file:
           input_arr += ["[" + file_call_double[i] + ", " + file_call_double[i + 1] + "]"]
        else:
           start = file.find(file_call_double[i])
           end = file.find("}", start)
           if not file_call_double[i + 1] in file[start: end]:
               input_arr += ["[" + file_call_double[i] + ", " + file_call_double[i + 1] + "]"]

#open and check against aether.json
missing_settings_aether = []
with open(aether_file, "r") as a:
    file_string = a.read()
    check_json_file(file_string, missing_settings_aether)

#open and check against defaults.json
missing_settings_defaults = []
with open(default_file, "r") as j:
    file_string = j.read()
    check_json_file(file_string, missing_settings_defaults)

#cross checks missing aether.json settings against missing defaults.json settings
missing_settings = []
for item in missing_settings_defaults:
    if item in missing_settings_aether:
        missing_settings += [item]

print("Missing settings: ")
print(missing_settings)