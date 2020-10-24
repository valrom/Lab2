import subprocess
import csv

subprocess.run(['mpic++', 'blocked.cpp'])

data_range = [64, 256, 1024, 4096, 8192]

mass = []
for j in data_range:
    print(subprocess.check_output(['mpiexec', '-np', '4', './a.exe', str(j)]).decode().split())

#myFile = open('results_blocked.csv', 'w')
#with myFile:
#    writer = csv.writer(myFile)
#    writer.writerows(mass)
