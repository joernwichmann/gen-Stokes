import csv

def read_datafile(location: str):
    with open(location, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        #skip header
        next(reader)

        output = []
        for row in reader:
            row_proc = [float(entry) for entry in row]
            output.append(row_proc) 
    return output
   
def organize_output(output):
    time = []
    stepsize = []
    mean = []
    eoc_mean = []
    sd = []
    for row in output:
        time.append(row[0])
        stepsize.append(row[1])
        mean.append(row[2])
        eoc_mean.append(row[3])
        sd.append(row[4])
    
    return time, stepsize, mean, eoc_mean, sd