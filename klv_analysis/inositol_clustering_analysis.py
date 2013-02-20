import csv

with open('clust_info.csv', 'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    for row in reader:
        indices = row[2].split(' ')
        print indices