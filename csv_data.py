from __future__ import print_function
import csv
import GTC as gtc

path = r'ubdict.csv' # Change this to suit
file = open(path,mode='rb') # Open in 'binary' mode
reader = csv.reader( file )
reader.next()
header = reader.next()
print(header)
table = gtc.number_strings.sequence_parser(reader)
file.close()
Nr = len(table)
Nc = len(header)
caldata = {}
acdc = {}
tempcoeffs = {}
stability = {}
ivd = {}
for r in range(Nr):
    dictionary = table[r][0]
    print(dictionary)
    assert dictionary in ['caldata', 'acdc','tempcoeffs', 'stability', 'ivd'], 'csv dictionary name not in list'
    key = table[r][1]
    value = gtc.ureal(table[r][2],table[r][3],table[r][4],table[r][5])
    if dictionary == 'caldata':
        caldata[key] = value
    elif dictionary == 'acdc':
        acdc[key] = value
    elif dictionary == 'tempcoeffs':
        tempcoeffs[key] = value
    elif dictionary == 'stability':
        stability[key] = value
    elif dictionary == 'ivd':
        ivd[key] = value

print('caldata', caldata)
print('acdc',acdc)
print('tempcoeffs',tempcoeffs)
print('stability',stability)
print('ivd', ivd)


