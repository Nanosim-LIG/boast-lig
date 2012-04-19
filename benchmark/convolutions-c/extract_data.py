#!/usr/bin/env python
import re, sys

conv_num = re.compile(r'.([0-9]+).+')
papi_data = re.compile(r'.*(PAPI\S+)\s*:\s+([0-9]+)')
time_data = re.compile(r'Finished.+')

names = ["dcm.dat", "dca.dat", "tc.dat", "icm.dat", "tlb_dm.dat", "tlb_im.dat"]

for j in range (6):

  fout = open("extracted_data/" + names[j], "w")
  fout.write("conv\tvalue\n")

  for i in range(1,int(sys.argv[1]) + 1 ):

    fin = open("data/" + str(i) + "/" + names[j], "r")

    for line in fin:

        m = conv_num.match(line)
        if m != None:
          current = m.group(1)	
          fout.write(current + '\t')	
          continue

        m = papi_data.match(line)
        if m != None: 
          counter_type = m.group(1)
          counter_value = m.group(2)
          fout.write(counter_value + '\n')
          continue

        m = time_data.match(line)
        if m != None: 
          continue 

fin.close()
fout.close()
