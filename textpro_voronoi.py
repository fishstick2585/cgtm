import numpy as np
import re
# Read list of resids; split into array of lists and convert to integers
input_file2 = open('processHop.out','r')
full_file_contents = input_file2.read()
rts_intermed = re.split('First shell: |\nSecond shell: |\n', full_file_contents)
#for i in lines_of_fie:
#    lines_of_file[i] = lines_of_file[i].split('Second shell: ')
#for i in lines_of_file, # here yo need to split these too by 'Second shell: '
#rts_intermed= np.array(lines_of_file)
#num_shellana = 2;
range_cuta = (len(rts_intermed) - 1)
range_cutb = range_cuta*2
range_cut = range_cutb/4
print(rts_intermed)
range_cut = range_cuta*2/4
print(range_cut)
rts_fs = [0 for i in range(range_cut/2)]
rts_ss = [0 for i in range(range_cut/2)]

for i in range(0, range_cuta, 4):
    j = i*1/4
    rts_fs[j] = rts_intermed[i + 2]
    rts_ss[j] = rts_intermed[i + 3]
#size_rts = len(rts);
#print(size_rts)
print(rts_fs)
print(rts_ss)
file = open("rts_fs.txt", "w")
for index in range(len(rts_fs)):
    file.write(str(rts_fs[index]) + " " + str(rts_ss[index]) + "\n") # + " " + str(b[index]) + "\n")
file.close()
