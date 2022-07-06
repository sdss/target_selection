# temp1 has output of 
# \d catalogdb.sdss_apogeeallstarmerge_r13

output_file = "arraycols_sdss_apogeeallstarmerge_r13.out"
fin = open("temp1", "r")
fout = open(output_file, "w")

for line in fin:
    if "[" not in line:
        tags = line.split()
        if(len(tags) > 4):
            fout.write("alpha." + tags[0] + "," + "\n")
        
fout.close
fin.close()
