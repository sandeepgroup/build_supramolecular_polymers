import re
with open("crest_conformer_10.xyz", 'r') as fp:
   natom = fp.readline().strip()
   line_counter = 0
   conformer_count = 0
   conformer_coordinates = []
   for line in fp:
     if(len(re.findall(r"[-+]?(?:\d*\.*\d+)", line))>1):
       line = line.strip()
       line_counter += 1
       conformer_coordinates.append(line)
       if(line_counter == int(natom)):
         conformer_count += 1
         line_counter =0
         file_name = f'conformer_{conformer_count}.xyz'
         with open(file_name, 'w') as file:
           file.write(natom)
           file.write("\n\n")
           for coord in conformer_coordinates:
             formatted_values = ["{:>1f}".format(value) if isinstance(value, (int, float)) else "{:>1}".format(value) for value in coord]
             line = "".join(formatted_values) + "\n"
             file.writelines(line)
         file.close()
         conformer_coordinates.clear()
