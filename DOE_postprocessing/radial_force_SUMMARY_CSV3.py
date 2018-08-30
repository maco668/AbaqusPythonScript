############################################################################################
#                                                                                          #
#  Script to copy and paste the fatigue results from the csv file to a summary file        #
#                                                                                          #
############################################################################################
import csv
import re
#------------------------------- MAX IFU --------------------------------------------------
Outputfile = open('RadialForce_SUMMARY_MaxIFU.csv', 'wb') # If you don't open the file in binary mode ('wb'), python will create a blank line between the data 
ofile=csv.writer(Outputfile)
realRow=re.compile("[0-9]+")
data=[]
for index in range(1,3):
    x_file = open('C:\\Users\\llargura\\Documents\\GORE_Projects\\ULP_Limb\\27mm\\Radial_Force_Results\\RUN_1_'+str(index)+'\\Radial_Force_IFUs_ForceData_MAXIFU.csv')#Update the file path
    results = csv.reader(x_file)
    i=0
    column=[]
    for row in results:
        if realRow.match(row[0])!=None:
            results2 = row[2]
            column.append(results2)
            print(results2)
    x_file.close()
    data.append(column)
n=map(list,zip(*data))
ofile.writerows(n)
Outputfile.close()
print('Done!')
