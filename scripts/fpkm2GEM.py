import os
import pandas as pd
import re


# Modify these parameters:
chunksize = 10

start_dir = '/scidas/arabidopsis/ath_PRJNA313379/sra2gev/'

#Test start_dir
# start_dir = './'

# Names file should be 2 column tab deliminated file with Run number in column
# one and run name in column 2, This is to name the columns.

names_file = '/scidas/arabidopsis/ath_PRJNA313379/sra2gev/PRJNA313379_Names.csv'

#Test names_file
# names_file = './Trial_Names'


# This is what you want your ematrix to be called at the end. You can select
# file path as well if you so desire.
ematrix_name = "GEM.txt"



os.chdir(start_dir)

#Finds all .fpkm files
fpkm_files = []
for root, dirs, files in os.walk(start_dir):
    #This causes our program to ignore the 'work' folder, so that we do not get
    #duplicates
    dirs[:] = [d for d in dirs if d not in 'work']
    for file in files:
        if file.endswith(".fpkm"):
             fpkm_files.append(os.path.join(root, file))


# extract fpkm column and add to new ematrix file
#ematrix file is established that only has "gene" as empty column.
ematrix = pd.DataFrame({'gene':[]})
for fpkm in fpkm_files:
    df = pd.read_csv(fpkm
                    #, chunksize=chunksize
                    , header=None
                    , sep='\t')

    #This gets the run name. It extracts it from the file name.
    file_basename = os.path.basename(fpkm)
    run_name = file_basename.split('_vs_')[0]

    # This splits up that grabage column that got made. I should talk to Dr. Ficklin about this one.
    df2 = df.join(df[8].str.split(' ', expand=True).rename(
                columns={0:10,1:11,2:12,3:"gene",4:14,5:15,6:16,7:run_name,8:18,9:19}))
    # print(df2)

    # Cleans up our columns of interest, we ignore all other columns becuse
    # They are not important:
    df2.iloc[:,12] = df2.iloc[:,12].str.replace(';', '')
    df2.iloc[:,12] = df2.iloc[:,12].str.replace('"', '')
    df2.iloc[:,16] = df2.iloc[:,16].str.replace(';', '')
    df2.iloc[:,16] = df2.iloc[:,16].str.replace('"', '')

    #Combine the file into the ematrix
    ematrix = pd.merge(ematrix, df2.iloc[:,[12,16]], on='gene', how='outer')
    print(file_basename)
#Rename rows of ematrix:
ematrix = ematrix.set_index('gene')
#Rename the ematrix:
names_file = pd.read_csv(names_file
               , header=None
               , sep='\t')
names_dict = names_file.set_index(0).to_dict()
print(names_dict)
#I did this command becasue the column remap was not working with dictionary for some reason
#pandas glitch or something
ematrix.columns = ematrix.columns.to_series().map(names_dict[1])

#This exports our ematrix
ematrix.to_csv(ematrix_name, sep = '\t', na_rep="NA")
