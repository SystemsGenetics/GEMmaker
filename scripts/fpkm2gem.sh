# FPKM to GEM
#
# This is for making GEMs from the output of GEMmaker
# 
# This is the directory where GEMmaker output directories for
# Each experiment
startdir='../'

# Name of the output file
out_file='GEM.txt'

# This allows you to rename columns. If you want to rename, you
# need a 2 column tab delimited file with experiment name in 
# column 1 and your new name in column 2.
#
# Set rename to "No" to keep the name of the experiment.
rename="No"



# This line is to make sure that we do not add onto an already
# exisiting GEM_File. It will delete any file with the 
# 'out_file' name

cd $startdir

rm $out_file


# This combines all files in the subdirectories that are fpkm
for fpkm in `ls */*.fpkm`
do 
    echo $fpkm
    # Checks to see if user wants to rename columns
    if [ $rename == "No" ]
    then
        fpkm_name=`echo $fpkm |sed -r 's/(.*)(\/.*)/\1/g'`
    else
        fpkm_name=`echo $fpkm |sed -r 's/(.*)(\/.*)/\1/g' | grep -f - $rename|sed -r 's/(.*	)(.*)/\2/g'`
    fi
    # Takes FPKM and adds to ematrix
    if [ -e $out_file ]
    then
	    sort $fpkm > tmp_fpkm
            echo gene	$fpkm_name | cat - tmp_fpkm > tmp_fpkm_header
            join $out_file tmp_fpkm_header  > tmp_out
	    mv tmp_out $out_file
    else
            sort $fpkm > tmp_fpkm
	    echo gene   $fpkm_name | cat - tmp_fpkm > $out_file
    fi

done

# Makes our file tab demiated and removes 'gene' from line 1
sed 's/ /\t/g' $out_file > tmp_out
sed 's/gene	//' tmp_out >  $out_file

# Clean_up
rm tmp_fpkm

