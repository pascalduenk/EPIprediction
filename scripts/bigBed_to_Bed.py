#!usr/bin/env  python3 
'''
a machine learning script

Author: yuhao wang

This script was used to transform bigbed files to bed files using bigBedToBed tool provided by UCSC

'''

# import statements
from sys import argv
import os
import subprocess

#functions
def file_name(file_dir):
    """ return name lists of ACAT,WGBS,CHIP files

    Attention: 
        for ATAC_files and others, ["ATAC", names...]
    """
    for root,dirs,files in os.walk(file_dir):
        #print(root)
        #print('---||')
        #print(dirs)
        #print("---||---")
        #print(files)
        #print("\n\n\n")

        if root == "./ATAC":
            files.insert(0,"ATAC/")
            ATAC_files = files
        if root == "./WGBS":
            files.insert(0,"WGBS/")
            WGBS_files = files
        if root == "./CHIP-seq":
            files.insert(0,"CHIP-seq/")
            CHIP_files = files
    return ATAC_files,WGBS_files,CHIP_files

def parse_file_names(list_files):
    """ extract the prefix names from the list of names

    """
    new_list = []
    for name in list_files:
        new_list.append(name[:-7])

    return new_list



def use_bigBedToBed(root_name,file_name):
    """  using bigBedToBed tool to transform the bigBed  to bed format
    
    """
    out_file = root_name + file_name[:-7] + ".bed"
    if os.path.exists(out_file):
        print("This file have been transformed! ")
        return out_file
    else:
        input_file = root_name + file_name
        cmd = "./bigBedToBed {} {}".format(input_file,out_file)
        e = subprocess.check_output(cmd, shell= True)
        print("Done")
        return out_file


#main
if __name__ == '__main__':

    #step1: return lists containing root name and file names in root
    atac_files,wgbs_files,chip_files = file_name("./")
    print("\n\n")

    #step2: transform bigBed of ATAC files to Bed files 
    for name in atac_files[1:]:
        use_bigBedToBed(atac_files[0],name)

    #step3: transform bigBed of CHIP-seq files to Bed files.
    for name in chip_files[1:]:
        use_bigBedToBed(chip_files[0],name)





   