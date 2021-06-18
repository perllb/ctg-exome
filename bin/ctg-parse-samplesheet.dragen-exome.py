#!/opt/conda/bin/python

# import libs
import csv

# import libs
import csv
import sys, getopt
import os

def read_arguments(argv):

    insheet = ''
    outsheet = ''
    indextype = ''

    usage="> Usage: ctg-parse-samplesheet.dragen-exome.py -s INPUT-SHEET -o OUTPUT-SHEET -i INDEX-TYPE [ -h HELP ] \n\n> Required columns (with the following header names): \n [Lane,Sample_ID,index,index2,Sample_Project]. \n - 'Lane' entries can be left blank if all lanes included. \n"

    try:
        opts, args = getopt.getopt(argv,"hs:o:i:",["insheet=", "outsheet=", "indextype="])
    except getopt.GetoptError:
        print(usage)
        sys.exit(2)
    if len(sys.argv) <= 3:
        print("> Error: please specify all arguments:")
        print(usage)
        sys.exit()
    for opt, arg in opts:
        if opt == '-h':
            print(usage)
            sys.exit()
        elif opt in ("-s", "--insheet"):
            insheet = arg
        elif opt in ("-o", "--outsheet"):
            outsheet = arg
        elif opt in ("-i", "--indextype"):
            indextype = arg   
            
    return insheet, outsheet, indextype

def main(argv):

    insheet, outsheet, index = read_arguments(argv)

    with open(outsheet, 'w') as outfile:
        writer = csv.writer(outfile)

        with open(insheet, 'r') as infile:
            my_reader = csv.reader(infile, delimiter=',')
            # row counter to define first line
            row_idx=0                # if first line - get index of the 3 columns needed
            datareached=0
            
            for row in my_reader:
                # read header
                if datareached == 0:
                    if 'Adapter' in row:
                        writer.writerow(['AdapterRead1',row[1]])
                    else:
                        writer.writerow(row)

                # if [Data] reached
                if datareached == 1:
                    datareached=2

                if datareached == 2:
                    if row_idx == 0:
                        sididx  = row.index('Sample_ID')
                        idxidx  = row.index('index')
                        idx2idx = row.index('index2')
                        projidx = row.index('Sample_Project')
                        row_idx = 1
                        if index == 'dual':
                            writer.writerow(['Sample_ID','Sample_Name','Sample_Plate','Sample_Well','index','index2','Sample_Project'])
                        else:
                            writer.writerow(['Sample_ID','index','Sample_Project'])
                    else:
                        currsid = row[sididx]
                        curridx = row[idxidx]
                        curridx2 = row[idx2idx]
                        currproj = row[projidx]

                        if index == 'dual':
                            writer.writerow([currsid,currsid,'','',curridx,curridx2,currproj])
                        else:
                            writer.writerow([currsid,curridx,currproj])
                if '[Data]' in row:
                    datareached=1             

if __name__ == "__main__":
    main(sys.argv[1:])
