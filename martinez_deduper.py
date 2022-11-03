#!/usr/bin/env python

#sort file w/ samtools in bash before inputing to this script

import argparse
import re

#create input arguments for input samfile, input umi file, output file and he for help, using "h" or "help" causes error
def get_args():
    parser = argparse.ArgumentParser(description = "client input program")
    parser.add_argument("-i", "--file", help = "input filename", type=str, required = True)
    parser.add_argument("-o", "--outfile", help = "output filename", type=str, required = False)
    parser.add_argument("-u", "--umi", help = "umi filename", type=str, required = True)
    return parser.parse_args()

#make arg input name easier to type
clinput = get_args()
i = clinput.file
o = clinput.outfile
u = clinput.umi
he = clinput.helpcommand


#create function to adjust pos by the amount of left clipping or right clipping
def parse_cigar(pos: int, cigar: str, strand:str) -> int:

    #make empty variabls we will populate
    s_count = 0
    d_count = 0
    n_count = 0
    m_count = 0
    rclip_int = 0

    if strand == "+": #for forward strand grab regex for first soft clipping
        lclip_list = re.findall('[0-9]+(?=S)', cigar)
        if len(lclip_list) == 2: #if there are 2 "S" in cigar string then add the first to list
            s_count = int(lclip_list[0])
        elif len(lclip_list) ==1: #if there is only one "s" and it is not the last character then add to list
            if cigar[-1] != "S": #DRE IF CHECK IF WORKS added brackets
                s_count = int(lclip_list[0])

        pos = pos - s_count #subtract the number before our isolated correct "s" from our position


    #for reverse strands, regex to find s, n, m counts and add them to their own list, convert that list to an int and add them to our starting pos
    elif strand == "-":
        rclip_list = re.findall('[0-9]+(?=S$)', cigar)
        if cigar[-1] == "S":
            rclip_str = rclip_list[0]
            rclip_int = int(rclip_str)

        d_list = re.findall('[0-9]+(?=d)', cigar)
        for value in d_list:
            d_count += int(value)

        m_list = re.findall('[0-9]+(?=m)', cigar)
        for value in m_list:
            m_count += int(value)

        n_list = re.findall('[0-9]+(?=n)', cigar)
        for value in n_list:
            n_count += int(value)
    
    pos = pos + d_count + m_count + n_count + rclip_int

    return pos

#create empty sets, one for umis, one that will be our output
umi_set = set()
out_set = set()

#add every entry in our input umi file to our umi set as long as it isn't in there already
with open (u, "r") as umi_file:
    for line in umi_file:
        if not line in umi_set:
            umi_set.add(line.strip())

# loop through input sam file and save all of the sam fields that we care about in to variables. save chromname, strandedness and position to set
iterate = 0
countumi = 0
with open (o, "w") as outfile:
    with open (i, "r") as samfile:
        for line in samfile:

            if line.startswith('@'): #if its a header line write it to the new file
                outfile.write(line)

            elif not line.startswith('@'): #if its not a header ine examine for the fileds we want
                line = (line.strip('\n')) #remove new lines
                cols = line.split('\t') #split in to new lines at each tab
                umis = cols[0] [-8:] #grab the UMI from the qname
                flag = int(cols[1])
                if ((flag & 16) == 16): #checks flag and assigns plus or minus conditionally
                    strand = "-"
                else: 
                    strand = '+'
                chrom = cols[2] 
                pos = int(cols[3])
                cigar = cols[5] #tells us if skipped or not
                iterate += 1
                adj_pos = parse_cigar(pos, cigar, strand) # call our function parse_cigar on the three inputs it expects and assign the functions return to adj_pos

                if umis not in umi_set: #if the umi isn't in the umi_set skip (continue) to the next iteration
                    countumi += 1
                    continue

                if (umis, chrom, adj_pos, strand) not in out_set: #if we haven't already written this umis+chrom+adj_pos+strand to our output set then add it to the set
                    out_set.add((umis, chrom, adj_pos, strand))
                
                    outfile.write(line+"\n") #write the set but add in new lines

print(countumi)
                
