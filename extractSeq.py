from symbol import except_clause
from Bio import SeqIO
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-i", "--input", required=True,
                    help="input fasta sequence")
parser.add_argument("-p", "--positions", required=True,
                    help='positions to extract. Must be in double quotes. Whitespaces doesn\'t matter. Example: "chr12: 123 - 321"')


args = parser.parse_args()

fasta = args.input
positions = args.positions

import sys
from Bio import SeqIO
from argparse import ArgumentParser

# try:
#     withoutChr = positions.split(":")[1].replace(" ","")
# except:
#     withoutChr = positions.replace(" ","")
# startpos = int(withoutChr.split("-")[0].replace(',',""))
# endpos = int(withoutChr.split("-")[1].replace(',',""))
# #chr = positions.split(":")[0]

# Extracting positions, by avoiding eventual chr 
try:
    withoutChr = positions.split(":")[1].replace(" ","")
except:
    withoutChr = positions.replace(" ","")
startpos = int(withoutChr.split("-")[0].replace(',',""))
endpos = int(withoutChr.split("-")[1].replace(',',""))

# Extracting chr name. If it is not indicated by 'chr' add it.
# If no chr indicator, assume single sequence id in fasta. Checking later. 

if ":" in positions:
    chr = positions.split(":")[0]
else:
    chr = ""



if len(chr) != 0:    
    # Checking that chr starts with a letter or fixing it by add "chr"
    if not str(chr)[0].isalpha():
        chr = "chr"+str(chr)
    else:
        chr = chr

print("Looking for: \n\n" +chr+ "\nStartposition: " +str(startpos)+ " (inclusive)\nEndposition: " +str(endpos)+ " (inclusive)") 

if int(startpos) > int(endpos):
	print("\nEnd position can't be less than start position.")
	exit

if len(list(SeqIO.parse(fasta, "fasta"))) > 1:
    
    if len(str(chr)) < 1:
        sys.exit("\nFasta contains multiple sequence IDs/chromosomes. Define one infront of positions followed by a colon.")


    print("All checks passed. Extracting...")

    records = list(SeqIO.parse(fasta, "fasta"))

    for item in records:
        if item.id == chr:
            sequence = item.seq[startpos-1:endpos]
            fullSeq = item

    if 'sequence' not in globals():
        sys.exit("\nSequence ID/chromosome does not exist or positions are out of bounds")

    if endpos > len(fullSeq[:]):
        sys.exit("\nEndposition is bigger than length of sequence (" +str(len(fullSeq[:]))+").")


    location = fasta.split('.')[0]
    
    with open(location+'_'+chr+'_'+str(startpos)+'-'+str(endpos)+".fasta", "w") as text_file:
        text_file.write(">"+chr+'_'+str(startpos)+'-'+str(endpos)+"\n")
        text_file.write(str(sequence)+"\n")
    print("Done. Output in fasta folder.")

else:
    if ":" in positions:
        print("\nOnly one sequence in fasta. Ignoring given ID.")
    print("Searching fasta...")

    item = SeqIO.parse(fasta, "fasta")

    # getting record
    for i in item:
        record = i

    sequence = record.seq[startpos-1:endpos]

    if 'sequence' not in globals():
        sys.exit("\nSequence ID/chromosome does not exist or positions are out of bounds")

    if endpos > len(record.seq[:]):
        sys.exit("\nEndposition is bigger than length of sequence (" +str(len(record.seq[:]))+").")

    location = fasta.split('.')[0]
    fastaname = fasta.split('.')[0].split("/")[-1]
    
    with open(location+'_'+str(startpos)+'-'+str(endpos)+".fasta", "w") as text_file:
        text_file.write(">"+fastaname+"_"+str(startpos)+'-'+str(endpos)+"\n")
        text_file.write(str(sequence)+"\n")

    print("Done. Output in fasta folder.")
