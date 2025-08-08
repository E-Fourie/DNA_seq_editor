#Import modules
from snapgene_reader import snapgene_file_to_seqrecord
import os
from os import walk


#Take DNA sequence and print reverse complimentary of seq 
def reverse_comp_seq(seq):
    comp_dict = {'A': 'T', 'T': 'A', 'U': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(comp_dict.get(base, 'N') for base in seq[::-1].upper())

#Takes DNA size in kb and calcuates weight in ng for 50fmol, for Gibson calculations
def moles_to_weight(moles=50, input_DNA_size=1):
    input_DNA_size = float(input_DNA_size)
    moles = float(moles)
    final_mass = round(input_DNA_size * 649000 * moles / 10**6, 2)
    return str(final_mass) + " ng"

#Finds non cutting enzymes in DNA sequence
def enzyme_noncutter(file):
    seqrecord = snapgene_file_to_seqrecord(file).seq
    enzymeDict = {
    "GACGTC":"AatII",
    "AGCGCT":"AfeI",
    "CTTAAG":"AflII",
    "ACRYGT":"AgeI-HF",
    "GGGCCC":"ApaI",
    "GGCGCGCC":"AscI",
    "GCGATCGC":"AsiSI",
    "GCGATCGC":"AvaI",
    "CCTAGG":"AvrII",
    "GGATCC":"BamHI-HF",
    "GRGCYC":"BbsI",
    "CCTCAGC":"BbvCI",
    "AGATCT":"BglII",
    "CACGTC":"BmgBI",
    "GRCGYC":"BsaHI",
    "GGTCTC":"BsaI-HFv2",
    "WCCGGW":"BsaWI",
    "CGTACG":"BsiWI-HF",
    "CGTCTC":"BsmBI-v2",
    "GAATGC":"BsmI",
    "CTCAG":"BspDI",
    "TCATGA":"BspQI",
    "CCGCTC":"BsrBI",
    "CGCG":"BstUI",
    "RGATCY":"BstZ17I-HF",
    "GATC":"DpnI",
    "CTCTTC":"EarI",
    "GAATTC":"EcoRI-HF",
    "GAATTC":"EcoRV-HF",
    "AAGCTT":"HindIII-HF",
    "GTTAAC":"HpaI",
    "GGCGCC":"KpnI-HF",
    "GAAGA":"MfeI-HF",
    "AATT":"MluI",
    "AATT":"MluI-HF",
    "GAGTC":"MlyI",
    "CCATGG":"NcoI-HF",
    "GCCGGC":"NheI-HF",
    "GCGGCCGC":"NotI-HF",
    "TCGCGA":"NsiI-HF",
    "TTAATTAA":"PacI",
    "GTTTAAAC":"PmeI",
    "CACGTG":"PmlI",
    "VCTCGAGB":"PstI-HF",
    "CGATCG":"PvuI-HF",
    "CGGWCCG":"SacI-HF",
    "CCGCGG":"SalI-HF",
    "CCCGGG":"SmaI",
    "TACGTA":"SpeI-HF",
    "GCATGC":"SphI-HF",
    "TCTAGA":"XbaI"
    }
    enzymeKeyList = list(enzymeDict)
    noncutters = []
    for i in enzymeKeyList:
        if i not in seqrecord:
            noncutters.append(enzymeDict[i])
    return noncutters

def fasta_converter():
    abspath = os.path.abspath(__file__)
    current_dir = os.path.dirname(abspath) 
    os.chdir(current_dir)

    # Walk directory and exclude __file__ from list 
    inputFileNames = next(walk(current_dir), (None, None, []))[2]  # [] if no file
    inputFileNames = [file for file in inputFileNames if file != os.path.basename(__file__)]  # Exclude the script file
    seqSizesList = []

    # Convert DNA file to fasta format and record DNA size
    for x in range(len(inputFileNames)):
        currentPlamsid = inputFileNames[x]
        if currentPlamsid[-3:] == "dna":
            readSeqRecord = snapgene_file_to_seqrecord(os.path.join(current_dir, currentPlamsid)) 
            seqRecord = readSeqRecord.seq
            with open(currentPlamsid[:-4] + ".fa", "w") as f:
                f.write(">" + currentPlamsid[:-4] + "\n" + str(seqRecord))
            seqSizesList.append(inputFileNames[x] + ": " + str(len(seqRecord)))

    # Print DNA size to screen
    for line in seqSizesList:
        print(line)
    if len(seqSizesList) == 0:
        print("No files converted")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("DNA sequence editing tools")
    parser.add_argument("--reverse_complentarity", "-r", type=str, help="Returns the reverse complimentarity of a DNA sequence. Usage: -r SEQUENCE")  
    parser.add_argument("--moles_to_weight", "-m2w", nargs='+', type=float, help="Returns the weight of a given size and moles of DNA. Usage -m2w MOLES(fmol) DNA_SIZE(kb)")
    parser.add_argument("--enzyme_noncutter", "-enc", type=str, help="Returns the enzymes which don't cut a DNA sequence. Usage -enc FILE")
    parser.add_argument("--fasta_converter", "-fc", action='store_true', help="Returns size of the DNA files in the same folder and converts them to a fa file format. Usage -fc")
    args = parser.parse_args()
    if args.reverse_complentarity:
        print(reverse_comp_seq(args.reverse_complentarity))
    elif args.moles_to_weight:
        print(moles_to_weight(args.moles_to_weight[0], args.moles_to_weight[1]))
    elif args.enzyme_noncutter:
        print(enzyme_noncutter(args.enzyme_noncutter))
    elif args.fasta_converter:
        fasta_converter()
    
