#!/home/marwene/.virtualenvs/BIO1/bin/python3

import requests
from bs4 import BeautifulSoup
import pandas as pd
from collections import Counter
import re
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-v", "--verbose", action="store_true",
                    help="By using the verbose option all the steps of the program "
                         "(downloading the genetic codes, "
                         "reading the file, translating the dna sequence to amino acids"
                         " and finding the sequence"
                         " with te lowest degeneration value will be shown and written"
                         " in the output file)")
parser.add_argument("-q", "--quite", action="store_true",
                    help="By choosing the quit option only the amino acid sequence with"
                         " the lowest degeneration value "
                         "will be printed")
parser.add_argument("-c", "--check", action="store_true",
                    help="This option will check if the entered dna sequence contains any "
                         "start/stop codons")
parser.add_argument("-o", "--output", help="By choosing this option the output will be printed"
                                           " in a file",
                    action="store_true")
parser.add_argument("input_filename",
                    help="Please write the full directory path leading to the file containing"
                         " the sequence to analyse")
args = parser.parse_args()


def create_file():
    """creating a file that will contain all the analysis of dna when choosing the flag -o"""
    with open("dna_analysis", "a") as out_file:
        out_file.write("*" * 50 + " " + "DNA analysis " + "*" * 50 + "\n")
        out_file.write("The used table in this analysis is downloaded from:\n{}\nIs the "
                       "following:\n{}\n".format(url, gc_table))
        out_file.write("The start/stop codons are {}\n".format(init_code))
        out_file.write("This dictionary shows the determined number of codons every amino"
                       " acid is encoded with\n{}\n".format(ct))
        out_file.write("After translation of the DNA, the obtained protein is {}\n".format(protein1))
        out_file.write("The sequence with the lowest degeneration number is {} with {} codons\n".format(
            protein[min_pos: min_pos + 15], min_value))
        out_file.write("The dna sequence with the lowest degeneration number is {}".format(dna[min_pos: min_pos + 45]))


def scrap_gc():
    global gc_dict
    global gc_table
    global url
    """This function takes the url containing the genetic code and organise it in a dictionary,
    where the keys are the codons and the values are the names of the amino acids,
    also it stores the initial/stop codons in the list init_code"""

    # url of the wikipedia page that contains the genetic code for microbes
    url = "https://en.wikipedia.org/wiki/Bacterial,_archaeal_and_plant_plastid_code"
    response = requests.get(url)
    soup = BeautifulSoup(response.text, "html.parser")
    genetic_code_html = soup.find_all("table")[1]
    genetic_code_list = pd.read_html(str(genetic_code_html), index_col=0)

    gc_table = pd.concat(genetic_code_list)

    # due to the form of the table some processing was needed to obtain the final result
    dict_u = {codon: AA for codon, AA in zip(gc_table[('2nd base', 'U')], gc_table[('2nd base', 'U.1')])}
    dict_c = {codon: AA for codon, AA in zip(gc_table[('2nd base', 'C')], gc_table[('2nd base', 'C.1')])}
    dict_a = {codon: AA for codon, AA in zip(gc_table[('2nd base', 'A')], gc_table[('2nd base', 'A.1')])}
    dict_g = {codon: AA for codon, AA in zip(gc_table[('2nd base', 'G')], gc_table[('2nd base', 'G.1')])}
    gc_dict = {**dict_u, **dict_c, **dict_a, **dict_g}
    if args.verbose:
        print("{} collecting the codons of amino acids from the genetic "
              "code database\n .....".format(list(gc_dict.values())[0:5]))
    return gc_dict, find_init_codons(gc_dict)


def find_init_codons(gc_dict):
    global init_code
    # storing the initial codons in a separate list
    init_code = []
    for key in list(gc_dict):
        if "[" in key:
            gc_dict[key[0:3]] = gc_dict.pop(key)
            init_code.append(key[0:3])
    if args.verbose:
        print("Collecting initiation codons...")
        print(init_code)
    return init_code


def read_file(file):
    file = args.input_filename
    global dna
    """this function read a .txt file containing dna and return it as a string containing the whole sequence of DNA"""
    with open(file, "r") as dna:
        dna = dna.read().strip("\n")
        return dna


def codon_number(gc_dict):
    global ct
    """Each amino acid is encoded by a determined number of codons, this function counts the number of codons
    for every amino acids and return the obtained values as dictionary"""
    ct = Counter()
    for key, value in gc_dict.items():
        ct[value] += 1
    return ct


def stop(dna):
    """This function checks if any of the stop codons exists in the dna fragment"""
    st_st = re.findall(r'init_code', dna)
    if len(st_st) > 0:
        print("This sequence contains {} stop/start codons")
    elif len(st_st) == 0:
        print("No start/stop codons were found in the entered dna sequence")


def translation(dna):
    """This function translates the dna into a sequence of amino acids,
    if the codon does not exist it will be replaced by a star"""
    global protein, protein1
    protein1 = []
    protein = []
    for base in range(0, len(dna) - 1):
        codon = dna[base: base + 3]
        if len(codon) == 3:
            protein1.append(gc_dict.get(codon, "*")[5:6])
            protein.append(gc_dict.get(codon, "*"))
    if args.verbose:
        print("The dna sequence has been translated to amino acids")
    return protein


def low_deg_zone(protein):
    global min_pos
    global min_value
    """This function finds the segment (made with 15 amino acids) with the lowest degeneration value,
    these regions are useful to find PCR primers the less codons an amino acid has the lower the degeneration value"""
    seg_values = []
    if args.verbose:
        print("Calculation of degeneration zones values has began...")
    for num in range(0, len(protein)):
        segment = protein[num: num + 15]
        degen = 0
        if len(segment) == 15:
            for aa in segment:
                degen += ct.get(aa, 3.05)
            seg_values.append(degen)
    min_value = min(seg_values)
    min_pos = seg_values.index(min_value)
    if args.verbose:
        print("The amino acid sequence with the lowest degeneration number is {} "
              "with {} codons\n".format(protein[min_pos: min_pos + 15], min_value))
        print("The dna sequence with the lowest degeneration number is {}".format(dna[min_pos: min_pos + 45]))
    elif args.quit:
        print(protein[min_pos: min_pos + 15] + "\n")
        print(dna[min_pos: min_pos + 45])


if __name__ == "__main__":
    scrap_gc()
    read_file("dna.txt")
    codon_number(gc_dict=gc_dict)
    if args.check:
        stop(dna)
    translation(dna)
    low_deg_zone(protein)
    if args.output:
        create_file()




