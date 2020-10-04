# Search_for_low_degeneration_zone_in_dna
This application searches for the lowest degeneration zone in a sequence of dna, the behavior of the application changes with the used flags in the command line.
To find the best regions in a dna sequence to use it as PCR primers, this application gives a better chance to find this region by detecting the sequences with 45 nucleotides from the dna and 15 amino acids from the protein characterized with the lowest number of codons.
Each amino acids is encoded by a determined number of codons, hence a region wich contains amino acids with less number of representative codons is more conservative and suitable for PCR primers.
This application scrap the genetic code from the wikipedia page (for Bacteria and Archaea)  store it as a dictionary, reads a .txt file containing the dna sequence,
searches for start/stop codons, translate the sequence to a proteine and output the sequences with the lowest values of codon variability.
A command line graphical interface has been created with different flags. (-o for creating an output file, -v for showing the results on the command line while the program is running, -q flag to obtain only the results and -c to check for the presence of start/stop codons)
