#!/home/data1/tools/bin/anaconda3/bin/python3
# ./toxfinder-1.0.py -i ../../../ToxFinderProto/test_data/Aspergillus_flavus_NRRL3357_small.fsa -o ../../../ToxFinderProto/test_data/ -p ../../database_archive/toxfinder_db/ -b /home/data1/tools/bin/blastn

import sys, os, re, math
from argparse import ArgumentParser

from cgecore.blaster import Blaster
import time

################################################################
#			FUNCTIONS
################################################################


################################################################
#			ARGUMENT PARSING
################################################################

parser = ArgumentParser(description="This program...")

# positional arguments
parser.add_argument("-i", "--inputfiles", nargs ='+', help="Input fasta file(s)") # Don't konw if default should be sat to anything
parser.add_argument("-o", "--out_path", help="Path to blast output") # er det runroot eller er det tmp directiory

# optional arguments
parser.add_argument("-p", "--databasePath", dest="db_path",help="Path to the databases", default='') #/home/data1/services/ResFinder/database_archive/database-3.0_Test_indels

#parser.add_argument("-m", "--method", dest="method",help="blastn or kma", choice={'kma','blastn'}, default='blastn')
#parser.add_argument("-m_p", "--method_path", dest="method_path",help="path to executable method, blastn or kma", default='blastn')

parser.add_argument("-b", "--blastPath", dest="blast_path",help="Path to blast", default='blastn')
parser.add_argument("-t", "--threshold", dest="threshold",help="Blast threshold for identity", default=98)
parser.add_argument("-l", "--min_cov", dest="min_cov",help="Minimum coverage", default=60)
parser.add_argument("-g", "--specific_genes", nargs ='+', dest="specific_genes",help="Specifie genes existing in the database to search for only - if non is specified all genes are used", default=None)

args = parser.parse_args()

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit()

# Check if valid input file is provided
for infile in args.inputfiles:
    if not os.path.exists(infile):
        sys.exit("Input Error: Input file, {:s} does not exist!\n".format(infile))
inputfile = args.inputfiles[0]

filename = "_".join(inputfile.split("/")[-1].split(".")[:-1])
#print(filename)

# Check if valid output directory is provided
if not os.path.exists(args.out_path):
    sys.exit("Input Error: Output dirctory does not exists!\n")
else:
    out_path = args.out_path

# Check if valid database path is provided
if not os.path.exists(args.db_path):
    sys.exit("Input Error: The specified database directory, %s, does not exist!\n"%args.db_path)
else:
    db_path = args.db_path
    gene_list = [db_file[:-4] for db_file in os.listdir(db_path) if db_file.endswith(".fsa")]
    if gene_list == []:
        sys.exit("Input Error: The specified database directory, %s, does not contain any '.fsa' files!\n"%args.db_path)

# Open gene list
# Getting phenotype hash
phenos = dict()
with open(db_path+"/notes.txt", 'r') as f:
   for line in f:
      if line.startswith("#"):
         continue
      else:
         tmp = line.split(":")
         phenos[tmp[0]] = [tmp[1], tmp[2]]


# Check if valid blast path is provided
if not os.path.exists(args.blast_path):
    sys.exit("Input Error: The specified blast path, %s, does not exist!\n"%args.blast_path)
else:
    blast = args.blast_path

gene_list = [db_file[:-4] for db_file in os.listdir(db_path) if db_file.endswith(".fsa")]
 

# Create user defined gene_list if applied 
if args.specific_genes:
    # Check that the genes are valid
    genes_specified = []
    for gene in args.specific_genes:
        if gene in gene_list:
            genes_specified.append(gene)
        else:
            print("Input Error: Provided database was not recognised! (%s)\
                   \nYou can between the following genes:\n"%(gene))
            for gene in gene_list:
                print(gene)
            sys.exit()
    # Change the gene_list to the user defined gene_list
    gene_list = genes_specified
 	
min_cov =  float(args.min_cov)/100.0
threshold = float(args.threshold)/100.0

################################################################
#			MAIN PROGRAM
################################################################

#print(inputfile, gene_list, db_path, out_path, min_cov, threshold, blast)

blast_run = Blaster(inputfile, gene_list, db_path, out_path, min_cov, threshold, blast, cut_off=False)
sbjct_align = blast_run.gene_align_sbjct
query_align = blast_run.gene_align_query
results = blast_run.results
count= 0

# Making output files
tab_file = []#open(out_path+"/results_tab.txt", 'w')
table_file = [] #open(out_path+"/results_table.txt", 'w')
txt_file = [] #open(out_path+"/results_table.txt", 'w')
ref_file = [] #open(out_path+"/Resistance_gene_seq.fsa", 'w')
hit_file = [] #open(out_path+"/Hit_in_genome_seq.fsa", 'w')

# Write the header for the tab file
tab_file.append("Toxin_gene\tIdentity\tHSP/Reference\tContig\tPosition_in_contig\tPhenotype\tAccession_no.")

# Getting and writing out the results
for db in results:
   if db == "excluded":
      continue
   profile = db[0].upper() + db[1:].lower()
   if results[db] == "No hit found":
      table_file.append("%s\n%s\n"%(profile, results[db]))
   else:
      table_file.append("%s"%(profile))
      table_file.append("Toxin_gene\tIdentity\tHSP/Reference\tContig\tPosition_in_contig\tPhenotype\tAccession_no.")
      for hit in results[db]:
         res_header = results[db][hit]["sbjct_header"]
         tmp = res_header.split("_")
         gene = tmp[0]
         acc = tmp[2]
         ID = float(results[db][hit]["perc_ident"])
         ID = round(ID, 2)
         sbjt_length = results[db][hit]["sbjct_length"]
         HSP = results[db][hit]["HSP_length"]
         contig_name = results[db][hit]["contig_name"]
         positions = "%s..%s"%(results[db][hit]["query_start"], results[db][hit]["query_end"])
         pheno = phenos[gene][0]
         acc = phenos[gene][1]         
         # Write tabels
         table_file.append("%s\t%s\t%s/%s\t%s\t%s\t%s\t%s"%(gene, ID, HSP, sbjt_length, contig_name, positions, pheno, acc))
         tab_file.append("%s\t%s\t%s/%s\t%s\t%s\t%s\t%s"%(gene, ID, HSP, sbjt_length, contig_name, positions, db.upper() + " (" + pheno +")", acc))
         
         # Writing subjet/ref sequence
         ref_seq = sbjct_align[db][hit]
         ref_file.append(">%s"%(gene))
         for i in range(0, len(ref_seq), 60):
            ref_file.append("%s"%(ref_seq[i:i + 60]))
         
         # Writing query/hit sequence
         hit_seq = query_align[db][hit]
         #>aac(2')-Ic: PERFECT MATCH, ID: 100.00%, HSP/Length: 546/546, Positions in reference: 1..546, Contig name: gi|375294201|ref|NC_016768.1|, Position: 314249..314794
         sbjct_start = results[db][hit]["sbjct_start"]
         sbjct_end = results[db][hit]["sbjct_end"]
         hit_file.append(">%s, ID: %s %%, HSP/Length: %s/%s, Positions in reference: %s..%s, Contig name: %s, Position: %s"%(gene, ID, HSP, sbjt_length, sbjct_start, sbjct_end, contig_name, positions))

         for i in range(0, len(hit_seq), 60):
            hit_file.append("%s"%(hit_seq[i:i + 60]))
      table_file.append("")

#print("\n".join(table_file))
#print("\n".join(tab_file))

#print("tab file")
with open(out_path+"/results_tab.txt", 'w') as tabfile:
    tabfile.write("\n".join(tab_file))
#print("**********")
#print("table file")
with open(out_path+"/results_table.txt", 'w') as tablefile:
    tablefile.write("\n".join(table_file))  
table_file = [] #open(out_path+"/results_table.txt", 'w')
txt_file = [] #open(out_path+"/results_table.txt", 'w')
ref_file = [] #open(out_path+"/Resistance_gene_seq.fsa", 'w')
hit_file = [] #open(out_path+"/Hit_in_genome_seq.fsa", 'w')


# Write header to output file
"""
for i in range(len(gene_list)):
    gene = gene_list[i]
    print(gene)
    for hit in results[locus]:
        allel_hit = results[locus][hit]['sbjct_header']
        allel = allel_hit.split("_")[-1]
        coverage  = results[locus][hit]['perc_coverage']
        evalue    = results[locus][hit]['evalue']
        identity  = results[locus][hit]['perc_ident']
        print("{}\t{}\t{}\t{}".format(allel_hit, coverage, identity, evalue))
        found_loci.append(locus)
        if identity == 100:
            if locus in best_alleles:
                if best_alleles[locus][1] < coverage:
                    best_alleles[locus] = [allel, coverage, evalue]
                elif best_alleles[locus][1] == coverage and coverage == 100:
                    # Test for different in evalue
                    evalue_before = best_alleles[locus][2]
                    if evalue > evalue_before:
                        small_evalue = evalue_before
                        large_evalue = evalue
                        best_alleles[locus] = [allel, coverage, evalue]
                    else:
                        small_evalue = evalue
                        large_evalue = evalue_before
                    try:
                        ratio =  small_evalue / large_evalue
                    except ZeroDivisionError:
                        print("Unknown ratio, set to 1 {}, {}".format(large_evalue, small_evalue))
                        ratio = 1
                    print(ratio) 
                    if ratio < 2:
                        best_alleles[locus] = [-1, coverage, evalue]
                        print(best_alleles[locus])
            else:
                best_alleles[locus] = [allel, coverage, evalue]

    print("")

# Get called alleles
allel_str = filename
for locus in gene_list:
    locus = locus.strip()
    if locus in best_alleles:
        allel_str += "\t%s" %(best_alleles[locus][0])
    elif locus in found_loci:
        allel_str += "\tNaN"
    else:
        allel_str += "\tN"
allel_output += [allel_str]
#for line in allel_output:
#    print(line)
#print(allel_output)


# Load ST-dict pickle
pickle_path = specie_path + "/%s_profile.p"%(species)

T0 = time.time()
if os.path.isfile(pickle_path):
    try:
        loci_allel_dict = pickle.load(open(pickle_path, "rb"))
        T1 = time.time()
        print("pickle_loaded: %d s"%(int(T1-T0)) )
    except IOError:
        sys.stdout.write("Error, pickle not found", pickle_path)
        quit(1)

    st_filename = args.out_path + "-st.txt"
    #if os.path.isfile(pickle_path):
    # Write header in output file
    st_output += st_typing(loci_allel_dict, allel_output)

    # Write ST-type output
    with open(st_filename, "w") as fh:
        fh.write(st_output)
    print(st_output)
with open(args.out_path + ".txt", "w") as fh:
    fh.write("\n".join(allel_output) + "\n")
"""
