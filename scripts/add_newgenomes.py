#!/usr/bin/python
import argparse
from Bio import SeqIO

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Append newly sequenced genomes to current genome dataset, and export metadata",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--genomes", required=True, help="FASTA file with latest genomes from GISAID")
    parser.add_argument("--new-genomes", required=True, help="FASTA file with newly sequenced genomes")
    parser.add_argument("--keep", required=True, help="TXT file with accession number of genomes to be included")
    parser.add_argument("--remove", required=True, help="TXT file with accession number of genomes to be removed")
    parser.add_argument("--output", required=True, help="FASTA file containing filtered sequences")
    args = parser.parse_args()

    genomes = args.genomes
    new_genomes = args.new_genomes
    keep = args.keep
    remove = args.remove
    outfile = args.output


    # genomes = path + "gisaid_hcov-19.fasta"
    # new_genomes = path + "new_genomes.fasta"
    # keep = path + 'keep.txt'
    # remove = path + "remove.txt"
    # outfile = path + "temp_sequences.fasta"
    # # outfile2 = path + "rename.tsv"


    # create a list of the existing sequences
    print('\n### Processing existing genomes...\n')
    all_sequences = {}
    for fasta in SeqIO.parse(open(genomes),'fasta'):
        id, seq = fasta.description, fasta.seq
        id = id.replace('hCoV-19/', '').split('|')[0].replace(' ', '')
        if id not in all_sequences: # avoid potential duplicates
            all_sequences[id] = str(seq)

    # store only new sequences in a dictionary, ignoring existing ones
    print('### Processing newly sequenced genomes...\n')
    newly_sequenced = {}
    for fasta in SeqIO.parse(open(new_genomes),'fasta'):
        id, seq = fasta.description, fasta.seq
        if id not in newly_sequenced.keys(): # avoid potential duplicates
            newly_sequenced[id] = str(seq)

    # create a list of sequences to be added in all instances
    print('### Searching for pre-selected genomes to be added...\n')
    keep_sequences = {}
    mismatch = []
    for id in sorted(open(keep, "r").readlines()):
        if id[0] not in ["#", "\n"]:
            id = id.strip()
            if id not in keep_sequences.keys():
                try:
                    keep_sequences[id] = all_sequences[id]
                except:
                    mismatch.append(id)

    # create a list of sequences to be ignored in all instances
    print('### Creating list of genomes to be ignored...\n')
    remove_sequences = []
    for id in open(remove, "r").readlines():
        if id[0] not in ["#", "\n"]:
            id = id.strip()
            remove_sequences.append(id)


    # export only sequences to be used in the nextstrain build
    c = 1
    sequences = {**keep_sequences, **newly_sequenced}
    print('### Exporting sequences\n')
    exported = []
    removed = []
    output = open(outfile, 'w')
    for id in sequences.keys():
        if id not in remove_sequences: # filter out unwanted sequences
            entry = ">" + id + "\n" + sequences[id].upper() + "\n"
            exported.append(id)
            output.write(entry)
            print(str(c) + '. ' + id)
        else:
            removed.append(id)
            c -= 1
        c += 1


    # mismatched sequence headers
    print('\n### Possible sequence header mismatches\n')
    m = 1
    for id in mismatch:
        print(str(m) + '. ' + id)
        m += 1
    if len(mismatch) < 1:
        print('\tNo mismatches found...')


    # excluding sequences
    if len(removed) > 0:
        print('\n### Excluding sequences ###\n')
        e = 1
        for id in removed:
            print(str(e) + '. ' + id)
            e += 1

    print('\n### Final result\n')

    print('GISAID file contains ' + str(len(all_sequences)) + ' sequences\n')

    print(str(len(mismatch)) + ' genomes in keep.txt were NOT FOUND on GISAID database')
    print(str(len(keep_sequences)) + ' genomes ADDED from GISAID dataset')
    print(str(len(newly_sequenced)) + ' newly sequenced genomes were added')
    print(str(len(removed)) + ' genomes were REMOVED according to remove.txt\n')
    print(str(len(exported)) + ' genomes included in FINAL dataset\n')