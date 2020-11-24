#!/usr/bin/python
import argparse
from Bio import SeqIO

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Append newly sequenced genomes to current genome dataset, and export metadata",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--genomes", required=True, help="FASTA file with latest genomes from GISAID")
    parser.add_argument("--keep", required=True, help="TXT file with accession number of genomes to be included")
    parser.add_argument("--remove", required=True, help="TXT file with accession number of genomes to be removed")
    parser.add_argument("--output1", required=True, help="FASTA file containing filtered sequences")
    parser.add_argument("--output2", required=True, help="Renaming file")

    args = parser.parse_args()

    genomes = args.genomes
    keep = args.keep
    remove = args.remove
    outfile1 = args.output1
    outfile2 = args.output2


    # genomes = path + "gisaid_hcov-19.fasta"
    # new_genomes = path + "new_genomes.fasta"
    # keep = path + 'keep.txt'
    # remove = path + "remove.txt"
    # outfile = path + "temp_sequences.fasta"



    # create a list of the existing sequences
    all_sequences = {}
    for fasta in SeqIO.parse(open(genomes),'fasta'):
        id, seq = fasta.description, fasta.seq
        id = id.replace('hCoV-19/', '').split('|')[0].replace(' ', '')
        all_sequences[id] = str(seq)

    # create a list of sequences to be added in all instances
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
    remove_sequences = []
    for id in open(remove, "r").readlines():
        if id[0] not in ["#", "\n"]:
            id = id.strip()
            remove_sequences.append(id)


    # export only sequences to be used in the nextstrain build
    c = 1
    sequences = {**keep_sequences}
    print('\n### Exporting sequences\n')
    exported = []
    output1 = open(outfile1, 'w')
    output2 = open(outfile2, 'w')
    for id in sequences.keys():
        if id not in remove_sequences: # filter out unwanted sequences
            entry = ">" + id + "\n" + sequences[id].upper() + "\n"
            exported.append(id)
            output1.write(entry)
            new_id = id.replace('/', '_')
            output2.write(new_id + '\t' + id + '\n')
            print(str(c) + '. ' + id)
        else:
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
    print('\n### Excluding sequences ###\n')
    e = 1
    for id in remove_sequences:
        print(str(e) + '. ' + id)
        e += 1

    print('\n### Final result\n')

    print('GISAID file contains ' + str(len(all_sequences)) + ' sequences\n')

    print(str(len(mismatch)) + ' genomes in keep.txt were NOT FOUND on GISAID database')
    print(str(len(keep_sequences)) + ' genomes ADDED from GISAID dataset')
    print(str(len(remove_sequences)) + ' genomes were REMOVED according to remove.txt\n')
    print(str(len(exported)) + ' genomes included in FINAL dataset\n')