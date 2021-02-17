#!/usr/bin/python

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Created by: Anderson Brito
# Email: andersonfbrito@gmail.com
#
# Release date: 2021-20-16
# Last update: 2021-20-17
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import argparse
import json

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Convert GISAID JSON file into a fasta sequence file, filtering by genome coverage",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--input", required=True, help="FASTA file with genome alignment")
    parser.add_argument("--max-missing", required=False, type=int,  default='30', help="Maximum percentage of Ns or gaps (int: 1-100)")
    parser.add_argument("--output", required=True, help="TSV file with variants")
    args = parser.parse_args()

    json_file = args.input
    genome_size = 29420
    max_gaps = args.max_missing
    min_size = genome_size - int(genome_size * max_gaps/100)
    output = args.output


    print('\n### Exporting sequences\n')
    outfile = open(output, 'w')
    with open(json_file) as infile:
        c = 0
        for line in infile:
            entry = json.loads(line)
            # print(entry)
            id = entry['covv_virus_name'].replace('hCoV-19/', '').replace(' ', '')
            seq = entry['sequence'].replace('\n','').replace('N','')
            if int(len(seq)) > min_size:
                c += 1
                print(str(c) + '. ' + id)
                entry = ">" + id + "\n" + seq.upper() + "\n"
                outfile.write(entry)
            else:
                print('- ' + id + ' is too small. Coverage = ' + str(round(float(len(seq)/genome_size) * 100, 2)) + '%')

    outfile.close()
    print('\n' + str(c) + ' genomes successfully saved in ' + output + ' (coverage >= ' + str(100 - max_gaps) + '%)\n')
