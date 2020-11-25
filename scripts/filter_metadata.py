# -*- coding: utf-8 -*-

import pycountry_convert as pyCountry
import pycountry
import argparse
from Bio import SeqIO
import pandas as pd
import numpy as np

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Filter nextstrain metadata files re-formmating and exporting only selected lines",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--genomes", required=True, help="FASTA file genomes to be used")
    parser.add_argument("--metadata1", required=True, help="Metadata file from NextStrain")
    parser.add_argument("--metadata2", required=False, help="Custom lab metadata file")
    parser.add_argument("--output1", required=True, help="Filtered metadata file")
    parser.add_argument("--output2", required=True, help="Reformatted, final FASTA file")
    parser.add_argument("--output3", required=False, help="IQTree renaming list")
    args = parser.parse_args()

    genomes = args.genomes
    metadata1 = args.metadata1
    metadata2 = args.metadata2
    output1 = args.output1
    output2 = args.output2
    output3 = args.output3

    # genomes = path + 'temp_sequences.fasta'
    # metadata1 = path + 'metadata_nextstrain.tsv'
    # metadata2 = path + 'COVID-19_sequencing.xlsx'
    # output1 = path + 'metadata_filtered.tsv'
    # output2 = path + 'sequences.fasta'
    # output3 = path + "rename.tsv"

    # create a dict of existing sequences
    sequences = {}
    for fasta in SeqIO.parse(open(genomes), 'fasta'):  # as fasta:
        id, seq = fasta.description, fasta.seq
        if id not in sequences.keys():
            sequences[id] = str(seq)

    # get ISO alpha3 country codes
    isos = {}
    def get_iso(country):
        global isos
        if country not in isos.keys():
            try:
                isoCode = pyCountry.country_name_to_country_alpha3(country, cn_name_format="default")
                isos[country] = isoCode
            except:
                try:
                    isoCode = pycountry.countries.search_fuzzy(country)[0].alpha_3
                    isos[country] = isoCode
                except:
                    isos[country] = ''
        return isos[country]


    # nextstrain metadata
    dfN = pd.read_csv(metadata1, encoding='utf-8', sep='\t', dtype='str')
    try:
        dfN.insert(4, 'iso', '')
    except:
        pass
    dfN['update'] = ''
    dfN.fillna('', inplace=True)
    list_columns = dfN.columns.values  # list of column in the original metadata file
    # print(dfN)

    # Lab genomes metadata
    dfL = pd.read_excel(metadata2, index_col=None, header=0, sheet_name=0,
                        converters={'Sample-ID': str, 'Collection-date': str, 'Update': str})
    dfL.fillna('', inplace=True)

    # dfL = dfL.rename(columns={'Sample-ID': 'strain', 'colB': 'columnB', 'colD': 'columnD', 'colE': 'columnE', 'colF': 'columnF'})
    dfL = dfL.rename(columns={'Sample-ID': 'id', 'Collection-date': 'date', 'Country': 'country', 'Division': 'division',
                              'State': 'code', 'Location': 'location', 'Lineage': 'pangolin_lineage', 'Source': 'originating_lab',
                              'Update': 'update'})

    # add inexistent columns
    for col in list_columns:
        if col not in dfL.columns:
            dfL[col] = ''

    # output dataframe
    outputDF = pd.DataFrame(columns=list_columns)
    found = []
    lab_label = {}

    # process metadata from excel sheet
    for idx, row in dfL.iterrows():
        id = dfL.loc[idx, 'id']
        if id in sequences:
            dict_row = {}
            for col in list_columns:
                dict_row[col] = ''
                if col in row:
                    # print(col)
                    dict_row[col] = dfL.loc[idx, col]

            # fix strain name
            if dfL.loc[idx, 'code'] == '':
                code = 'XX'
            else:
                code = dfL.loc[idx, 'code']
            strain = 'USA/' + code + '-' + dfL.loc[idx, 'id'] + '/2020' # set the strain name
            dict_row['strain'] = strain

            if len(str(dict_row['date'])) > 1:
                dict_row['date'] = dict_row['date'].split(' ')[0].replace('.', '-').replace('/', '-')

            # fix exposure
            columns_exposure = ['country_exposure', 'division_exposure']
            for level_exposure in columns_exposure:
                level = level_exposure.split('_')[0]
                if dict_row[level_exposure] in ['', None]:
                    dict_row[level_exposure] = dict_row[level]

            dict_row['iso'] = get_iso(dict_row['country_exposure'])
            dict_row['submitting_lab'] = 'Grubaugh Lab - Yale School of Public Health'
            dict_row['authors'] = 'Fauver et al'
            dict_row['update'] = 'Update' + str('0' * (2 - len(dfL.loc[idx, 'update']))) + dfL.loc[idx, 'update']

            found.append(strain)
            lab_label[id] = strain

            outputDF = outputDF.append(dict_row, ignore_index=True)


    # process metadata from TSV
    dfN = dfN[dfN['strain'].isin(sequences.keys())]
    for idx, row in dfN.iterrows():
        strain = dfN.loc[idx, 'strain']
        if strain in sequences:
            if strain in outputDF['strain'].to_list(): # skip metadata line if already sampled from lab metadata
                continue
            dict_row = {}
            for col in list_columns:
                dict_row[col] = ''
                if col in row:
                    dict_row[col] = dfN.loc[idx, col]

            # fix exposure
            columns_exposure = ['country_exposure', 'division_exposure']
            for level_exposure in columns_exposure:
                level = level_exposure.split('_')[0]
                if dict_row[level_exposure] in ['', None]:
                    dict_row[level_exposure] = dict_row[level]

            dict_row['iso'] = get_iso(dict_row['country_exposure'])
            found.append(strain)

            outputDF = outputDF.append(dict_row, ignore_index=True)


    # write new metadata files
    outputDF = outputDF.drop(columns=['region', 'region_exposure'])
    outputDF.to_csv(output1, sep='\t', index=False)


    # write sequence file
    exported = []
    # with open(output2, 'w') as outfile2:
    outfile2 = open(output2, 'w')
    outfile3 = open(output3, 'w')

    # export new metadata lines
    for id, sequence in sequences.items():
        if id in lab_label:
            if lab_label[id] not in exported:
                entry = '>' + lab_label[id] + '\n' + sequence + '\n'
                outfile2.write(entry)
                print('* Exporting newly sequenced genome and metadata for ' + id)
                new_id = lab_label[id].replace('/', '_')
                outfile3.write(new_id + '\t' + lab_label[id] + '\n')
                exported.append(lab_label[id])
        else:
            if id not in exported:
                entry = '>' + id + '\n' + sequence + '\n'
                outfile2.write(entry)
                new_id = id.replace('/', '_')
                outfile3.write(new_id + '\t' + id + '\n')
                exported.append(id)

print('\nMetadata file successfully reformatted and exported!\n')
