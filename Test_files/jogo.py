#!/usr/bin/env python

# python Laboratorio_Final.py -o "Psammodromus algirus" -l "['Psammodromus microdactylus', 'Psammodromus occidentalis','Psammodromus blanci', 'Psammodromus edwarsianus', 'Psammodromus hispanicus', 'Psammodromus algirus']"

from Bio import Entrez
from Bio import SeqIO
import time
import pandas as pd
import argparse


def get_sequence_from_accession_number(gene):
    list_syn = []
    hand = Entrez.esearch(db="gene", term=gene + "[sym]", usehistory="n", retmax=10000)
    rec = Entrez.read(hand)
    hand.close()
    with Entrez.efetch(db="gene", retmode="text", usehistory="y", webenv=rec["WebEnv"],
                       query_key=rec["QueryKey"]) as handle:
        sequence_fasta = str(handle.read())[1:][:-2]
        # print(sequence_fasta)
        if "<ERROR>" not in sequence_fasta:
            for piece_data in sequence_fasta.split("\n\n"):
                row_piecedata = piece_data.split("\n")
                if row_piecedata[0] == "":
                    row_piecedata.pop(0)

                if row_piecedata[0] != "":
                    list_syn.append(row_piecedata[0].split(". ")[1])
    # print(list_syn)
    name = ""
    counter = 0
    for gene_name in set(list_syn):
        if list_syn.count(gene_name) > counter:
            print("---")
            counter = list_syn.count(gene_name)
            name = gene_name
    print(counter, name)
    return name.lower()


# Parte 1
# recebe strig e returna set
def get_sequence_from_accession_number_i(gene):
    list_syn = []
    hand = Entrez.esearch(db="gene", term=gene + "[sym]", usehistory="n", retmax=10000)
    rec = Entrez.read(hand)
    hand.close()
    with Entrez.efetch(db="gene", retmode="text", usehistory="y", webenv=rec["WebEnv"],
                       query_key=rec["QueryKey"]) as handle:
        sequence_fasta = str(handle.read())[1:][:-2]
        # print(sequence_fasta)
        if "<ERROR>" not in sequence_fasta:
            for piece_data in sequence_fasta.split("\n\n"):
                row_piecedata = piece_data.split("\n")
                if row_piecedata[0] == "":
                    row_piecedata.pop(0)

                if row_piecedata[0] != "":
                    list_syn.append(row_piecedata[0].split(". ")[1])
        for i in set(list_syn):
            print()
        # print(set(list_syn))
        return set(list_syn)


def menor_string(string1, string2):
    if len(string1) < len(string2):
        return string1, string2
    else:
        return string2, string1


def reduzir_keys(new_di):
    dict_unique = {}
    for k, v in new_di.items():
        for key, value in new_di.items():
            if k != key and v == value:
                menor_str, maior_str = menor_string(k, key)
                if menor_str in dict_unique.keys() or maior_str in dict_unique.keys():
                    pass
                else:
                    dict_unique[menor_string(k, key)[0]] = value
    return dict_unique


def unir_set_com_key_iguais(di_set):
    # apenas um iteracao para cada dicionario

    new_di = {}
    # print("---------------------------------------------------------------------------")
    # print("di_set.keys()", di_set.keys())
    for k, v in di_set.items():
        for key, value in di_set.items():
            if k != key:
                if len(v.intersection(value)) >= 1:
                    # print("KEYYYYYSSSS===", k, key, "----->  v.intersection(value)===", v.intersection(value))
                    new_di[k] = value.union(v)
            else:
                new_di[k] = value

    print(new_di.keys())
    dict_unique = new_di.copy()
    l = []

    for k, v in new_di.items():
        for key, value in new_di.items():
            if k != key and v == value:
                if k not in l and k in dict_unique.keys():
                    l.append(key)

                    del dict_unique[k]

    return list(dict_unique.keys())

# Parte 3
def get_sequence(acc_number, beg_sequence, end_sequence):
    if beg_sequence != 0:
        beg_sequence = beg_sequence - 1

    hand = Entrez.esearch(db="nucleotide", term=acc_number, usehistory="y")
    rec = Entrez.read(hand)
    hand.close()

    with Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", usehistory="y", webenv=rec["WebEnv"],
                       query_key=rec["QueryKey"]) as handle:
        seq_record = handle.read()
        nome = "_".join(seq_record.split("\n")[0].split()[1:])
        sequence = "".join(seq_record.split("\n")[1:])
        acc_number.replace(" ", "")
        return ">" + acc_number + "_" + str(beg_sequence) + "_" + str(end_sequence), nome, sequence[
                                                                                           beg_sequence:end_sequence]


# ##########################################################

# Parte 1
def add_gene(ind, dict_gene_local, list_names_features, seq_record):
    index = list_names_features.index(ind)
    dict_qualifiers = seq_record.features[index].qualifiers
    try:
        gene = dict_qualifiers["gene"]
    except:
        try:
            gene = dict_qualifiers["product"]  # + "p"
        except:
            try:
                gene = dict_qualifiers["locus_tag"]
            except:
                # print("ERRO no TRY: (dict_qualifiers)")
                # gene = ""  # ainda por mudificar
                return dict_gene_local
    print(gene[0])
    # gene_lower = get_sequence_from_accession_number(gene[0])
    # print("-------------------------------------------------------------------")
    # print("----------------------->", gene_lower)

    if gene[0] not in list(dict_gene_local.keys()):
        dict_gene_local[gene[0]] = []
    location = seq_record.features[index].location

    start_local = int(str(location.start).replace("<", ""))
    end_local = int(str(location.end).replace(">", ""))

    # print(start_local, end_local)
    dict_gene_local[gene[0]].append([seq_record.id, start_local, end_local])

    return dict_gene_local


def get_all_gene_organism(organism_name, option):
    hand = Entrez.esearch(db="nucleotide", term=organism_name + "[organism]", usehistory="y")
    rec = Entrez.read(hand)
    hand.close()
    limit = 100
    dict_gene_local = {}

    for start in range(0, int(rec["Count"]), limit):
        with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", usehistory="y", retstart=start,
                           retmax=limit, webenv=rec["WebEnv"], query_key=rec["QueryKey"]) as handle:

            for seq_record in SeqIO.parse(handle, "gb"):
                list_names_features = [seq_record.features[j].type
                                       for j in range(len(seq_record.features))]

                n_gene = list_names_features.count("gene")
                if n_gene == 0 or n_gene == 1:
                    if "gene" in list_names_features:
                        dict_gene_local = add_gene("gene", dict_gene_local, list_names_features, seq_record)
                    elif "CDS" in list_names_features:
                        dict_gene_local = add_gene("CDS", dict_gene_local, list_names_features, seq_record)
                    elif "mRNA" in list_names_features:
                        dict_gene_local = add_gene("mRNA", dict_gene_local, list_names_features, seq_record)
                    elif "rRNA" in list_names_features:
                        dict_gene_local = add_gene("rRNA", dict_gene_local, list_names_features, seq_record)
                    elif "tRNA" in list_names_features:
                        dict_gene_local = add_gene("tRNA", dict_gene_local, list_names_features, seq_record)
                    elif "misc_feature" in list_names_features:
                        dict_gene_local = add_gene("misc_feature", dict_gene_local, list_names_features, seq_record)
                    else:
                        pass

                if n_gene > 1:
                    list_t = [indice for indice in range(len(list_names_features))
                              if list_names_features[indice] == "gene"]

                    for index in list_t:
                        dict_gene_local = add_gene(list_names_features[index], dict_gene_local, list_names_features,
                                                   seq_record)

    def escolher_sequencia_b_s(dict_gene_local, location):
        dict_b_s = {}

        for key in dict_gene_local.keys():
            list_b_s = []
            for value in dict_gene_local[key]:
                start_local = value[1]

                end_local = value[2]
                list_b_s.append(end_local - start_local)
            index_p = list_b_s.index(sorted(list_b_s)[location])
            dict_b_s[key] = dict_gene_local[key][index_p]

        return dict_b_s

    def escolher_sequencia_f_l(dict_gene_local, local_list):

        dict_new = {}
        for key in dict_gene_local.keys():
            dict_new[key] = dict_gene_local[key][local_list]

        return dict_new

    if option == "b":
        return escolher_sequencia_b_s(dict_gene_local, 0)
    elif option == "s":
        return escolher_sequencia_b_s(dict_gene_local, -1)
    elif option == "f":
        return escolher_sequencia_f_l(dict_gene_local, 0)
    else:  # option == "-l":
        return escolher_sequencia_f_l(dict_gene_local, -1)


# Parte 2
def delete_missing_data(dataframe, percentage_org_gene):
    for name in dataframe.columns:
        percentage_null = dataframe[name].isnull().sum() / len(dataframe[name])
        if percentage_org_gene < percentage_null:
            del dataframe[name]
    return dataframe


if __name__ == '__main__':
    Entrez.email = "rubenaalbquerque@hotmail.com"
    time_a = time.time()
    parser = argparse.ArgumentParser(description='Process some integers.')

    parser.add_argument('--organism', '-o', type=str, help='nome do organismo')
    parser.add_argument('--list_org', '-l', type=str, help='nome do organismo')
    parser.add_argument('--proportion_org', '-p_o', type=float, help='proporcao de missing organismo')
    parser.add_argument('--proportion_gene', '-p_g', type=float, help='proporcao de missing gene')
    parser.add_argument('--option_sequence', '-o_s', type=str, help='Opcao')
    args = parser.parse_args()
    list_with_organism_names = eval(args.list_org)

    # Parte 1
    dict_gene_unique = {}
    dict_gene_alias = {}  # dicionario com o nome do gene e os respetivos alias, designacao...
    new_dict_name = {}
    dict_name = {}

    for name in list_with_organism_names:
        print(name.upper())
        dict_name[name] = get_all_gene_organism(name, option=args.option_sequence)

    print("----------------------------------------------------------------")
    print(dict_name)


    def get_sequence_from_accession_number(gene):
        list_syn = []
        hand = Entrez.esearch(db="gene", term=gene + "[sym]", usehistory="n", retmax=10000)
        rec = Entrez.read(hand)
        hand.close()
        with Entrez.efetch(db="gene", retmode="text", usehistory="y", webenv=rec["WebEnv"], query_key=rec["QueryKey"]) as handle:
            sequence_fasta = str(handle.read())[1:][:-2]
            # print(sequence_fasta)
            if "<ERROR>" not in sequence_fasta:
                for piece_data in sequence_fasta.split("\n\n"):
                    row_piecedata = piece_data.split("\n")
                    if row_piecedata[0] == "":
                        row_piecedata.pop(0)

                    if row_piecedata[0] != "":
                        list_syn.append(row_piecedata[0].split(". ")[1])
        # print(list_syn)
        name = ""
        counter = 0
        for gene_name in set(list_syn):
            if list_syn.count(gene_name) > counter:
                counter = list_syn.count(gene_name)
                name = gene_name
        # print(counter, name)
        return name.lower()
    d = {}
    for name_org in list_with_organism_names:
        for name_gene in dict_name[name_org].keys():
            print("--->", name_gene)
            nam = get_sequence_from_accession_number(name_gene)
            print(nam)
            d[nam] = dict_name[name_org][name_gene]

    print("----------------------------------------------------------------")
    print("\n\n\n", d, "\n\n")

    print("-----------------------------------------------------------------")

    # Parte 2
    df = pd.DataFrame(dict_name)
    df = delete_missing_data(df, args.proportion_gene)
    dataframe = delete_missing_data(df.T, args.proportion_gene)

    print(dataframe.T)

    # Parte 3
    for col in dataframe.columns:
        for line in list(dataframe[col]):
            if str(line) != str("nan"):
                accession, name, sequence = get_sequence(line[0], int(line[1]), int(line[2]))
                if str(sequence) != "":
                    with open(col + ".fasta", "a") as new_file:
                        new_file.write(accession + "\n" + sequence + "\n\n")

                    with open("Accession_numbers__name.txt", "a") as file:
                        file.write(accession[1:] + "__" + name + "\n")

    print("Program in", time.time() - time_a)
    print("Minutes:", (time.time() - time_a) / 60)
