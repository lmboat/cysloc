# Created Date: 201116
# Modified Date: 230110
# Author: Lisa Boatner
# Technical Updates: LFQ ratio quant multiplexed peptides.
# Notes: Assumes numbering is off for mulitplexed labeled peptides, therefore this version adds +10 on the second amino acid postition for mouse fasta mapping.

import os, sys, csv, re, argparse, re, collections, glob, statistics, numbers, string, time, math
import numpy as np
from os import listdir, path
from os.path import isfile, isdir, join
from Bio import SeqIO
# from quoters import Quote

class UniprotDatabase:
  def __init__ (self, args, cd):
    self.game_dir = cd
    self.uniprot_dir = args.uf
    self.uniprot_file = args.uniref
    self.uniprot_names = []
    self.uniprot_proteins = []

    self.run()

  def run(self):
    os.chdir(self.uniprot_dir)
    self.read_file(self.uniprot_file)
    os.chdir(self.game_dir)

  # read Uniprot reviewed protein sequence database fasta file
  def read_file(self, file):
    self.check_file_exists(file)
    for seq_record in SeqIO.parse(file, "fasta"):
      # get string sequence
      protein_data = self.get_parsed_protein_data(seq_record.description)
      # create UniprotFastaProteins
      pro = UniprotFastaProtein(protein_data[0], protein_data[1], str(seq_record.seq), protein_data[2])
      self.uniprot_names.append(protein_data[0])
      self.uniprot_proteins.append(pro)

  # get protein id and gene name
  def get_parsed_protein_data(self, line):
    # remove un/reviewed prefixes
    if ("sp|" in line):
      line = line.replace("sp|", "").split("|")
    if ("tr|" in line):
      line = line.replace("tr|", "").split("|")  
    protein_id = line[0]
    name = line[1].split(" ")
    name = name[0]
    description = line[1]

    return [protein_id, name, description]

class UniprotFastaProtein:
  def __init__ (self, protein, name, sequence, description):
    self.protein = protein
    self.name = name
    self.sequence = sequence
    self.no_cysteine = False
    self.total_c_count = self.sequence.count("C")
    if self.sequence.count("C") == 0:
      self.no_cysteine = True
    self.description = description

  # convert data to string for writing csv
  def get_individual_data(self):
    return str(self.protein) + "," + str(self.name) + "," + str(self.sequence) + "," + str(self.no_cysteine) + "," + str(self.total_c_count)

class Game:
  def __init__ (self, args, cd, uniprot_db):
    self.args = args

    self.uniprot_db = uniprot_db

    self.game_dir = cd
    self.data_dir = cd + '/' + args.df
    self.data_folder = args.df
    self.results_dir = cd + '/' + args.rf
    self.results_folder = args.rf
    self.experiment_folders = []
    self.experiments = []
    self.unique_identifiers_dict = {}
    self.unique_proteins_dict = {}
    self.unique_peptides_dict = {}
    self.experiment_no = 0
    self.unique_identifiers_list = []
    self.unique_proteins_list = []
    self.unique_peptides_list = []

    self.total_proteins = []
    self.proteins_names = []
    self.total_identifiers = []
    self.total_peptides = []
    self.identifiers_names = []
    self.peptide_sequences = []
    self.result_files = []
    self.protein_comparison_dict, self.identifier_comparison_dict, self.peptide_comparison_dict = {}, {}, {}

    self.distinct_proteins_dict, self.distinct_identifiers_dict, self.distinct_peptides_dict = {}, {}, {}
    self.protein_counts_dict, self.identifier_counts_dict, self.peptide_counts_dict = {}, {}, {}

    self.protein_avgs, self.protein_stdevs = {}, {}
    self.identifier_avgs, self.identifier_stdevs = {}, {}
    self.peptide_avgs, self.peptide_stdevs = {}, {}
    self.replicate_folders = []

    if (args.entire_process == 'False') & (args.compile_experiments == 'False'):
      self.setup_game()
      self.run_game()
    elif (args.entire_process == 'False') & (args.compile_experiments == 'True'):
      self.finish_game()
    elif args.entire_process == 'True':
      self.setup_game()
      self.run_game()
      self.finish_game()
    else:
      print("I don't know what you want from me. Email lisaboatner@ucla.edu - I am a weeny, help me.")

  def setup_game(self):
    print('\n' + "On your mark... get set...")
    os.makedirs(self.results_dir, exist_ok=True)

    os.chdir(self.game_dir)
    os.chdir(self.data_folder)

    self.experiment_folders = [f for f in listdir(os.getcwd()) if isdir(join(os.getcwd(), f))]
    self.experiment_folders = sorted(self.experiment_folders)
    
    if (".DS_Store") in self.experiment_folders:
      self.experiment_folders.remove(".DS_Store")
    if (self.results_folder) in self.experiment_folders:
      self.experiment_folders.remove(self.results_folder)

    self.experiment_no = len(self.experiment_folders)

    self.empty_matrix = np.zeros((len(self.experiment_folders), len(self.experiment_folders)))
    self.protein_sim_matrix, self.identifier_sim_matrix, self.peptide_sim_matrix, self.protein_diff_matrix, self.identifier_diff_matrix, self.peptide_diff_matrix = np.copy(self.empty_matrix), np.copy(self.empty_matrix), np.copy(self.empty_matrix), np.copy(self.empty_matrix), np.copy(self.empty_matrix), np.copy(self.empty_matrix)

  def run_game(self):
    print('\n' + "GOOOOOOOO ┏(• ◡•)┛ RUN FOREST RUN ┏(• ◡•)┛ " + '\n')

    for i in range(len(self.experiment_folders)):
      os.chdir(self.data_dir)
      os.chdir(self.experiment_folders[i])
      experiment = Experiment(self.args, self.uniprot_db, self.experiment_folders[i], self.data_dir + '/' + self.experiment_folders[i], self.data_dir, self.results_dir)
      self.experiments.append(experiment)
      self.replicate_folders += experiment.replicate_folders
      self.unique_identifiers_dict[self.experiment_folders[i]] = experiment.experiment_identifiers
      self.unique_proteins_dict[self.experiment_folders[i]] = experiment.experiment_proteins
      self.unique_peptides_dict[self.experiment_folders[i]] = experiment.experiment_peptides
      self.protein_counts_dict.update(experiment.protein_counts)
      self.identifier_counts_dict.update(experiment.identifier_counts)
      self.peptide_counts_dict.update(experiment.peptide_counts)
      self.protein_avgs[self.experiment_folders[i]], self.protein_stdevs[self.experiment_folders[i]] = experiment.exp_protein_avg, experiment.exp_protein_stdev
      self.identifier_avgs[self.experiment_folders[i]], self.identifier_stdevs[self.experiment_folders[i]] = experiment.exp_identifier_avg, experiment.exp_identifier_stdev
      self.peptide_avgs[self.experiment_folders[i]], self.peptide_stdevs[self.experiment_folders[i]] = experiment.exp_peptide_avg, experiment.exp_peptide_stdev

      os.chdir(self.game_dir)

    self.compare_experiments()

  def finish_game(self):
    print('\n' + "Down to the wire, don't give up!" + '\n')

    os.chdir(self.results_dir)

    self.experiment_folders = [f for f in listdir(os.getcwd()) if isdir(join(os.getcwd(), f))]

    for i in range(len(self.experiment_folders)):
      os.chdir(self.experiment_folders[i])
      result_files = [f for f in listdir(os.getcwd()) if isfile(join(os.getcwd(), f))]
      self.result_files = self.result_files + result_files
      os.chdir(self.results_dir)

    if (self.args.cmc == "True"):
      print('\n' + "Move it or lose it." + '\n')

    for j in range(len(self.experiment_folders)):
      os.chdir(self.experiment_folders[j])
      result_files = [f for f in listdir(os.getcwd()) if isfile(join(os.getcwd(), f))]
      for k in range(len(result_files)):
        self.read_results_file(result_files[k], self.experiment_folders[j], self.experiment_folders)
      os.chdir(self.results_dir)

    experiments = ','.join([str(elem) for elem in self.experiment_folders]) 
    datasets = (','.join([str(elem) for elem in self.result_files])).replace('_results.csv', '') 
    experiment_counts = experiments.replace(',', '_experiment_count,') + '_experiment_count,'
    replicate_counts = datasets.replace(',', '_replicate_count,') + '_replicate_count'
    protein_header = "protein,description,gene,identifiers,peptides,peptide_count,spectral_count,modifications,modification_count,amino_acids,amino_acid_count,modification_masses,no_of_experiments_count,no_of_replicates_count," + experiment_counts + replicate_counts
    identifier_header = "identifier,protein,description,gene,peptides,peptide_count,spectral_count,modifications,modification_count,amino_acids,amino_acid_count,modification_masses,no_of_experiments_count,no_of_replicates_count," + experiment_counts + replicate_counts
    peptide_header = "peptide,unmodified_peptide,charge_state,probability,localization,proteins,identifiers,description,spectral_count,intensity,miscleavages,modifications,modification_count,amino_acids,amino_acid_count,modification_masses,no_of_experiments_count,no_of_replicates_count," + experiment_counts + replicate_counts
    experiment_medians = experiments.replace(',', '_avg_of_medians,') + '_avg_of_medians,'
    experiment_stdevs = experiments.replace(',', '_stdev_of_medians,') + '_stdev_of_medians,'
    replicate_medians = datasets.replace(',', '_median,') + '_median,'

    # if (self.args.experiment != "identification") & (self.args.experiment != "iodination"):
    if (self.args.experiment != "identification"):  
      experiment_medians = experiments.replace(',', '_avg_of_medians,') + '_avg_of_medians,'
      experiment_stdevs = experiments.replace(',', '_stdev_of_medians,') + '_stdev_of_medians,'
      protein_header += ',' + 'aggregate_avg_of_avg_of_medians,aggregate_stdev_of_avg_of_medians,' + experiment_medians + experiment_stdevs + 'total_quant_list,' + replicate_medians.replace('_avg_of_medians,', '_medians,')
      identifier_header += ',' + 'aggregate_avg_of_avg_of_medians,aggregate_stdev_of_avg_of_medians,' + experiment_medians + experiment_stdevs + 'total_quant_list,' + replicate_medians.replace('_avg_of_medians,', '_medians,')
      peptide_header += ',' + 'aggregate_avg_of_avg_of_medians,aggregate_stdev_of_avg_of_medians,' + experiment_medians + experiment_stdevs + 'total_quant_list,' + replicate_medians.replace('_avg_of_medians,', '_medians,')

    if ((self.args.c == "True") & (self.args.experiment == "isotop")):
      self.args.cpn = 1

    if (self.args.c == "True"):
      print('\n' + "Just one more." + '\n')
      new_proteins = self.get_new_team(self.total_proteins, True)
      new_identifiers = self.get_new_team(self.total_identifiers, False)
      new_peptides = self.get_new_team(self.total_peptides, False)
    else:
      new_proteins = self.total_proteins
      new_identifiers = self.total_identifiers
      new_peptides = self.total_peptides

    if self.args.write_outfile == 'True':
      self.write_file("compiled_proteins.csv", protein_header, new_proteins, True)  
      self.write_file("compiled_identifiers.csv", identifier_header, new_identifiers, True)
      self.write_file("compiled_peptides.csv", peptide_header, new_peptides, True)

    print('\n' + "Good game mate." + '\n')

  def check_file_exists(self, file):
    if (path.exists(file) == False):
      print(file + " does not exist in " + str(os.getcwd()))
      sys.exit()

  def get_new_team(self, current_team, peptide):
    new_team = []
    for i in range(len(current_team)):
      current_player = current_team[i]
      if (peptide == True):
        if (int(current_player.peptide_count) >= int(self.args.cpn)):
          new_team.append(current_player)
      else:
        if (int(current_player.experiment_reps) >= int(self.args.cen)):
          new_team.append(current_player)
    return new_team

  def read_results_file(self, file, experiment_id, experiments):
    self.check_file_exists(file)

    in_file = open(file, 'r')
    header_title = in_file.readline().strip()
    header = header_title.split(",")

    for line in in_file:
      line = line.strip().split(",")
      if ("contaminant" in line[header.index("Identifier")]):
        continue

      # TODO: Enzymatic termini
      if (self.args.cmc == "True"):
        if ((int(line[header.index("Number of Missed Cleavages")]) > int(self.args.cmcn)) == True):
          continue

      if (self.args.experiment == 'intensity'):
        col_name = 'Intensity'
      else:
        col_name = 'Area Ratio'

      self.get_proteins(line, header, file, experiment_id, experiments, col_name)
      self.get_identifiers(line, header, file, experiment_id, experiments, col_name)

      if (self.args.experiment == 'iodination'):
        self.get_iodination_peptides(line, header, file, experiment_id, experiments)
      else:
        self.get_peptides(line, header, file, experiment_id, experiments, col_name)

    in_file.close()

  def get_proteins(self, line, header, file, experiment_id, experiments, col_name):
    # protein, name, description, identifier, peptide, spectral_count, ratio, dataset, pep_files
    if line[header.index("Protein ID")] not in self.proteins_names:
      protein = Protein(line[header.index("Protein ID")], line[header.index("Protein Description")], line[header.index("Gene")], line[header.index("Identifier")], line[header.index("Peptide")], line[header.index("Spectral Count")], line[header.index("Modifications")], line[header.index("Amino Acids")], line[header.index("Modification Masses")], line[header.index(col_name)], file, self.result_files, self.args, experiment_id, experiments)
      self.total_proteins.append(protein)
      self.proteins_names.append(line[header.index("Protein ID")])
    # identifier, peptide, ratio, dataset
    else:
      idx = self.proteins_names.index(line[header.index("Protein ID")])
      current = self.total_proteins[idx]
      current.update(line[header.index("Identifier")], line[header.index("Peptide")], line[header.index("Spectral Count")], line[header.index("Modifications")], line[header.index("Amino Acids")], line[header.index("Modification Masses")], line[header.index(col_name)], file, experiment_id)

  def get_identifiers(self, line, header, file, experiment_id, experiments, col_name):
    # protein, description, identifier, peptide, ratio, dataset, pep_files
    if line[header.index("Identifier")] not in self.identifiers_names:
      identifier = Identifier(line[header.index("Protein ID")], line[header.index("Protein Description")], line[header.index("Identifier")], line[header.index("Gene")], line[header.index("Peptide")], line[header.index("Spectral Count")], line[header.index("Modifications")], line[header.index("Amino Acids")], line[header.index("Modification Masses")], line[header.index(col_name)], file, self.result_files, self.args, experiment_id, experiments)
      self.total_identifiers.append(identifier)
      self.identifiers_names.append(line[header.index("Identifier")])
    # peptide, ratio, dataset
    else:
      idx = self.identifiers_names.index(line[header.index("Identifier")])
      current = self.total_identifiers[idx]
      current.update(line[header.index("Peptide")], line[header.index("Spectral Count")], line[header.index("Modifications")], line[header.index("Amino Acids")], line[header.index("Modification Masses")], line[header.index(col_name)], file, experiment_id)

  def get_peptides(self, line, header, file, experiment_id, experiments, col_name):
    # peptide, charge_state, prob, localization, protein, identifier, description, spectral_count, dataset, pep_files
    if line[header.index("Peptide")] not in self.peptide_sequences:
      # Uncomment for updated results
      if (self.args.cmc == "False"):        
        peptide = Peptide(line[header.index("Peptide")], line[header.index("Charge")], line[header.index("Probability")], line[header.index("Best Localization")], line[header.index("Protein ID")], line[header.index("Identifier")], line[header.index("Protein Description")], line[header.index("Spectral Count")], 0, "", line[header.index("Modifications")], line[header.index("Amino Acids")], line[header.index("Modification Masses")], line[header.index(col_name)], file, self.result_files, self.args, experiment_id, experiments)
      else:
        peptide = Peptide(line[header.index("Peptide")], line[header.index("Charge")], line[header.index("Probability")], line[header.index("Best Localization")], line[header.index("Protein ID")], line[header.index("Identifier")], line[header.index("Protein Description")], line[header.index("Spectral Count")], 0, line[header.index("Number of Missed Cleavages")], line[header.index("Modifications")], line[header.index("Amino Acids")], line[header.index("Modification Masses")], line[header.index(col_name)], file, self.result_files, self.args, experiment_id, experiments)
      self.total_peptides.append(peptide)
      self.peptide_sequences.append(line[header.index("Peptide")])
    # charge_state, prob, protein, identifier, description, dataset
    else:
      idx = self.peptide_sequences.index(line[header.index("Peptide")])
      current = self.total_peptides[idx]
      current.update(line[header.index("Charge")], line[header.index("Probability")], line[header.index("Best Localization")], line[header.index("Protein ID")], line[header.index("Identifier")], line[header.index("Protein Description")], line[header.index("Spectral Count")], 0, line[header.index("Modifications")], line[header.index("Amino Acids")], line[header.index("Modification Masses")],  line[header.index(col_name)], file, experiment_id)

  def get_iodination_peptides(self, line, header, file, experiment_id, experiments):
    # peptide, charge_state, prob, localization, protein, identifier, description, spectral_count, dataset, pep_files
    identifier = line[header.index("Peptide")] + "_" + line[header.index("Charge")]
    if identifier not in self.peptide_sequences: 
      peptide = Peptide(line[header.index("Peptide")], line[header.index("Charge")], line[header.index("Probability")], line[header.index("Best Localization")], line[header.index("Protein ID")], line[header.index("Identifier")], line[header.index("Protein Description")], line[header.index("Spectral Count")], line[header.index("Intensity")], line[header.index("Number of Missed Cleavages")], line[header.index("Modifications")], line[header.index("Amino Acids")], line[header.index("Modification Masses")], line[header.index("Area Ratio")], file, self.result_files, self.args, experiment_id, experiments)
      self.total_peptides.append(peptide)
      self.peptide_sequences.append(identifier)
    else:
      idx = self.peptide_sequences.index(identifier)
      current = self.total_peptides[idx]
      current.update(line[header.index("Charge")], line[header.index("Probability")], line[header.index("Best Localization")], line[header.index("Protein ID")], line[header.index("Identifier")], line[header.index("Protein Description")], line[header.index("Spectral Count")], 0, line[header.index("Modifications")], line[header.index("Amino Acids")], line[header.index("Modification Masses")],  line[header.index("Area Ratio")], file, experiment_id)

  def compare_experiments(self):
    list_of_pairs = [(self.experiment_folders[p1], self.experiment_folders[p2]) for p1 in range(len(self.experiment_folders)) for p2 in range(p1+1,len(self.experiment_folders))]
    symbol = '|'

    for i in range(len(list_of_pairs)):
      current = list_of_pairs[i]

      protein_intersection = len(list(set(self.unique_proteins_dict[current[0]]).intersection(self.unique_proteins_dict[current[1]])))
      identifier_intersection = len(list(set(self.unique_identifiers_dict[current[0]]).intersection(self.unique_identifiers_dict[current[1]])))
      peptide_intersection = len(list(set(self.unique_peptides_dict[current[0]]).intersection(self.unique_peptides_dict[current[1]])))      
      self.protein_comparison_dict[current[0] + symbol + current[1]] = protein_intersection
      self.identifier_comparison_dict[current[0] + symbol + current[1]] = identifier_intersection
      self.peptide_comparison_dict[current[0] + symbol + current[1]] = peptide_intersection      

    self.unique_proteins_list = self.create_unique_list(self.unique_proteins_dict)
    self.unique_identifiers_list = self.create_unique_list(self.unique_identifiers_dict)
    self.unique_peptides_list = self.create_unique_list(self.unique_peptides_dict)

    self.compare_matrices(symbol)

    self.distinct_proteins_dict = self.create_total_list(self.unique_proteins_dict)
    self.distinct_identifiers_dict = self.create_total_list(self.unique_identifiers_dict)
    self.distinct_peptides_dict = self.create_total_list(self.unique_peptides_dict)

    os.chdir(self.results_dir)
    experiments = ','.join([str(elem) for elem in self.experiment_folders]) 
    replicates = ','.join([str(elem) for elem in self.replicate_folders]) 
    if self.args.write_outfile == "True":
      self.write_comparison_file(self.data_folder + '_comparison_results.csv', experiments, replicates)
    os.chdir(self.game_dir)

  def create_total_list(self, unique_dict):
    discinct_dict = {}
    for i in range(len(self.experiment_folders)):
      current = self.experiment_folders[i]
      mega_list = []
      for k in unique_dict:
        if k != current:
          mega_list = mega_list + unique_dict[k]
      distinct = len(unique_dict[current]) - len(set(mega_list).intersection(unique_dict[current]))
      discinct_dict[current] = distinct
    return discinct_dict

  def create_unique_list(self, unique_dict):
    unique_list = []
    for k in unique_dict:
      unique_list = list(set(unique_list + unique_dict[k]))
    return unique_list

  def update_matrix_diagonal(self, unique_dict, matrix):
    for k in unique_dict:
      matrix[self.experiment_folders.index(k), self.experiment_folders.index(k)] = len(unique_dict[k])
    return matrix

  def update_matrix_nondiagonal(self, compare_dict, unique_dict, matrix, similarity, symbol):
    for k in compare_dict:
      current = k.split(symbol)
      if similarity == True:
        matrix[self.experiment_folders.index(current[0]), self.experiment_folders.index(current[1])] = compare_dict[k]
        matrix[self.experiment_folders.index(current[1]), self.experiment_folders.index(current[0])] = compare_dict[k]
      else:
        matrix[self.experiment_folders.index(current[0]), self.experiment_folders.index(current[1])] = len(unique_dict[current[0]]) - compare_dict[k]
        matrix[self.experiment_folders.index(current[1]), self.experiment_folders.index(current[0])] = len(unique_dict[current[1]]) - compare_dict[k]
    return matrix    

  def compare_matrices(self, symbol):
    self.protein_sim_matrix = self.update_matrix_diagonal(self.unique_proteins_dict, self.protein_sim_matrix)
    self.protein_sim_matrix = self.update_matrix_nondiagonal(self.protein_comparison_dict, {}, self.protein_sim_matrix, True, symbol)
    self.identifier_sim_matrix = self.update_matrix_diagonal(self.unique_identifiers_dict, self.identifier_sim_matrix)
    self.identifier_sim_matrix = self.update_matrix_nondiagonal(self.identifier_comparison_dict, {}, self.identifier_sim_matrix, True, symbol)
    self.peptide_sim_matrix = self.update_matrix_diagonal(self.unique_peptides_dict, self.peptide_sim_matrix)
    self.peptide_sim_matrix = self.update_matrix_nondiagonal(self.peptide_comparison_dict, {}, self.peptide_sim_matrix, True, symbol)

    self.protein_diff_matrix = self.update_matrix_nondiagonal(self.protein_comparison_dict, self.unique_proteins_dict, self.protein_diff_matrix, False, symbol)
    self.identifier_diff_matrix = self.update_matrix_nondiagonal(self.identifier_comparison_dict, self.unique_identifiers_dict, self.identifier_diff_matrix, False, symbol)
    self.peptide_diff_matrix = self.update_matrix_nondiagonal(self.peptide_comparison_dict, self.unique_peptides_dict, self.peptide_diff_matrix, False, symbol)

  def write_comparison_file(self, file, experiments, replicates):
    out_file = open(file, 'w')

    header = ','.join([str(elem) for elem in self.experiment_folders])
    out_file.write('Similarity Protein Matrix,' + header + '\n')
    
    self.matrix_to_string(self.protein_sim_matrix, out_file)
    out_file.write('\n')
    
    out_file.write('Difference Protein Matrix,' + header + '\n')
    self.matrix_to_string(self.protein_diff_matrix, out_file)
    out_file.write('\n')
    
    out_file.write('Similarity Amino Acid Matrix,' + header + '\n')
    self.matrix_to_string(self.identifier_sim_matrix, out_file)
    out_file.write('\n')
    
    out_file.write('Difference Amino Acid Matrix,' + header + '\n')
    self.matrix_to_string(self.identifier_diff_matrix, out_file)
    out_file.write('\n')

    out_file.write('Similarity Peptide Matrix,' + header + '\n')
    self.matrix_to_string(self.peptide_sim_matrix, out_file)
    out_file.write('\n')
    
    out_file.write('Difference Peptide Matrix,' + header + '\n')
    self.matrix_to_string(self.peptide_diff_matrix, out_file)
    out_file.write('\n')

    out_file.write('Distinct Labeled Proteins,' + header + '\n')
    out_file.write(',' + self.dict_to_string(self.distinct_proteins_dict) + '\n' + '\n')

    out_file.write('Distinct Labeled Amino Acids,' + header + '\n')
    out_file.write(',' + self.dict_to_string(self.distinct_identifiers_dict) + '\n' + '\n')

    out_file.write('Distinct Labeled Peptides,' + header + '\n')
    out_file.write(',' + self.dict_to_string(self.distinct_peptides_dict) + '\n' + '\n')

    out_file.write('Total Labeled Proteins,Total Labeled Amino Acids,Total Labeled Peptides' + '\n')
    out_file.write(str(len(self.unique_proteins_list)) + ',' + str(len(self.unique_identifiers_list)) + ',' + str(len(self.unique_peptides_list)) + '\n' + '\n')

    out_file.write('Protein Counts,' + replicates + '\n')
    out_file.write(',' + self.dict_to_string(self.protein_counts_dict) + '\n')

    out_file.write('Protein Avgs,' + experiments + '\n')
    out_file.write(',' + self.dict_to_string(self.protein_avgs) + '\n')

    out_file.write('Protein Stdevs,' + experiments + '\n')
    out_file.write(',' + self.dict_to_string(self.protein_stdevs) + '\n' + '\n')

    out_file.write('Identifier Counts,' + replicates + '\n')
    out_file.write(',' + self.dict_to_string(self.identifier_counts_dict) + '\n')

    out_file.write('Identifier Avgs,' + experiments + '\n')
    out_file.write(',' + self.dict_to_string(self.identifier_avgs) + '\n')

    out_file.write('Identifier Stdevs,' + experiments + '\n')
    out_file.write(',' + self.dict_to_string(self.identifier_stdevs) + '\n' + '\n')

    out_file.write('Peptide Counts,' + replicates + '\n')
    out_file.write(',' + self.dict_to_string(self.peptide_counts_dict) + '\n')

    out_file.write('Peptide Avgs,' + experiments + '\n')
    out_file.write(',' + self.dict_to_string(self.peptide_avgs) + '\n')

    out_file.write('Peptide Stdevs,' + experiments + '\n')
    out_file.write(',' + self.dict_to_string(self.peptide_stdevs) + '\n')

    out_file.close()

  def matrix_to_string(self, matrix, out_file):
    for i in range(len(matrix)):
      current_row = ','.join([str(elem) for elem in matrix[i]])
      out_file.write(self.experiment_folders[i] + ',' + current_row + '\n')

  def dict_to_string(self, dct):
    st = ''
    for k in dct:
      st += str(dct[k]) + ','
    return st[:-1]

  def list_to_string(self, lst, symbol):
    return (symbol.join([str(elem) for elem in lst]))

  def write_file(self, file, header, misc, special):
    out_file = open(file, 'w')
    out_file.write(header + '\n')

    if isinstance(misc, dict):
      self.write_dict_file(out_file, misc, special)

    elif isinstance(misc, list):
      self.write_list_file(out_file, misc, special)

    else:
      print("Cannot write file, unrecognized format. " + str(misc))

    out_file.close()

  def write_dict_file(self, out_file, misc):
    for k in misc:
      out_file.write(k + ',' + str(misc[k]) + '\n')

  def write_list_file(self, out_file, misc, special):
    for i in range(len(misc)):
      current = misc[i]
      if special == True:
        st = current.output()
        out_file.write(st + '\n')
      else:
        st = ''
        for j in range(len(current)):
          st +=  str(current[j]) + ','
        out_file.write(st[:-1] + '\n')

class Experiment(Game):
  def __init__ (self, args, uniprot_db, experiment_id, cd, pd, rd):
    self.args = args
    self.uniprot_db = uniprot_db
    self.exp_dir = cd
    self.game_dir = pd
    self.results_dir = rd
    self.experiment_id = experiment_id
    self.experiment_db = args.database
    self.replicate_folders, self.replicates = [], []
    self.experiment_proteins, self.experiment_identifiers, self.experiment_peptides = [], [], []
    self.protein_counts, self.identifier_counts, self.peptide_counts = {}, {}, {}
    self.exp_protein_avg, self.exp_protein_stdev, self.exp_identifier_avg, self.exp_identifier_stdev, self.exp_peptide_avg, self.exp_peptide_stdev = 0, 0, 0, 0, 0, 0
    self.replicate_counts = 0
    self.setup_exp()
    self.replicate_no = len(self.replicate_folders)
    self.run_exp()

  def setup_exp(self):
    # ToDo: mixed database case

    if self.experiment_db == 'ip2':
      self.replicate_folders = [f for f in os.listdir() if os.path.isfile(f)]

    # single experiment
    elif (len([f for f in listdir(self.exp_dir) if isdir(join(self.exp_dir, f))]) == 0):

      # manual single replicate
      if (self.args.single_replicate == "True"):
        self.replicate_folders = [self.experiment_id]
        self.args.single_replicate = "True"
      
      # programmatic single replicate
      elif ((len([f for f in listdir(os.getcwd()) if isdir(join(os.getcwd(), f))]) == 0) | (self.args.single_replicate == "True")):
        self.replicate_folders = [self.experiment_id]
        self.args.single_replicate = "True"

     # # mulitple replicates
      else:
        self.replicate_folders = [f for f in listdir(os.getcwd()) if isfile(join(os.getcwd(), f))]
    
    # mulitple experiments
    else:
      self.replicate_folders = [f for f in listdir(self.exp_dir) if isdir(join(self.exp_dir, f))]
    
    if (".DS_Store") in self.replicate_folders:
      self.replicate_folders.remove(".DS_Store")

    if ((self.args.ptm == "True") & (self.args.experiment != "isotop")):
      self.args.pepref = "psm.tsv"
    
    if (self.args.experiment == "isotop") | (self.args.experiment == "silac") | (self.args.experiment == "silac_iodination") | (self.args.experiment == "lfq"):
      if self.args.database_version == '17':
        self.args.pepref = "peptide_label_quant.tsv"
      else:
        self.args.pepref = "modified_peptide_label_quant.tsv"

  def run_exp(self):
    for i in range(len(self.replicate_folders)):
      try:
        replicate = self.replicate_folders[i].split('-')
        replicate_no = replicate[-1]
      except:
        replicate_no = i+1

      if (self.experiment_db == "msf") & (self.args.single_replicate == "False"):
        os.chdir(self.exp_dir)
        os.chdir(self.replicate_folders[i])
      
      replicate = Replicate(self.args, self.uniprot_db, self.experiment_id, self.experiment_db, replicate_no, self.replicate_folders[i], self.game_dir, self.results_dir, self.replicate_folders)
      self.replicates.append(replicate)
      self.protein_counts[self.replicate_folders[i]] = replicate.protein_count
      self.identifier_counts[self.replicate_folders[i]] = replicate.identifier_count
      self.peptide_counts[self.replicate_folders[i]] = replicate.peptide_count      
      self.experiment_proteins = list(set(self.experiment_proteins + replicate.proteins_list))
      self.experiment_identifiers = list(set(self.experiment_identifiers + replicate.identifiers_list))      
      self.experiment_peptides = list(set(self.experiment_peptides + replicate.peptides_list))
      os.chdir(self.exp_dir)
    
    if len(self.protein_counts.keys()) > 1:
      self.exp_protein_avg, self.exp_protein_stdev = statistics.mean(self.protein_counts.values()), statistics.stdev(self.protein_counts.values())
      self.exp_identifier_avg, self.exp_identifier_stdev = statistics.mean(self.identifier_counts.values()), statistics.stdev(self.identifier_counts.values())
      self.exp_peptide_avg, self.exp_peptide_stdev = statistics.mean(self.peptide_counts.values()), statistics.stdev(self.peptide_counts.values())
    else:
      self.exp_protein_avg = self.protein_counts[self.replicate_folders[i]]
      self.exp_identifier_avg = self.identifier_counts[self.replicate_folders[i]]
      self.exp_peptide_avg = self.peptide_counts[self.replicate_folders[i]]

class Replicate(Experiment):
  def __init__ (self, args, uniprot_db, experiment_id, experiment_db, replicate_id, cd, pd, rd, replicate_folders):
    self.args = args
    self.uniprot_db = uniprot_db
    self.rep_dir = cd
    self.exp_dir = pd
    self.results_dir = rd
    self.experiment_id = experiment_id
    self.experiment_db = experiment_db
    self.replicate_no = str(replicate_id)
    self.replicate_folders = replicate_folders
    self.uniprot_dict = {}
    self.ms_list = []
    self.proteins_list = []
    self.identifiers_list = []
    self.peptides_list = []
    self.descriptions_dict = {}
    self.protein_count, self.identifier_count, self.peptide_count = 0, 0, 0
    self.header = ''
    self.run_rep()

  def run_rep(self):
    extra_folders = [f for f in listdir(os.getcwd()) if isdir(join(os.getcwd(), f))]
    if (len(extra_folders) > 0):
      os.chdir(extra_folders[0])

    if self.experiment_db == "msf":
      self.read_reference_fasta_file(self.args.proref)

      if ((self.args.ptm == "True") & (self.args.experiment != "isotop")) | (self.args.pepref == "psm.tsv"):
        self.read_peptide_file(self.args.pepref, False)
      elif (self.args.experiment == "lfq"):
        self.read_msf_quantitative_file(self.args.pepref, False, self.args.raw_ratio)
      elif (self.args.experiment == "isotop"):
        self.read_msf_quantitative_file(self.args.pepref, False, self.args.raw_ratio)
      elif (self.args.experiment == "silac") | (self.args.experiment == "silac_iodination"):
        self.read_msf_quantitative_file(self.args.pepref, True, self.args.raw_ratio)
      else:
        self.read_peptide_file(self.args.pepref, True)
    elif self.experiment_db == "ip2":
      self.read_ip2_quantitative_file(self.rep_dir)
    elif self.experiment_db == "skyline":
      self.read_skyline_quantitative_file(self.rep_dir)
    else:
      print("Which database?")

    self.protein_count = len(self.proteins_list)
    self.identifier_count = len(self.identifiers_list)
    self.peptide_count = len(self.peptides_list)

    print(self.experiment_id + " replicate-" + self.replicate_no + " Protein count: " + str(self.protein_count) + " Identifier count: " + str(self.identifier_count) + " Peptide count: " + str(self.peptide_count))
    
    os.chdir(self.results_dir)
    os.makedirs(self.experiment_id, exist_ok=True)
    os.chdir(self.experiment_id)

    if self.args.write_outfile == 'True':
      if self.replicate_no in self.rep_dir:
        self.write_file(self.experiment_id + '_' + self.rep_dir  + '_results.csv', self.header, self.ms_list, False)
      else:
        self.write_file(self.experiment_id + '_' + self.rep_dir + '-' + self.replicate_no + '_results.csv', self.header, self.ms_list, False)

  def read_reference_fasta_file(self, file):
    self.check_file_exists(file)
    for record in SeqIO.parse(file, "fasta"):
      self.uniprot_dict[record.id] = record.seq

  def read_peptide_file(self, file, pep):
    self.check_file_exists(file)
    with open(file) as fd:
      in_file = csv.reader(fd, delimiter="\t", quotechar='"')
      
      header = next(in_file, None)

      # check localization present
      if (pep == True):
        psm = self.get_peptide_header(header)
      else:
        psm = self.get_psm_header(header)

      for row in in_file:
        row = [item.replace(',', ';') for item in row]

        if ";" in self.args.aa:
          self.get_pep_mods(header, row, self.args.aa.split(";"), pep, psm)
        else:
          self.get_pep_mods(header, row, [self.args.aa], pep, psm)

  def get_peptide_header(self, header):
    modified_header = ["Peptide"]
    for i in range(len(header)):
      if (header[i] == "Peptide"):
        modified_header.append("Unmodified Peptide")
      elif (header[i] == "Charges"):
        modified_header.append("Charge")         
      else:
        modified_header.append(header[i])

    modified_header = modified_header + ["Best Localization", "Area Ratio", "Ratio", "Number of Missed Cleavages", "Modifications", "Modification Count", "Amino Acids", "Amino Acid Count", "Modification Masses", "PTM","Identifier"]
    self.header = ','.join([str(elem) for elem in modified_header])
    return False

  def get_psm_header(self, header):
    modified_header = ["Peptide"]
    for i in range(len(header)):
      if (header[i] == "Peptide"):
        modified_header.append("Unmodified Peptide")
      elif ("PeptideProphet Probability" in header[i]):
        new_probability = header[i].replace("PeptideProphet Probability", "Probability")
        modified_header.append(new_probability)
      elif ("Best Localization" in header[i]):
        modified_header.append("Best Localization")
      else:
        modified_header.append(header[i])

    if ("Best Localization" not in modified_header):
      localization = False
      modified_header = modified_header + ["Best Localization", "Spectral Count", "Area Ratio", "Ratio", "Number of Missed Cleavages", "Modifications", "Modification Count", "Amino Acids", "Amino Acid Count", "Modification Masses", "PTM","Identifier"]
    else:
      localization = True
      modified_header = modified_header + ["Spectral Count", "Area Ratio", "Ratio", "Number of Missed Cleavages", "Modifications", "Modification Count", "Amino Acids", "Amino Acid Count", "Modification Masses", "PTM","Identifier"]
      
    self.header = ','.join([str(elem) for elem in modified_header])
    return localization

  def get_pep_mods(self, header, row, aas, pep, psm):
    modification_locations = []
    identifiers = []

    if (len(row[header.index("Assigned Modifications")]) == 0) & (self.args.probe_mass != '0'):
      return
    
    protein = row[header.index("Protein ID")]
    
    modifications = row[header.index("Assigned Modifications")].replace(" ", "")
    if ";" in modifications:
      modifications = modifications.split(";")
    else:
      modifications = [modifications]
    
    miscleavages = row[header.index("Peptide")].count("K", 1, -1) + row[header.index("Peptide")].count("R", 1, -1)

    argument_aas, argument_mods = self.get_argument_mods()

    if (psm == True):
      target_mods, total_mods, total_aas, total_masses  = self.get_target_pep_mods(self.args, row[header.index('Peptide')], argument_aas, argument_mods, modifications, "[", "]", False)
    else:
      target_mods, total_mods, total_aas, total_masses  = self.get_target_pep_mods(self.args, row[header.index('Peptide')], argument_aas, argument_mods, modifications, "(", ")", False)
 
    if (self.args.experiment == "iodination"):
      misc_info = self.get_misc_iodination_info(pep, psm, row[header.index("Intensity")])
    else:
      misc_info = self.get_misc_pep_info(pep, psm)

    # skip unmodified peptides
    if len(target_mods) == 0:
      return

    modified_peptide = self.get_pep_modified_peptide(header, row, target_mods)
    identifiers = self.get_singleplex_identifiers(header, row, "Peptide", target_mods)

    if (len(identifiers) == 0):
      return
    elif (len(identifiers) == 1):
      singleplex_identifier = identifiers[0]
      protein_identifier = protein + "_" + singleplex_identifier
      ms_info = [modified_peptide] + row + misc_info + [str(miscleavages), self.list_to_string(total_mods, ';'), str(len(total_mods)), self.list_to_string(total_aas, ';'), str(len(total_aas)), self.list_to_string(total_masses, ';'), singleplex_identifier, protein_identifier]
    elif (len(identifiers) > 1):
      multiplex_identifier = self.get_multiplex_identifiers(identifiers)
      protein_identifier = protein + "_" + multiplex_identifier
      ms_info = [modified_peptide] + row + misc_info + [str(miscleavages), self.list_to_string(total_mods, ';'), str(len(total_mods)), self.list_to_string(total_aas, ';'), str(len(total_aas)), self.list_to_string(total_masses, ';'), multiplex_identifier, protein_identifier]
    else:
      print("Unexpected identifier case: " + identifiers)
      return

    self.get_update_ms_list(ms_info, protein, protein_identifier, modified_peptide)

  def get_misc_pep_info(self, pep, psm):
    if ((pep == True) & (psm == False)):
      misc_info = ['0', '0', '0']
    elif ((pep == True) & (psm == True)):
      misc_info = ['0', '0', '0']
    elif ((pep == False) & (psm == True)):
      misc_info = ['0', '1', '1']
    elif ((pep == False) & (psm == False)):
      misc_info = ['0', '0', '0', '0']
    else:
      print("Unexpected pep and psm case.")

    return misc_info

  def get_misc_iodination_info(self, pep, psm, intensity):
    if ((pep == True) & (psm == False)):
      misc_info = ['0', str(intensity), '0']
    elif ((pep == True) & (psm == True)):
      misc_info = ['0', str(intensity), '0']
    elif ((pep == False) & (psm == True)):
      misc_info = ['0', str(intensity), '1']

    #TODO: check case
    elif ((pep == False) & (psm == False)):
      misc_info = ['0', '0', str(intensity), '0']
    else:
      print("Unexpected pep and psm case.")

    return misc_info

  def get_pep_modified_peptide(self, header, row, target_mods):
    unmodified_peptide = row[header.index("Peptide")]
    modification_locations = []
    modified_peptide = ""

    for i in range(len(target_mods)):
      current_mod, current_aa = target_mods[i].split("_")
      modification_locations.append(int(current_mod)-1)

    for j in range(len(unmodified_peptide)):
      if j in modification_locations:
        modified_peptide += unmodified_peptide[j] + "*"
      else:
        modified_peptide += unmodified_peptide[j]

    return modified_peptide

  def get_residue_aa(self, identifier, peptide, uniprot_dict, uniprot_aa):
    if (identifier in uniprot_dict.keys()) & (peptide in str(uniprot_dict[identifier])):
      correct_aa = str(uniprot_dict[identifier]).index(peptide) + (int(uniprot_aa[:-1]))
      return correct_aa

  def update_lists(self, protein, identifier, peptide):
    if protein not in self.proteins_list:
      self.proteins_list.append(protein)
    if identifier not in self.identifiers_list:
      self.identifiers_list.append(identifier)   
    if peptide not in self.peptides_list:
      self.peptides_list.append(peptide) 

  def read_ip2_quantitative_file(self, file):
    self.check_file_exists(file)

    in_file = open(file, 'r')
    header = in_file.readline().strip().replace("\"", "").split(',')
    # "Peptide,Peptide Length,Charge,Probability,Spectral Count,Protein,Protein ID,Entry Name,Gene,Protein Description,Mapped Genes,Mapped Proteins,Area Ratio,Ratio,PTM,Identifier"
    modified_header = ["Peptide,Unmodified Peptide,Peptide Length"]

    for i in range(len(header)):
      if (header[i] == "CS"):
        modified_header.append("Charge")
      elif (header[i] == "Xcorr"):
        modified_header.append("Probability")
      elif (header[i] == "Protein Description"):
        modified_header.append("Protein Descriptions")
      else:
        modified_header.append(header[i])

    modified_header = modified_header + ["Area Ratio", "Modifications", "Modification Count", "Amino Acids", "Amino Acid Count", "Modification Masses", "Best Localization", "Spectral Count", "Gene", "Protein Description", "Protein ID", "PTM", "Identifier"]
    self.header = ','.join([str(elem) for elem in modified_header])

    for line in in_file:
      line = line.strip().split("\"")
      new_line = []
      for i in range(len(line)):
        if (line[i] != '') & (line[i] != ','):
          new_line.append(line[i].replace(',', ';'))

      proteins, identifiers, peptides, descriptions = self.get_ids(new_line[header.index("Proteins")], new_line[header.index("PTM index")], new_line[header.index("Sequence")], new_line[header.index("Protein Description")])
      
      argument_aas, argument_mods = self.get_argument_mods()
      numbered_mods =  self.get_assigned_mods(new_line[header.index("Sequence")], '(', ')', True)
      target_mods, total_mods, total_aas, total_masses = self.get_target_mods(argument_aas, argument_mods, numbered_mods, "(", ")", True)

      # inf, and ratios larger than 20 cases
      #TODO: consider removing peptides with inf ratios, singletons, and peptides not found in all datasets
      if ("inf" in new_line[header.index("areaRatio")]) | ("INF" in new_line[header.index("areaRatio")]):
        ratio = 20
      elif (float(new_line[header.index("areaRatio")]) > 20):
        ratio = 20
      else:
        ratio = float(new_line[header.index("areaRatio")])

      # inverted ratios case
      if (self.args.light_heavy_ratio == "heavy_light"):
        if ratio > 0:
          ratio = 1/ratio
        else:
          ratio = ratio

      # peptide, peptide length, charge, probability, spectral count, protein, protein id, entry name, gene, description, mapped genes, mapped proteins, ratio, area ratio, ptm, identifier
      for j in range(len(identifiers)):
        if "contaminant" in identifiers[j]:
          continue
        elif "Reverse" in identifiers[j]:
          continue
        elif identifiers[j].split('_')[0] in descriptions.keys():
          #self.ms_list.append([peptides[0], peptides[0].replace("*", ""), str(len(peptides[0]))] + new_line + [ratio, '', descriptions[identifiers[j].split('_')[0]], identifiers[j].split('_')[0], identifiers[j].split('_')[1], identifiers[j]])

          self.ms_list.append([peptides[0], peptides[0].replace("*", ""), str(len(peptides[0]))] + new_line + [ratio, self.list_to_string(total_mods, ';'), str(len(total_mods)), self.list_to_string(total_aas, ';'), str(len(total_aas)), self.list_to_string(total_masses, ';')] + ['0', '0', ''] + [descriptions[identifiers[j].split('_')[0]], identifiers[j].split('_')[0], identifiers[j].split('_')[1], identifiers[j]])
        else:
          #self.ms_list.append([peptides[0], peptides[0].replace("*", ""), str(len(peptides[0]))] + new_line + [ratio, '', '', identifiers[j].split('_')[0], identifiers[j].split('_')[1], identifiers[j]])

          self.ms_list.append([peptides[0], peptides[0].replace("*", ""), str(len(peptides[0]))] + new_line + [ratio, self.list_to_string(total_mods, ';'), str(len(total_mods)), self.list_to_string(total_aas, ';'), str(len(total_aas)), self.list_to_string(total_masses, ';')] + ['0', '0', '', ''] + [identifiers[j].split('_')[0], identifiers[j].split('_')[1], identifiers[j]])

      self.proteins_list = list(set(self.proteins_list + proteins))
      self.identifiers_list = list(set(self.identifiers_list + identifiers))
      self.peptides_list = list(set(self.peptides_list + peptides))
      self.descriptions_dict.update(descriptions)

    in_file.close() 

  def read_msf_quantitative_file(self, file, silac, raw):
    self.check_file_exists(file)
    with open(file) as fd:
      in_file = csv.reader(fd, delimiter="\t", quotechar='"')
      header = next(in_file, None)
      
      if (self.args.experiment == "silac") | (self.args.experiment == "silac_iodination"):
        self.args.light_probe_mass = self.args.heavy_probe_mass
        self.get_silac_header(header)
      else:
        self.get_msf_header(header)

      for row in in_file:
        row = [item.replace(',', ';') for item in row]

        # skip rows without ratios for silac
        if (len(row[header.index("Log2 Ratio HL")]) == 0) & (len(row[header.index("Heavy Modified Peptide")]) == 0) & (silac == True):
          continue

        protein = row[header.index("Protein ID")]
        modified_peptide, light_peptide, heavy_peptide, flagged_light, flagged_heavy, ratio = self.get_reformatted_mods(header, row, silac, raw)
        miscleavages = row[header.index("Peptide Sequence")].count("K", 1, -1) + row[header.index("Peptide Sequence")].count("R", 1, -1)
        
        # ['463.2366', '467.2529']
        argument_aas, argument_mods = self.get_argument_mods()

        # TODO: Check for non-mouse data multiplexed peptides do not have +10 on +1 labeled sites
        #['20_C[463.2366]', '20_C[467.2529]', '2_C[467.2529]', '2_C[463.2366]']
        numbered_mods = self.get_numbered_mods(light_peptide, heavy_peptide)

        # 25_C[463.2366]
        target_mods, total_mods, total_aas, total_masses = self.get_target_mods(argument_aas, argument_mods, numbered_mods, "[", "]", False)
        reformatted_probability = self.get_reformatted_probability(row[header.index("Best Light PeptideProphet Probability")], row[header.index("Best Heavy PeptideProphet Probability")])

        # skip unmodified peptides
        if (len(target_mods) == 0):
          continue

        identifiers = self.get_singleplex_identifiers(header, row, "Peptide Sequence", target_mods)

        # skip peptides that are both modified but do not have a ratio
        check_ratio = isinstance(ratio, str)
        if (check_ratio == True):
          if (len(ratio) == 0):
            continue

        # skip peptides that do not have matching uniprot identifiers/sequences
        if (len(identifiers) == 0):
          continue
        elif (len(identifiers) == 1):
          singleplex_identifier = identifiers[0]
          protein_identifier = protein + "_" + singleplex_identifier
          ms_info = [modified_peptide, flagged_light, flagged_heavy] + row + [self.list_to_string(numbered_mods, ';'), str(ratio), '0', reformatted_probability, str(miscleavages), self.list_to_string(total_mods, ';'), str(len(total_mods)), self.list_to_string(total_aas, ';'), str(len(total_aas)), self.list_to_string(total_masses, ';')] + [singleplex_identifier, protein_identifier]
        elif (len(identifiers) > 1):
          multiplex_identifier = self.get_multiplex_identifiers(identifiers)
          protein_identifier = protein + "_" + multiplex_identifier
          ms_info = [modified_peptide, flagged_light, flagged_heavy] + row + [self.list_to_string(numbered_mods, ';'), str(ratio), '0', reformatted_probability, str(miscleavages), self.list_to_string(total_mods, ';'), str(len(total_mods)), self.list_to_string(total_aas, ';'), str(len(total_aas)), self.list_to_string(total_masses, ';')] + [multiplex_identifier, protein_identifier]          
        else:
          print("Unexpected identifier case: " + identifiers)
          continue

        self.get_update_ms_list(ms_info, protein, protein_identifier, modified_peptide)

  def get_msf_header(self, header):
    modified_header = ["Peptide", "Light Peptide", "Heavy Peptide"]
    for i in range(len(header)):
      if (header[i] == "Peptide Sequence"):
        modified_header.append("Unmodified Peptide")
      elif (header[i] == "Charges"):
        modified_header.append("Charge")
      elif (header[i] == "Label Count"):
        modified_header.append("Spectral Count")
      else:
        modified_header.append(header[i])
    modified_header = modified_header + ["Assigned Modifications", "Area Ratio", "Best Localization", "Probability", "Number of Missed Cleavages", "Modifications", "Modification Count", "Amino Acids", "Amino Acid Count", "Modification Masses", "PTM","Identifier"]
    self.header = ','.join([str(elem) for elem in modified_header])

  def get_silac_header(self, header):
    # modified_header = ["Peptide", "Light Peptide", "Heavy Peptide"]
    modified_header = ["Unmodified Peptide", "Light Peptide", "Heavy Peptide"]
    for i in range(len(header)):
      if (header[i] == "Peptide Sequence"):
        modified_header.append("Peptide")
      elif (header[i] == "Charges"):
        modified_header.append("Charge")
      elif (header[i] == "Label Count"):
        modified_header.append("Spectral Count")
      else:
        modified_header.append(header[i])
    modified_header = modified_header + ["Assigned Modifications", "Area Ratio", "Best Localization", "Probability", "Number of Missed Cleavages", "Modifications", "Modification Count", "Amino Acids", "Amino Acid Count", "Modification Masses", "PTM","Identifier"]
    self.header = ','.join([str(elem) for elem in modified_header])

  # TODO: Note to self, mouse multiplexed have +10 on modificiation for second labeled site on a single peptide for ratio_quant
  # Check that this numbering is consistent for human data one day
  def get_assigned_mods(self, peptide, symbol1, symbol2, numbering):
    idxs = self.get_paired_idxs(peptide, symbol1, symbol2)

    inside_data = []

    for i in range(0, len(idxs), 2):
      if numbering == True:
        if (i >=2 ):
          # st = str(idxs[i]) + '_' + ((peptide[idxs[i]-1:idxs[i+1]+1]) - 10)
          st = str(idxs[i]-10) + '_' + peptide[idxs[i]-1:idxs[i+1]+1]

        else:
          st = str(idxs[i]) + '_' + peptide[idxs[i]-1:idxs[i+1]+1]
      else:
        st = peptide[idxs[i]-1:idxs[i+1]+1]
      
      if (st not in inside_data):
        inside_data.append(st)
  
    return inside_data

  def get_assigned_mods_V0(self, peptide, symbol1, symbol2, numbering):
    idxs = self.get_paired_idxs(peptide, symbol1, symbol2)

    inside_data = []

    for i in range(0, len(idxs), 2):
      if numbering == True:
        if (i >=2 ):
          st = str(idxs[i]-10) + '_' + peptide[idxs[i]-1:idxs[i+1]+1]

        else:
          st = str(idxs[i]) + '_' + peptide[idxs[i]-1:idxs[i+1]+1]
      else:
        st = peptide[idxs[i]-1:idxs[i+1]+1]
      
      if (st not in inside_data):
        inside_data.append(st)
  
    return inside_data

  def get_reformatted_mods(self, header, row, silac, raw):

    # if silac != True:
    #   light_mods, heavy_mods, light_pep, heavy_pep = self.get_original_mods(header, row)
    # else:
    light_mods, heavy_mods, light_pep, heavy_pep = self.get_original_mods(header, row)      

    found_light, found_heavy, flagged_light, flagged_heavy = self.get_quant_modified_peptides(header, row, light_mods, heavy_mods, light_pep, heavy_pep)
  
    if silac != True:
      peptide, ratio = self.get_modified_isotop_ratio(header, row, found_light, found_heavy, flagged_light, flagged_heavy, raw)
    else:
      peptide, ratio = self.get_modified_silac_ratio(header, row, found_light, found_heavy, flagged_light, flagged_heavy, raw)

    return peptide, light_pep, heavy_pep, flagged_light, flagged_heavy, ratio

  def get_original_mods_V0(self, header, row):
    extra_mods = self.get_assigned_mods(row[header.index("Modified Peptide")], '[', ']', False)
    light_mods = self.get_assigned_mods(row[header.index("Light Modified Peptide")],'[', ']', False)
    heavy_mods = self.get_assigned_mods(row[header.index("Heavy Modified Peptide")], '[', ']',False)

    light_pep = flagged_light = row[header.index("Light Modified Peptide")]
    heavy_pep = flagged_heavy = row[header.index("Heavy Modified Peptide")]

    for i in range(len(extra_mods)):
      current_mod = extra_mods[i]
      if current_mod in light_pep:
        if "n" in current_mod:
          light_pep = light_pep.replace(current_mod[current_mod.index("[")-1:], "")
        else:
          light_pep = light_pep.replace(current_mod[current_mod.index("["):], "")
        light_mods.remove(current_mod)     
      if current_mod in heavy_pep:
        if "n" in heavy_pep:
          heavy_pep = heavy_pep.replace(current_mod[current_mod.index("[")-1:], "")
        else:
          heavy_pep = heavy_pep.replace(current_mod[current_mod.index("["):], "")
          heavy_mods.remove(current_mod)

    return light_mods, heavy_mods, light_pep, heavy_pep

  def get_original_mods(self, header, row):

    argument_aas, argument_mods = self.get_argument_mods()

    # extra_mods = self.get_assigned_mods(row[header.index("Modified Peptide")], '[', ']', False)
    light_mods = self.get_assigned_mods(row[header.index("Light Modified Peptide")],'[', ']', False)
    heavy_mods = self.get_assigned_mods(row[header.index("Heavy Modified Peptide")], '[', ']',False)

    light_pep = flagged_light = row[header.index("Light Modified Peptide")]
    heavy_pep = flagged_heavy = row[header.index("Heavy Modified Peptide")]

    tot_mods = list(set(light_mods + heavy_mods))

    for i in range(len(tot_mods)):
      current_mod = tot_mods[i]
      current_mass =  current_mod.split('[', 1)[1].split(']')[0]
      if current_mass not in argument_mods:

        light_pep, light_mods = self.check_mod(current_mod, light_pep, light_mods)
        heavy_pep, heavy_mods = self.check_mod(current_mod, heavy_pep, heavy_mods)

    return light_mods, heavy_mods, light_pep, heavy_pep

  def check_mod(self, current_mod, pep, mods):
    if current_mod in pep:
      if "n" in pep:
        pep = pep.replace(current_mod[current_mod.index("[")-1:], "")
      else:
        pep = pep.replace(current_mod[current_mod.index("["):], "")
      mods.remove(current_mod)

    return pep, mods

  def get_quant_modified_peptides(self, header, row, light_mods, heavy_mods, light_pep, heavy_pep):
    found_light = found_heavy = True 
    if len(light_mods) == 0:
      found_light = False 
      flagged_light = light_pep
    else:
      flagged_light = self.get_flagged_sequence(light_pep)

    if len(heavy_mods) == 0:
      found_heavy = False 
      flagged_heavy = heavy_pep
    else:
      flagged_heavy = self.get_flagged_sequence(heavy_pep)

    return found_light, found_heavy, flagged_light, flagged_heavy

  def get_modified_isotop_ratio(self, header, row, found_light, found_heavy, flagged_light, flagged_heavy, raw):

    if (found_heavy == True) & (found_light == True):
      # skip peptides that are both modified but do not have a ratio
      if (row[header.index("Log2 Ratio HL")] == ''):
        return "", ""
      if (self.args.light_heavy_ratio == "heavy_light"):
        old_ratio = 2 ** float(row[header.index("Log2 Ratio HL")])
        raw_ratio = 1/old_ratio
        if raw == 'True':
          if raw_ratio > 20:
            ratio = 20
          else:
            ratio = raw_ratio
        else:
          ratio = str(math.log(raw_ratio,2))
      else:
        if raw == 'True':
          ratio = 2 ** float(row[header.index("Log2 Ratio HL")])
          if ratio > 20:
            ratio = 20
        else:
          ratio = row[header.index("Log2 Ratio HL")]
      peptide = flagged_light
    elif (found_heavy == True) & (found_light == False):
      if raw == 'True':
        ratio = 20
      else:
        ratio = math.log(20, 2)
      peptide = flagged_heavy
    elif (found_heavy == False) & (found_light == True):
      if raw == 'True':
        ratio = 1/20
      else:
        ratio = math.log((1/20), 2)
      peptide = flagged_light
    else:
      return "", ""

    return peptide, ratio

  def get_modified_silac_ratio(self, header, row, found_light, found_heavy, flagged_light, flagged_heavy, raw):

    if (found_heavy == True) & (len(row[header.index("Log2 Ratio HL")]) == 0):
      if raw == 'True':
        ratio = 20
      else:
        ratio = math.log(20,2)
      peptide = flagged_heavy
    elif (found_heavy == False) & (len(row[header.index("Log2 Ratio HL")]) == 0):
      ratio = 0
      peptide = flagged_light
    elif (found_heavy == True) & (len(row[header.index("Log2 Ratio HL")]) != 0):
      if (self.args.light_heavy_ratio == "heavy_light"):
        old_ratio = 2 ** float(row[header.index("Log2 Ratio HL")])
        raw_ratio = 1/old_ratio
        if raw == 'True':
          ratio = raw_ratio
        else:
          ratio = str(math.log(raw_ratio,2))
      else:
        if raw == 'True':
          ratio = 2 ** float(row[header.index("Log2 Ratio HL")])
        else:
          ratio = row[header.index("Log2 Ratio HL")]
      peptide = flagged_heavy
    else:
      # print("Unidentified ratio: " + str(row[header.index("Log2 Ratio HL")]))
      return None, None

    return peptide, ratio

  def get_argument_mods(self):
    if (";" in self.args.aa):
      args_aa = self.args.aa.split(";")
    else:
      args_aa = [self.args.aa]

    if (";" in self.args.probe_mass):
      args_pm = self.args.probe_mass.split(";")
    else:
      args_pm = [self.args.probe_mass]

    if (";" in self.args.light_probe_mass):
      args_lpm = self.args.light_probe_mass.split(";")
    else:
      args_lpm = [self.args.light_probe_mass]

    if (";" in self.args.heavy_probe_mass):
      args_hpm = self.args.heavy_probe_mass.split(";")
    else:
      args_hpm = [self.args.heavy_probe_mass]
    
    if (self.args.experiment == "identification") | (self.args.experiment == "iodination"):
      args_mod = args_pm
    else:
      args_mod = list(set(args_lpm + args_hpm))

    return args_aa, args_mod

  def get_numbered_mods(self, light_peptide, heavy_peptide):
    numbered_light_mods = self.get_assigned_mods(light_peptide, '[', ']', True)
    numbered_heavy_mods = self.get_assigned_mods(heavy_peptide, '[', ']', True)

    total_mods = list(set(numbered_light_mods + numbered_heavy_mods))

    return total_mods

  def get_target_pep_mods(self, args, peptide, args_aa, args_mod, modifications, symbol1, symbol2, ip2):
    aas, mods = self.get_argument_mods()

    target_mods = []
    total_mods = []
    total_aas = []
    total_masses = []

    if self.args.probe_mass == '0':
      num = '1'
      aa = peptide[0]
      new_mod = num + "_" + aa
      target_mods.append(new_mod) 
      if (aa not in total_aas):
        total_aas.append(aa)
      return target_mods, total_mods, total_aas, total_masses

    #4C(463.2366)
    for h in range(len(modifications)): 
      outside, mod = self.get_inside_data(modifications[h], symbol1, symbol2)

      # N-term case
      if (args.aa == 'N-term'):
        num = '1'
        aa = peptide[0]
      else:
        # num = outside.replace(aa, '')
        num = outside[:-1]
        aa = ''.join([i for i in outside if i.isalpha()])

      if args.aa == 'N-term':
        if (mod in args_mod):
          if (aa not in target_mods):
            new_mod = num + "_" + aa
            target_mods.append(new_mod) 
          if (aa not in total_aas):
            total_aas.append(aa)
          if (mod not in total_masses):
            total_masses.append(mod)
          if (modifications[h] not in total_mods):
            total_mods.append(modifications[h])
      else:
        if (mod in args_mod) & (aa in args_aa):
          if (aa not in target_mods):
            new_mod = num + "_" + aa
            target_mods.append(new_mod) 
          if (aa not in total_aas):
            total_aas.append(aa)
          if (mod not in total_masses):
            total_masses.append(mod)
          if (modifications[h] not in total_mods):
            total_mods.append(modifications[h])

    return target_mods, total_mods, total_aas, total_masses

  def get_target_mods(self, args_aa, args_mod, pre_modifications, symbol1, symbol2, ip2):
    aas, mods = self.get_argument_mods()

    target_mods = []
    total_mods = []
    total_aas = []
    total_masses = []

    modifications = pre_modifications

    for h in range(len(pre_modifications)):

      if ("[" in symbol1):
        pattern = r'\[[^()]*\]'
      elif ("(" in symbol1):
        pattern = r'\((.*?)\)'
      else:
        print("Symbol for pattern not identified.")
        sys.exit()

      separate = re.sub(pattern, '', pre_modifications[h])

      if "_" in separate:
        aa = separate.split("_")[1]
        num = separate.split("_")[0]
      elif (ip2 == True):
        aa = separate[0]
        num = ''
      else:
        aa = separate[-1]
        num = separate[:-1]
      
      try:
        mod = re.search(pattern, modifications[h]).group(0)[1:-1]
      except AttributeError:
        if ("[" in symbol1):
          mod = re.search(r'\((.*?)\)', modifications[h]).group(0)[1:-1]
        elif ("(" in symbol1):
          mod = re.search(r'\((.*?)\)', modifications[h]).group(0)[1:-1]
      except:
        print("Unexpected error in modification: " + modifications[h])
        sys.exit()

      if (mod in args_mod) & (aa in args_aa):
        if (aa not in target_mods):
          new_mod = num + "_" + aa
          target_mods.append(new_mod) 
        if (aa not in total_aas):
          total_aas.append(aa)
        if (mod not in total_masses):
          total_masses.append(mod)
        if (modifications[h] not in total_mods):
          total_mods.append(modifications[h])

    return target_mods, total_mods, total_aas, total_masses
 
  # reformat embedded lists separated by commas
  def get_inside_data(self, st, symbol1, symbol2):
    # get information inside ()
    idxs = self.get_idxs(st, symbol1, symbol2)

    inside_data = ""
    for i in range(0, len(idxs), 2):
      inside_data += st[idxs[i]:idxs[i+1]+1]

    # get information outside ()
    outside_data = st[:idxs[0]]
    return outside_data, str(inside_data[1:-1])

  # get indeces of pairs of [] or ()
  def get_idxs(self, st, ch1, ch2):
    idxs = [i for i, a in enumerate(st) if (a == ch1) or (a == ch2)]
    return idxs

  def get_reformatted_probability(self, light, heavy):

    new_probability = ""
    if len(light) == 0:
      new_probability += "L:0"
    else:
      new_probability += "L:" + str(light)

    if len(heavy) == 0:
      new_probability += ";H:0"
    else:
      new_probability += ";H:" + str(heavy)
    
    return new_probability

  def get_singleplex_identifiers(self, header, row, peptide_column, target_mods):
    single_identifiers = []

    for i in range(len(target_mods)):
      current_num, current_aa = target_mods[i].split("_")
      # Check 
      identifier = self.get_residue_aa(row[header.index("Protein")].strip(), row[header.index(peptide_column)].strip(), self.uniprot_dict, current_num + current_aa)
      if identifier != None:
        full_identifier = current_aa + str(identifier)
        single_identifiers.append(full_identifier)
    return single_identifiers

  def get_multiplex_identifiers(self, single_identifiers):
    sorted_identifiers = sorted(single_identifiers)
    
    full_identifier = ""
    for k in range(len(sorted_identifiers)):
      if sorted_identifiers[k] not in full_identifier:
        full_identifier += sorted_identifiers[k] + '_'
    
    full_identifier = full_identifier[:-1]
    return full_identifier
          
  def get_update_ms_list(self, ms_info, protein, protein_identifier, modified_peptide):
    self.ms_list.append(ms_info)
    self.update_lists(protein, protein_identifier, modified_peptide)

  def read_skyline_quantitative_file(self, file):
    self.check_file_exists(file)
    in_file = open(file, 'r')
    header_title = in_file.readline().strip()
    header = header_title.split(',')
    modified_header = ["Peptide,Unmodified Peptide,Peptide Length,Probability,Charge"] + header +  ["Area Ratio", "Best Localization", "Protein ID", "Protein Description","PTM","Identifier"]
    self.header = ','.join([str(elem) for elem in modified_header])

    for line in in_file:
      line = line.strip().split(",")

      # empty ratio case
      if ("NAN" in line[header.index("RatioLightToHeavy")]) | ("#N/A" in line[header.index("RatioLightToHeavy")]):
        continue

      # contaminated case
      if "contaminant" in line[header.index("Protein Name")]:
        continue

      if "Library Peptides" in line[header.index("Protein Name")]:
        continue

      protein = line[header.index("Protein Name")].split("|")[1]
      gene = line[header.index("Protein Gene")]
      peptide = line[header.index("Peptide Modified Sequence Three Letter Codes")]

      # oxidated peptides case
      if ("Oxi" in line[header.index("Peptide Modified Sequence Three Letter Codes")]) | ("1Ac" in line[header.index("Peptide Modified Sequence Three Letter Codes")]):
        continue

      flagged_sequence = self.get_flagged_sequence(peptide)

      # empty tag case
      if "*" not in flagged_sequence:
        continue

      # ratios larger than 20 cases
      #TODO: consider removing peptides with inf ratios, singletons, and peptides not found in all datasets
      
      if (line[header.index("RatioLightToHeavy")] == 'NaN') | (line[header.index("RatioLightToHeavy")] == '#N/A'):
        ratio = 0
      elif (float(line[header.index("RatioLightToHeavy")]) > 20):
        ratio = 20
      else:
        ratio = float(line[header.index("RatioLightToHeavy")])

      # inverted ratios case
      if (self.args.light_heavy_ratio == "heavy_light"):
        if ratio > 0:
          ratio = 1/ratio
        else:
          ratio = ratio

      # get identifier
      if protein in self.uniprot_db.uniprot_names:
        uniprot_idx = self.uniprot_db.uniprot_names.index(protein)
        uniprot_protein = self.uniprot_db.uniprot_proteins[uniprot_idx]
        
        identifiers = self.get_skyline_identifiers(protein, flagged_sequence, uniprot_protein)
      else:
        print("Protein " + protein + " not in Uniprot Reference File.")
        continue

      #ToDo Update
      for j in range(len(identifiers)):
        data = [flagged_sequence, flagged_sequence.replace("*",""), len(flagged_sequence)-1, '', ''] + line + [str(ratio), '', str(identifiers[j].split('_')[0]), '', str(identifiers[j].split('_')[1]), identifiers[j]]
        if data not in self.ms_list:
          self.ms_list.append(data)

      self.proteins_list = list(set(self.proteins_list + [protein]))
      self.identifiers_list = list(set(self.identifiers_list + identifiers))
      self.peptides_list = list(set(self.peptides_list + [flagged_sequence]))

    in_file.close() 

  def get_ids(self, protein, identifier, sequence, description):
    identifiers = []
    proteins = []
    peptides = []
    descriptions = {}
    
    idxs = []
    for i in range(len(sequence)):
      if sequence[i] == ".":
        idxs.append(i)

    flagged_sequence = re.sub(r'\([^)]*\)', '*', sequence[idxs[0]+1:idxs[-1]])
    if flagged_sequence not in peptides:
      peptides.append(flagged_sequence)

    if ";" in protein:
      # multiple protein case
      pros = protein.replace('[', '').replace(']', '').split(';')
      
      # multiple identifiers for multiple proteins
      sep_identifiers = identifier.split(' ; ')

      for i in range(len(pros)):
        ptms = self.get_new_identifier(sep_identifiers[i].split(';'))
        if pros[i].strip() + '_' + ptms.strip() not in identifiers:
          identifiers.append(pros[i].strip() + '_' + ptms.strip())
        if pros[i].strip() not in proteins:
          proteins.append(pros[i].strip())
    else:
      # single protein case
      protein = protein.replace('[', '').replace(']', '')
      if protein.strip() not in proteins:
        proteins.append(protein.strip())

      # mulitple identifiers for single protein
      if ";" in identifier:
        ptms = self.get_new_identifier(identifier.split(';'))
        if protein.strip() + '_' + ptms.strip() not in identifiers:
          identifiers.append(protein.strip() + '_' + ptms.strip())
      # single identifier for single protein
      else:
        if protein.strip() + '_' + identifier[1:].strip() not in identifiers:
          identifiers.append(protein.strip() + '_' + self.args.aa + identifier[1:].strip())
    
    if ";" in protein:
      desc = description.replace('[', '').replace(']', '').replace(',', ' ').split(';')
      for k in range(len(desc)):
        if desc[k].strip()[:6] not in descriptions.keys():
          descriptions[desc[k].strip()[:6]] = desc[k].strip()
    else:
      desc = description.replace('[', '').replace(']', '').replace(';', '')
      if desc.strip()[:6] not in descriptions.keys():
        descriptions[desc.strip()[:6]] = desc.strip()

    return proteins, identifiers, peptides, descriptions

  def get_new_identifier(self, sep_identifiers):
    ptms = ''
    # sep_identifiers = identifier.split(';')
    for i in range(len(sep_identifiers)):
      ptms += sep_identifiers[i].strip() + '_'
    return ptms[:-1]

  # get indeces of pairs of [] or ()
  def get_paired_idxs(self, st, ch1, ch2):
    idxs = [i for i, a in enumerate(st) if (a == ch1) or (a == ch2)]
    return idxs

  # get sequence with tagged C's flagged by an "*"
  def get_flagged_sequence(self, st):
    flagged_sequence = st

    # get information inside ()
    idxs = self.get_paired_idxs(st, "[", "]")

    # get location of tag, replace tag with an "*"
    inside_data = []
    for i in range(0, len(idxs), 2):
      if (st[idxs[i]:idxs[i+1]+1] not in inside_data):
        inside_data.append(st[idxs[i]:idxs[i+1]+1])

    for j in range(len(inside_data)):
      flagged_sequence = flagged_sequence.replace(inside_data[j], "*")
    
    return flagged_sequence

  def get_skyline_identifiers(self, protein, flagged_sequence, uniprot_protein):
    flag_count = flagged_sequence.count("*")
    flag_idxs = []
    identifiers_list = []

    for i in range(len(flagged_sequence)):
      if flagged_sequence[i] == "*":
        if flagged_sequence[i-1] == "C":
          flag_idxs.append(i)

    if len(flag_idxs) > 1:
      for j in range(len(flag_idxs)):    
        flag_idx = flag_idxs[j]
        unflagged_sequence = flagged_sequence.replace("C*", "C")
        if unflagged_sequence in uniprot_protein.sequence:
          start_idx = uniprot_protein.sequence.index(unflagged_sequence)
          # index of cysteine + 1
          end_idx = start_idx + int(flag_idx)
          identifier = protein + "_C" + str(end_idx)
          if identifier not in identifiers_list:
            identifiers_list.append(identifier)
        else:
          print(uniprot_protein.sequence, unflagged_sequence)
    
    return identifiers_list

class Protein:
  def __init__ (self, protein, description, name, identifier, peptide, spectral_count, modifications, amino_acids, masses, ratio, dataset, pep_files, args, experiment_id, experiments):
    self.args = args
    self.protein = protein
    self.name = name
    self.description = description
    self.spectral_count = spectral_count
    self.modifications, self.modification_count, self.masses, self.amino_acids, self.amino_acid_count = [], [], [], [], []
    self.identifiers, self.peptides, self.ratios = [], [], []
    self.peptide_count = 0
    self.dataset_dict, self.ratio_dict = {}, {}
    self.ratio_mean, self.ratio_median, self.ratio_stdev, self.ratio_median_mean_list = [], [], [], []
    self.ratio_median_mean = 0
    self.experiment_dataset_dict, self.dataset_dict, self.experiment_rep_dict = {}, {}, {}
    self.experiment_reps, self.dataset_reps = 0, 0
    self.experiment_ratio_dict, self.experiment_stdev_dict = {}, {}
    self.aggregate_median_mean, self.aggregate_median_mean_stdev = 0, 0
    for h in range(len(experiments)):
      self.experiment_dataset_dict[experiments[h]] = []
      self.experiment_rep_dict[experiments[h]] = 0
    for i in range(len(pep_files)):
      self.dataset_dict[pep_files[i]] = 0
      self.ratio_dict[pep_files[i]] = []
    self.ratio_mean_dict, self.ratio_median_dict, self.ratio_stdev_dict = self.ratio_dict.copy(), self.ratio_dict.copy(), self.ratio_dict.copy()

    self.update(identifier, peptide, spectral_count, modifications, amino_acids, masses, ratio, dataset, experiment_id)

  def update(self, identifier, peptide, spectral_count, modifications, amino_acids, masses, ratio, dataset, experiment_id):
    self.identifiers = self.update_list(identifier, self.identifiers)
    self.peptides = self.update_list(peptide, self.peptides)
    self.peptide_count = len(self.peptides)

    self.spectral_count = self.update_count(spectral_count, self.spectral_count)
    self.modifications = self.update_semicolon_list(modifications, self.modifications)
    self.modification_count  = len(self.modifications)
    self.masses = self.update_semicolon_list(masses, self.masses)
    self.amino_acids = self.update_semicolon_list(amino_acids, self.amino_acids)
    self.amino_acid_count = len(self.amino_acids)

    self.update_datasets(dataset, experiment_id)
    self.update_experiments(experiment_id)

    if (self.args.experiment != "identification"):
      self.update_ratios(ratio, dataset)
      self.update_median_means()
      self.update_aggregate_median_mean()

  def update_list(self, item, old_list):
    new_list = old_list
    if (item not in new_list) & (item != ''):
      new_list.append(item)
    return new_list

  def update_semicolon_list(self, item, old_list):
    if (";" in item):
      current = item.split(";")
    else:
      current = [item]

    for i in range(len(current)):
      old_list = self.update_list(current[i], old_list)
    return old_list

  def update_count(self, item, old_count):
    return int(item) + int(old_count)

  def update_datasets(self, dataset, experiment_id):
    self.dataset_dict[dataset] += 1

    if dataset not in self.experiment_dataset_dict[experiment_id]:
      self.experiment_dataset_dict[experiment_id] = self.experiment_dataset_dict[experiment_id] + [dataset]
    
    for k in self.experiment_dataset_dict:
      self.experiment_rep_dict[k] = len(self.experiment_dataset_dict[k])

    self.dataset_reps = 0
    for l in self.dataset_dict:
      if self.dataset_dict[l] >= 1:
        self.dataset_reps += 1

  def update_experiments(self, experiment_id):
    self.experiment_reps = 0
    for m in self.experiment_rep_dict:
      if (self.experiment_rep_dict[m] > 0):
        self.experiment_reps += 1

  def update_ratios(self, ratio, dataset): 
    # update ratio dictionary
    if float(ratio.strip()) != 0:
      self.ratios.append(float(ratio.strip()))
      self.ratio_dict[dataset] = self.ratio_dict[dataset] + [float(ratio.strip())]
    
    # update ratio median dictionary
    for k in self.ratio_dict:
      if len(self.ratio_dict[k]) == 0:
        continue
      elif len(self.ratio_dict[k]) == 1:
        self.ratio_median_dict[k] = self.ratio_dict[k][0]
      else:
        self.ratio_median_dict[k] = statistics.median(self.ratio_dict[k])
  
  def update_median_means(self):
    # update ratio mean of medians dictionary
    for j in self.experiment_dataset_dict:
      current_experiment = self.experiment_dataset_dict[j]
      experiment_ratios = []

      for i in range(len(current_experiment)):
        if (current_experiment[i] in self.ratio_median_dict.keys()):
          if self.ratio_median_dict[current_experiment[i]]:
            experiment_ratios.append(self.ratio_median_dict[current_experiment[i]])

      # [] case
      if not experiment_ratios:
        self.experiment_ratio_dict[j] = 0
        self.experiment_stdev_dict[j] = 0
      # single ratio case
      elif (len(experiment_ratios) == 1):
        self.experiment_ratio_dict[j] = experiment_ratios[0]
        self.experiment_stdev_dict[j] = 0
      # multiple ratio case
      elif len(experiment_ratios) > 1:
        self.experiment_ratio_dict[j] = statistics.mean(experiment_ratios)
        self.experiment_stdev_dict[j] = statistics.stdev(experiment_ratios)
      else:
        self.experiment_ratio_dict[j] = 0
        self.experiment_stdev_dict[j] = 0

  def update_aggregate_median_mean(self):
    aggregate_ratios = []
    for k in self.experiment_ratio_dict:
      if self.experiment_ratio_dict[k]:
        aggregate_ratios.append(self.experiment_ratio_dict[k])

    if len(aggregate_ratios) == 1:
      self.aggregate_median_mean = aggregate_ratios[0]
      self.aggregate_median_mean_stdev = 0
    elif len(aggregate_ratios) > 1:
      self.aggregate_median_mean = statistics.mean(aggregate_ratios)
      self.aggregate_median_mean_stdev = statistics.stdev(aggregate_ratios)  
    else:
      self.aggregate_median_mean = 0
      self.aggregate_median_mean_stdev = 0

  def replace(self, misc, dictionary):
    if dictionary == True:
      tmp = {}
      for k in misc:
        if misc[k] == 0:
          tmp[k] = float("nan")
        else:
          tmp[k] = misc[k]
    else:
      tmp = []
      for i in range(len(misc)):
        if misc[i] == 0:
          tmp.append(float('nan'))
        else:
          tmp.append(misc[i])

    return tmp

  def list_to_string(self, lst):
    return (';'.join([str(elem) for elem in lst]))

  def dict_to_string(self, dct):
    st = ''
    for k in dct:
      st += str(dct[k]) + ','
    return st[:-1]

  def output(self):

    if (self.args.experiment == "identification"):
      return self.protein + ',' + self.description + ',' + self.name + ',' + self.list_to_string(self.identifiers)  + ',' + self.list_to_string(self.peptides) + ',' + str(self.peptide_count) + ',' + str(self.spectral_count) + ',' + self.list_to_string(self.modifications) + ',' + str(self.modification_count) + ',' + self.list_to_string(self.amino_acids) + ',' + str(self.amino_acid_count) + ',' + self.list_to_string(self.masses) + ',' + self.list_to_string(self.replace([self.experiment_reps], False)) + ',' + self.list_to_string(self.replace([self.dataset_reps], False)) + ',' + self.dict_to_string(self.replace(self.experiment_rep_dict, True)) + ',' + self.dict_to_string(self.replace(self.dataset_dict, True)) + ','
    
    else:

      experiment_ratio_dict_str = self.dict_to_string(self.replace(self.experiment_ratio_dict, True))
      if '[]' in experiment_ratio_dict_str:
        experiment_ratio_dict_str = experiment_ratio_dict_str.replace('[]', '')

      return self.protein + ',' + self.description + ',' + self.name + ',' + self.list_to_string(self.identifiers)  + ',' + self.list_to_string(self.peptides) + ',' + str(self.peptide_count) + ',' + str(self.spectral_count) + ',' + self.list_to_string(self.modifications) + ',' + str(self.modification_count) + ',' + self.list_to_string(self.amino_acids) + ',' + str(self.amino_acid_count) + ',' + self.list_to_string(self.masses) + ',' + self.list_to_string(self.replace([self.experiment_reps], False)) + ',' + self.list_to_string(self.replace([self.dataset_reps], False)) + ',' + self.dict_to_string(self.replace(self.experiment_rep_dict, True)) + ',' + self.dict_to_string(self.replace(self.dataset_dict, True)) + ','+ self.list_to_string(self.replace([self.aggregate_median_mean], False)) + ',' + self.list_to_string(self.replace([self.aggregate_median_mean_stdev], False)) + ',' + experiment_ratio_dict_str + ',' + self.dict_to_string(self.replace(self.experiment_stdev_dict, True)) + ','+ self.list_to_string(self.ratios) + ',' + self.dict_to_string(self.replace(self.ratio_median_dict, True)).replace('[]', '') + ','

class Identifier(Protein):
  def __init__ (self, protein, description, identifier, name, peptide, spectral_count, modifications, amino_acids, masses, ratio, dataset, pep_files, args, experiment_id, experiments):
    self.args = args
    self.identifier = identifier
    self.protein = protein
    self.description = description
    self.name = name
    self.spectral_count = spectral_count
    self.modifications, self.modification_count, self.masses, self.amino_acids, self.amino_acid_count = [], [], [], [], []
    self.identifiers, self.peptides, self.ratios = [], [], []
    self.peptide_count = 0
    self.dataset_dict, self.ratio_dict = {}, {}
    self.ratio_mean, self.ratio_median, self.ratio_stdev, self.ratio_median_mean_list = [], [], [], []
    self.ratio_median_mean = 0
    self.experiment_dataset_dict = {}
    self.dataset_dict = {}
    self.experiment_rep_dict = {}
    self.experiment_reps = 0
    self.dataset_reps = 0
    self.experiment_ratio_dict = {}
    self.experiment_stdev_dict = {}
    self.aggregate_median_mean, self.aggregate_median_mean_stdev = 0, 0
    for h in range(len(experiments)):
      self.experiment_dataset_dict[experiments[h]] = []
      self.experiment_rep_dict[experiments[h]] = 0
    for i in range(len(pep_files)):
      self.dataset_dict[pep_files[i]] = 0
      self.ratio_dict[pep_files[i]] = []
    self.ratio_mean_dict, self.ratio_median_dict, self.ratio_stdev_dict = self.ratio_dict.copy(), self.ratio_dict.copy(), self.ratio_dict.copy()
    self.update(peptide, spectral_count, modifications, amino_acids, masses, ratio, dataset, experiment_id)

  def update(self, peptide, spectral_count, modifications, amino_acids, masses, ratio, dataset, experiment_id):
    self.peptides = self.update_list(peptide, self.peptides)
    self.peptide_count = len(self.peptides)

    self.spectral_count = self.update_count(spectral_count, self.spectral_count)
    self.modifications = self.update_semicolon_list(modifications, self.modifications)
    self.modification_count  = len(self.modifications)
    self.masses = self.update_semicolon_list(masses, self.masses)
    self.amino_acids = self.update_semicolon_list(amino_acids, self.amino_acids)
    self.amino_acid_count = len(self.amino_acids)

    self.update_datasets(dataset, experiment_id)
    self.update_experiments(experiment_id)
    
    if (self.args.experiment != "identification"):
      self.update_ratios(ratio, dataset)
      self.update_median_means()
      self.update_aggregate_median_mean()

  def output(self):
    if (self.args.experiment == "identification"):
      return self.identifier + ',' + self.protein + ',' + self.description + ',' + self.name + ',' + self.list_to_string(self.peptides) + ',' + str(self.peptide_count) + ',' + str(self.spectral_count) + ',' + self.list_to_string(self.modifications) + ',' + str(self.modification_count) + ',' + self.list_to_string(self.amino_acids) + ',' + str(self.amino_acid_count) + ',' + self.list_to_string(self.masses) + ',' + self.list_to_string(self.replace([self.experiment_reps], False)) + ',' + self.list_to_string(self.replace([self.dataset_reps], False)) + ',' + self.dict_to_string(self.replace(self.experiment_rep_dict, True)) + ',' + self.dict_to_string(self.replace(self.dataset_dict, True)) + ','
    
    else:

      experiment_ratio_dict_str = self.dict_to_string(self.replace(self.experiment_ratio_dict, True))
      if '[]' in experiment_ratio_dict_str:
        experiment_ratio_dict_str = experiment_ratio_dict_str.replace('[]', '')

      return self.identifier + ',' + self.protein + ',' + self.description + ',' + self.name + ',' + self.list_to_string(self.peptides) + ',' + str(self.peptide_count) + ',' + str(self.spectral_count) + ',' + self.list_to_string(self.modifications) + ',' + str(self.modification_count) + ',' + self.list_to_string(self.amino_acids) + ',' + str(self.amino_acid_count) + ',' + self.list_to_string(self.masses) + ',' + self.list_to_string(self.replace([self.experiment_reps], False)) + ',' + self.list_to_string(self.replace([self.dataset_reps], False)) + ',' + self.dict_to_string(self.replace(self.experiment_rep_dict, True)) + ',' + self.dict_to_string(self.replace(self.dataset_dict, True)) + ','+ self.list_to_string(self.replace([self.aggregate_median_mean], False)) + ',' + self.list_to_string(self.replace([self.aggregate_median_mean_stdev], False)) + ',' + experiment_ratio_dict_str + ',' + self.dict_to_string(self.replace(self.experiment_stdev_dict, True)) + ','+ self.list_to_string(self.ratios) + ',' + self.dict_to_string(self.replace(self.ratio_median_dict, True)).replace('[]', '') + ','

class Peptide(Protein):
  def __init__ (self, peptide, charge_state, prob, localization, protein, identifier, description, spectral_count, intensity, miscleavages, modifications, amino_acids, masses, ratio, dataset, pep_files, args, experiment_id, experiments):
    self.args = args
    self.peptide = peptide
    self.unmodified_peptide = peptide.replace("*", "")
    self.spectral_count = spectral_count
    self.intensity = 0
    self.intensities = []
    self.miscleavages = miscleavages
    self.modifications, self.modification_count, self.masses, self.amino_acids, self.amino_acid_count = [], [], [], [], []
    self.charge_states, self.probabilities, self.localizations, self.proteins, self.identifiers, self.descriptions, self.ratios = [], [], [], [], [], [], []
    self.ratio_mean, self.ratio_median, self.ratio_stdev, self.ratio_median_mean_list = [], [], [], []
    self.ratio_median_mean = 0
    self.experiment_dataset_dict, self.dataset_dict, self.experiment_rep_dict = {}, {}, {}
    self.experiment_reps, self.dataset_reps, self.dataset_reps = 0, 0, 0
    self.ratio_dict, self.experiment_ratio_dict, self.experiment_stdev_dict = {}, {}, {}
    self.aggregate_median_mean, self.aggregate_median_mean_stdev = 0, 0
    for h in range(len(experiments)):
      self.experiment_dataset_dict[experiments[h]] = []
      self.experiment_rep_dict[experiments[h]] = 0
    for i in range(len(pep_files)):
      self.dataset_dict[pep_files[i]] = 0
      self.ratio_dict[pep_files[i]] = []
    self.ratio_mean_dict, self.ratio_median_dict, self.ratio_stdev_dict = self.ratio_dict.copy(), self.ratio_dict.copy(), self.ratio_dict.copy()
    self.update(charge_state, prob, localization, protein, identifier, description, spectral_count, intensity, modifications, amino_acids, masses, ratio, dataset, experiment_id)

  def update(self, charge_state, prob, localization, protein, identifier, description, spectral_count, intensity, modifications, amino_acids, masses, ratio, dataset, experiment_id):
    
    self.charge_states = self.update_semicolon_list(charge_state, self.charge_states)
    self.probabilities = self.update_semicolon_list(prob, self.probabilities)
    self.proteins = self.update_list(protein, self.proteins)
    self.descriptions = self.update_list(description, self.descriptions)
    self.identifiers = self.update_list(identifier, self.identifiers)
    self.localizations = self.update_list(localization, self.localizations)

    self.spectral_count = self.update_count(spectral_count, self.spectral_count)
    self.modifications = self.update_semicolon_list(modifications, self.modifications)
    self.modification_count  = len(self.modifications)
    self.masses = self.update_semicolon_list(masses, self.masses)
    self.amino_acids = self.update_semicolon_list(amino_acids, self.amino_acids)
    self.amino_acid_count = len(self.amino_acids)

    self.update_datasets(dataset, experiment_id)
    self.update_experiments(experiment_id)

    if (self.args.experiment != "identification"):
      self.update_ratios(ratio, dataset)
      self.update_median_means()
      self.update_aggregate_median_mean()
    
    if (self.args.experiment == "identification") | (self.args.experiment == "iodination"):
      self.intensities = self.update_semicolon_list(str(intensity), self.intensities)
    else:
      self.intensities = intensity

  def output(self):
    if (self.args.experiment == "identification"):
      return self.peptide + ',' + self.unmodified_peptide + ',' + self.list_to_string(self.charge_states) + ',' + self.list_to_string(self.probabilities) + ',' + self.list_to_string(self.localizations) + ',' + self.list_to_string(self.proteins) + ',' + self.list_to_string(self.identifiers) + ',' + self.list_to_string(self.descriptions) + ',' + str(self.spectral_count) + ',' + self.list_to_string(self.intensities) + ',' + str(self.miscleavages) + ',' + self.list_to_string(self.modifications) + ',' + str(self.modification_count) + ',' + self.list_to_string(self.amino_acids) + ',' + str(self.amino_acid_count) + ',' + self.list_to_string(self.masses) + ',' + self.list_to_string(self.replace([self.experiment_reps], False)) + ',' + self.list_to_string(self.replace([self.dataset_reps], False)) + ',' + self.dict_to_string(self.replace(self.experiment_rep_dict, True)) + ',' + self.dict_to_string(self.replace(self.dataset_dict, True)) + ','
    else:
      experiment_ratio_dict_str = self.dict_to_string(self.replace(self.experiment_ratio_dict, True))
      if '[]' in experiment_ratio_dict_str:
        experiment_ratio_dict_str = experiment_ratio_dict_str.replace('[]', '')

      if (self.args.experiment == "iodination"):
        return self.peptide + ',' + self.unmodified_peptide + ',' + self.list_to_string(self.charge_states) + ',' + self.list_to_string(self.probabilities) + ',' + self.list_to_string(self.localizations) + ',' + self.list_to_string(self.proteins) + ',' + self.list_to_string(self.identifiers) + ',' + self.list_to_string(self.descriptions) + ',' + str(self.spectral_count) + ',' + self.list_to_string(self.intensities) + ',' + str(self.miscleavages) + ',' + self.list_to_string(self.modifications) + ',' + str(self.modification_count) + ',' + self.list_to_string(self.amino_acids) + ',' + str(self.amino_acid_count) + ',' + self.list_to_string(self.masses) + ',' + self.list_to_string(self.replace([self.experiment_reps], False)) + ',' + self.list_to_string(self.replace([self.dataset_reps], False)) + ',' + self.dict_to_string(self.replace(self.experiment_rep_dict, True)) + ',' + self.dict_to_string(self.replace(self.dataset_dict, True)) + ','+ self.list_to_string(self.replace([self.aggregate_median_mean], False)) + ',' + self.list_to_string(self.replace([self.aggregate_median_mean_stdev], False)) + ',' + experiment_ratio_dict_str + ',' + self.dict_to_string(self.replace(self.experiment_stdev_dict, True)) + ','+ self.list_to_string(self.ratios) + ',' + self.dict_to_string(self.replace(self.ratio_median_dict, True)).replace('[]', '') + ','
      else:
        return self.peptide + ',' + self.unmodified_peptide + ',' + self.list_to_string(self.charge_states) + ',' + self.list_to_string(self.probabilities) + ',' + self.list_to_string(self.localizations) + ',' + self.list_to_string(self.proteins) + ',' + self.list_to_string(self.identifiers) + ',' + self.list_to_string(self.descriptions) + ',' + str(self.spectral_count) + ',' + str(self.intensity) + ',' + str(self.miscleavages) + ',' + self.list_to_string(self.modifications) + ',' + str(self.modification_count) + ',' + self.list_to_string(self.amino_acids) + ',' + str(self.amino_acid_count) + ',' + self.list_to_string(self.masses) + ',' + self.list_to_string(self.replace([self.experiment_reps], False)) + ',' + self.list_to_string(self.replace([self.dataset_reps], False)) + ',' + self.dict_to_string(self.replace(self.experiment_rep_dict, True)) + ',' + self.dict_to_string(self.replace(self.dataset_dict, True)) + ','+ self.list_to_string(self.replace([self.aggregate_median_mean], False)) + ',' + self.list_to_string(self.replace([self.aggregate_median_mean_stdev], False)) + ',' + experiment_ratio_dict_str + ',' + self.dict_to_string(self.replace(self.experiment_stdev_dict, True)) + ','+ self.list_to_string(self.ratios) + ',' + self.dict_to_string(self.replace(self.ratio_median_dict, True)).replace('[]', '') + ','

def start_game(args):
  start_time = time.time()

  cd = os.getcwd()

  if (args.database == "skyline"):
    uniprot_db = UniprotDatabase(args, cd)
    Game(args, cd, uniprot_db)
  else: 
    Game(args, cd, '')

  print("--- %s seconds ---" % (time.time() - start_time))
  # print('\n' + Quote.print() + '\n')

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-df', '--data_folder', dest='df', nargs='?', default="data", type=str, help='default data')
  parser.add_argument('-rf', '--results_folder', dest='rf', nargs='?', default="results", type=str, help='default results')
  parser.add_argument('-uf', '--uniprot_folder', dest='uf', nargs='?', default="/Users/lisamarieboatner/dropbox/backus/master", type=str, help='default "/Users/lisamarieboatner/dropbox/backus/master"')
  parser.add_argument('-sr', '--single_replicate', dest='single_replicate', nargs='?', default="False", type=str, help='default False')  
  parser.add_argument('-db', '--database', dest='database', nargs='?', default="msf", type=str, help='default msf; options = ip2, fp, skyline, mix')
  parser.add_argument('-dbv', '--database_version', dest='database_version', nargs='?', default="18", type=str, help='default 18; options = 17')
  parser.add_argument('-aa', '--amino_acid', dest='aa', nargs='?', default="C", type=str, help='default C; options = K, C;K, N-term')
  parser.add_argument('-pm', '--probe_mass', dest='probe_mass', nargs='?', default="464.28595", type=str, help='default 464.28595; options = 0')  
  parser.add_argument('-lpm', '--light_probe_mass', dest='light_probe_mass', nargs='?', default="463.2366", type=str, help='default 463.2366')  
  parser.add_argument('-hpm', '--heavy_probe_mass', dest='heavy_probe_mass', nargs='?', default="467.2529", type=str, help='default 467.2529')  
  parser.add_argument('-uniref', '--uniprot_ref', dest='uniref', nargs='?', default="072621_uniprot_reviewed.fasta", type=str, help='default 072621_uniprot_reviewed.fasta')
  parser.add_argument('-proref', '--protein_ref', dest='proref', nargs='?', default="protein.fas", type=str, help='default protein.fas')
  parser.add_argument('-pepref', '--pep_ref', dest='pepref', nargs='?', default="peptide.tsv", type=str, help='default peptide.tsv, peptide_label_quant.tsv, psm.tsv')
  parser.add_argument('-exp', '--experiment', dest='experiment', nargs='?', default="identification", type=str, help='default identification; options = isotop, silac, silac_iodination, iodination, intensity, lfq')
  parser.add_argument('-lhr', '--light_heavy_ratio', dest='light_heavy_ratio', nargs='?', default="light_heavy", type=str, help='default light_heavy; options = heavy_light')
  parser.add_argument('-rr', '--raw_ratio', dest='raw_ratio', nargs='?', default="False", type=str, help='default False; options = True')
  parser.add_argument('-ptm', '--ptm_prophet', dest='ptm', nargs='?', default="False", type=str, help='default False; options = True')
  parser.add_argument('-c', '--clean', dest='c', nargs='?', default="True", type=str, help='default True; options = False')
  parser.add_argument('-cpn', '--clean_peptide_number', dest='cpn', nargs='?', default="3", type=str, help='default 3; options = 1, 2, 3, etc.')
  parser.add_argument('-cen', '--clean_experiment_number', dest='cen', nargs='?', default="1", type=str, help='default 1; options = 1, 2, 3, etc.')
  parser.add_argument('-cmc', '--clean_miscleavages', dest='cmc', nargs='?', default="False", type=str, help='default False; options = True')
  parser.add_argument('-cmcn', '--clean_miscleavage_number', dest='cmcn', nargs='?', default="1", type=str, help='default 1; options = 1, 2, 3, etc.')
  parser.add_argument('-wo', '--write_outfile', dest='write_outfile', nargs='?', default="True", type=str, help='default True')
  parser.add_argument('-ce', '--compile_experiments', dest='compile_experiments', nargs='?', default="True", type=str, help='default True; options = False')
  parser.add_argument('-fin', '--entire_process', dest='entire_process', nargs='?', default="True", type=str, help='default True; options = False')

  args = parser.parse_args()
  start_game(args)

if __name__ == "__main__":
  main()
