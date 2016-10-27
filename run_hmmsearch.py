import os
###### this script should be run via the terminal by the code:
############## python /Users/lab/Documents/Neta_Livne/run_hmmsearch.py

## GENERAL EXPLANATION:
##open the file with list of hmm files
##read it into a list
##iterate the list, for each item redefine the variables
##for each item run the os commanD
##saves a file with a list of the result files paths.

file_names_list = []
#hmm_file_list = open("/Users/lab/Documents/Neta_Livne/tf_domain_analysis/hmm_file_names_list_dna_binding.txt", "r")
hmm_file_list = open("/Users/lab/Documents/Neta_Livne/tf_domain_analysis/hmm_file_names_list_tf_new.txt", "r") #ADJUST FOLDER
######the file 'hmm_file_names_list_SOME RAW INPUT".txt' was
######created in the the script "fetching_hmm_from_http.py"
###### by using some input value/term given by the user


###### the var "hmm_file_list" should be a list of hmm *models* file names,
######e.g "PF00010_HLH_Helix-loop-helix_DNA-binding_domain.hmm"

for line in hmm_file_list:
    hmm_file = line.replace("\n", "")
    file_names_list.append(hmm_file)
hmm_file_list.close()
#print file_names_list

list_res_files = []

counter=1

algorithm = "/Users/lab/Downloads/hmmer-3.1b2-macosx-intel/src/hmmsearch"
parm1= " -o "
parm2 = " --incT 3.0"
ofas_fasta = " /Users/lab/Documents/Neta_Livne/tf_domain_analysis/cds_to_peptide_seq_10A_10P.fasta" ### ADJUST
#ofas_fasta = " /Users/lab/Downloads/prd_scaff2098_nt_seqs_for_motiff_search.fasta" ### ADJUST

for item in file_names_list:
    hmm_model = " /Users/lab/Documents/Neta_Livne/tf_domain_analysis/hmm_tf_search_term_files/"+item #ADJUST
    res_name = item.replace(".hmm", "")
    result_file = "/Users/lab/Documents/Neta_Livne/tf_domain_analysis/hmm_test_search_results_14_feb_16/res"+res_name + ".txt" #ADJUST FOLDER
    list_res_files.append(result_file) #ADJUST FOLDER
    command = algorithm+ parm1+ result_file + parm2 + hmm_model+ ofas_fasta
    print item
    print counter
    counter= counter+1
    os.system(command)

file_list_paths= open( "/Users/lab/Documents/Neta_Livne/tf_domain_analysis/file_res_list_paths_14_feb_16.txt", "w") #ADJUST
file_list_paths.write(str(list_res_files))
file_list_paths.close()
#####the file produced by the var: "file_list_paths" will be used in the
######script "sum_hmmsearch_res_to_rnaseq_try.py" as var: res_list_file


#hmm_model = " /Users/lab/Documents/Neta_Livne/tf_domain_analysis/hmm_dna_binding_files/PF00751_DM_DM_DNA_binding_domain.hmm"
##
##command = algorithm+ parm1+ result_file + parm2 + hmm_model+ ofas_fasta

print "hi"
