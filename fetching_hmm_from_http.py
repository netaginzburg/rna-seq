import urllib

#hmm_id_list_full = open("/Users/lab/Documents/Neta_Livne/tf_domain_analysis/dna_binding_domain_id_list.txt", "r")
term=raw_input("what was the hmm search term?")
hmm_id_list_full = open("/Users/lab/Documents/Neta_Livne/tf_domain_analysis/term_tf_domain_hmm_list.txt", "r")
#hmm_id_list= []
#hmm_file_name_list = []
hmm_file_name_list = open("/Users/lab/Documents/Neta_Livne/tf_domain_analysis/hmm_file_names_list_"+term+".txt", "w")
## the output of hmm_file_name_list will be used in the script "run_hmmsearch.py" as var "hmm_file_list"
hmm_id_list_txt = open("/Users/lab/Documents/Neta_Livne/tf_domain_analysis/hmm_id_list_"+term+".txt", "w")
counter=1
for line in hmm_id_list_full:
##    print line
    pf_id=line[:line.find("\t")]
    print pf_id
    print counter
    counter = counter+1
#    hmm_id_list.append(pf_id)
    pf_name= line[8:line.find("\t", 9)]
##    print pf_name
    hmm_id_list_txt.write(pf_id)
    hmm_id_list_txt.write("\n")
    pf_ful_name=line[line.find("\t", 9):]
    pf_full_name= pf_ful_name.replace(" ", "_")
    pfa_full_name=pf_full_name.replace("\n", "")
    pfam_full_name=pfa_full_name.replace("/", "_slash_")
#    print pf_full_name
#    hmm_file_name_list.append(pf_full_name)
    pf_file_nam=pf_id+"_"+pf_name+"_"+pfam_full_name+".hmm"
    pf_file_name= pf_file_nam.replace("\t", "")
##    hmm_file_name_list.write(pf_file_name)
##    hmm_file_name_list.write("\n")
    source_url="http://pfam.xfam.org/family/"+pf_id+"/hmm"
    output_hmm_file="/Users/lab/Documents/Neta_Livne/tf_domain_analysis/hmm_tf_search_term_files/"+pf_id+".hmm" ##ADJUST FOLDER
    hmm_file_name_list.write(pf_id+".hmm")
    hmm_file_name_list.write("\n")
#    print output_hmm_file
    testfile = urllib.URLopener() #extracting the hmm online

    testfile.retrieve(source_url, output_hmm_file)




#hmm_id_list_txt.write(hmm_id_list)
hmm_id_list_txt.close()

#hmm_file_name_list.write(hmm_file_name_list)
hmm_file_name_list.close()


#pf_id=
#pf_name=

##
##pf_id="PF1"
##pf_name="PAX"

##source_url="http://pfam.xfam.org/family/"+pf_id+"/hmm"
##output_hmm_file="/Users/lab/Documents/Neta_Livne/"+pf_id+"_"+pf_name+".hmm"
##
##testfile = urllib.URLopener()
##testfile.retrieve(source_url, output_hmm_file)


######
#testfile.retrieve("http://pfam.xfam.org/family/PF00096/hmm", "/Users/lab/Documents/Neta_Livne/zf-C2H2.hmm")


#./hmmsearch -o /Users/lab/Documents/Neta_Livne/tf_domain_analysis/hmm_search_results/outputfile1.txt --incT 3.0 /Users/lab/Documents/Neta_Livne/tf_domain_analysis/hmm_dna_binding_files/ /Users/lab/Downloads/prd_scaff2098_nt_seqs_for_motiff_search.fasta
