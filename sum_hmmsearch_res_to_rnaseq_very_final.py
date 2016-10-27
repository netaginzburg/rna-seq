##open the rna-seq excel
##read the files from the hmm-results folder
##for every file, save of the results

xloc_dict= {}

res_list_file= open("/Users/lab/Documents/Neta_Livne/tf_domain_analysis/file_res_list_paths_14_feb_16.txt", "r")
#the file "file_list_paths_test.txt" was created in the script "run_hmmsearch.py" as var: file_list_paths
    
a=res_list_file.read()
res_list=eval(a)
res_list_file.close()

count=1
model_txt_list=[]
print "open hmmsearch results files and copy it into a list, each item in the list is a list of all the text in a certain res file"
for model_obj in res_list: #open hmmsearch results files
##    print model_obj
    model_item = open(model_obj, "r")
    model_lines = model_item.readlines()
    model_txt_list.append(model_lines)
##    print model_lines
    print count
    count= count+1
#    print "\n"
    model_item.close()

no_hit= '   [No hits detected that satisfy reporting thresholds]\n'
end_table_mark ='Domain annotation for each sequence (and alignments):\n'
print "looking into the results. if there are hits these are added to a dict"
####ADD COUNTER!!!
for model in model_txt_list:
    if no_hit in model:
        continue
    else:
        hmm_query = model[11].split()[1]
        hmm_acces = model[12].split()[1]
        hmm_description = (model[13].split(":")[1]).replace(" ", "_")
        print hmm_description
##        hmm_description = model[13].split()[1]+"_"+model[13].split()[2]
        hmm_title = hmm_query +"~"+ hmm_acces + "~"+ hmm_description ### Maybe add revachim? effects the way 
                                                            ### the file will be imported to excel
        
        beg_table1 = '  ------ inclusion threshold ------\n'
        beg_table2 = '    ------- ------ -----    ------- ------ -----   ---- --  --------                                      -----------\n'
        if beg_table1 in model:
            hits_lines= model[model.index(beg_table1)+1:model.index(end_table_mark)]
        else:
            hits_lines= model[18:model.index(end_table_mark)]

##        beg_hits_table1 = model.index('  ------ inclusion threshold ------\n')
        
##        end_hits_table = model.index(end_table_mark)
##        hits_lines= model[beg_hits_table+1:end_hits_table]
##        print hits_lines
        for hit_item in hits_lines:
            hit_params= hit_item.split()
            hit_str= hmm_title+hit_item #var to be added to the dict- contains parameters related to the hit
                                        ##and hmm model details. Saved as str that can be later separeted to
                                        ###list by str.split() with empty paranthesis
            print hit_params
##            print "\n"
            if len(hit_params)>0:
                hit_name = hit_params[8] ### fasta files might have different names. WHY? because I was lazy... 
                xloc = hit_name[:11]#hit_name[hit_name.find("XLOC"):hit_name.find("|")]
                xloc_det = hit_name[12:]
                hit_str= xloc_det+hit_str
                print hit_str
                print xloc
                if xloc in xloc_dict:
                    xloc_dict[xloc].append(hit_str)
                if xloc not in xloc_dict:
                    xloc_dict[xloc]=[hit_str]
##                print xloc_dict
##
##                    
print "going over the xloc list excel file"
import xlrd
import xlwt
from tempfile import TemporaryFile

xloc_input = xlrd.open_workbook('/Users/lab/Documents/Neta_Livne/10a10p_comp_gene_exp_hmmsearch.xlsx')

rel_sheet= xloc_input.sheet_by_index(0)
print rel_sheet

book = xlwt.Workbook()
sheet1 = book.add_sheet('Sheet 1')

for i in range(rel_sheet.nrows):
    print i
    xloc_in= str(rel_sheet.cell_value(i, 0))
##    print xloc_in
    sheet1.write(i,0,xloc_in)
    amuda_num=1
    if xloc_in in xloc_dict:
        for j in xloc_dict[xloc_in]: ###Adjust to fasta title
            
            print j
            all_params = j.split()
            hited_hmm= all_params[0]
            fs_e_val=float(all_params[1])
            fs_score= float(all_params[2])
            fs_bias= float(all_params[3])
            reading_frame1= all_params[9].split("|")[-1]
            reading_frame= reading_frame1[reading_frame1.find("_", 5):]
            if fs_bias==0:
                order_fs="ok"
            else:
                order_fs= fs_score/fs_bias
##            all_params
            print all_params
            sheet1.write(i,amuda_num, hited_hmm)
            sheet1.write(i,amuda_num+1, fs_e_val)
            sheet1.write(i,amuda_num+2, reading_frame)

            amuda_num= amuda_num+3

###### (SHURA, AMUDA)

book.save('/Users/lab/Documents/Neta_Livne/hmm_res_10a10p.xls')
book.save(TemporaryFile())



############################################################

##xloc_input = xlrd.open_workbook('/Users/lab/Documents/Neta_Livne/tf_domain analysis input_test_file.xlsx')
                
##for i in some_list[40]:
##    print i
##    if "[No hits detected that satisfy reporting thresholds]" in i:
##        print "********NO HITS!!********"
##        continue
##    else:
##        print "go on"
##
##
##count = 0
##for j in model_txt_list[40]:
##    print count
##    count = count+1
##    print j
##    
