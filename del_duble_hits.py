import xlrd

#expdata= xlrd.open_workbook('/Users/lab/Documents/Neta_Livne/test10.xlsx') #test file
#expdata= xlrd.open_workbook('/Users/lab/Documents/Neta_Livne/local_blast_results2.xlsx')
expdata= xlrd.open_workbook('/Users/lab/Documents/Neta_Livne/tfs_withdoub.xlsx')
sheet = expdata.sheet_by_index(0)
locus_list=[]
summ_file=open('/Users/lab/Documents/Neta_Livne/no_double_best_hits_tfs.txt', 'w')
for j in range(sheet.nrows):
    print j
    Dmel_gene = str(sheet.cell_value(j, 0))
#    Dmel_isoform_gene_bank_ID1 = sheet.cell_value(j, 1)
#    if type(Dmel_isoform_gene_bank_ID1) == float:
#        Dmel_isoform_gene_bank_ID2 =int(Dmel_isoform_gene_bank_ID1)
#        Dmel_isoform_gene_bank_ID =str(Dmel_isoform_gene_bank_ID2)
#    else:
#        Dmel_isoform_gene_bank_ID = str(Dmel_isoform_gene_bank_ID1)
    e_value = sheet.cell_value(j, 1)
    Ofas_scaffold = str(sheet.cell_value(j, 2))
    Ofas_start1 = sheet.cell_value(j, 3)
    if type(Ofas_start1)==float:
        Ofas_start2 = int(Ofas_start1)
        Ofas_start = str(Ofas_start2)
    else:
        Ofas_start = str(Ofas_start1)
    Ofas_end1 = sheet.cell_value(j, 4)
    if type(Ofas_end1)==float:
        Ofas_end2 = int(Ofas_end1)
        Ofas_end = str(Ofas_end2)
    else:
        Ofas_end = str(Ofas_end1)        
#    best_hit = sheet.cell_value(j, 6)
#    print e_value
    strand = str(sheet.cell_value(j, 5))
    gene_locus = str(Ofas_scaffold+','+ Ofas_start+ ','+Ofas_end+','+strand)
#        gene_locus=str(sheet.cell_value(j, 4)+sheet.cell_value(j, 5)+sheet.cell_value(j, 6))
    if gene_locus not in locus_list:
            locus_list.append(gene_locus)
#            print locus_list
   #         print e_value
            save_to_table=['Dmel_gene', Dmel_gene ,'e_value', e_value, 'Ofas_scaffold', Ofas_scaffold, 'Ofas_start', Ofas_start, 'Ofas_end', Ofas_end,'strand', strand]#, best_hit] 
            summ_file.write(str(save_to_table))
            summ_file.write('\n')

summ_file.close()
