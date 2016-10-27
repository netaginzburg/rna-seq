import xlrd

#expdata= xlrd.open_workbook('/Users/lab/Documents/Neta_Livne/test9.xlsx') #test file
#expdata= xlrd.open_workbook('/Users/lab/Documents/Neta_Livne/local_blast_results2.xlsx')
expdata= xlrd.open_workbook('/Users/lab/Documents/Neta_Livne/tfs_no_dou_rf.xlsx') 
sheet = expdata.sheet_by_index(0)
locus_list=[]
locus_list2=[]
locus_list3=[{'j':-1,'Dmel_gene':1, 'Dmel_isoform_gene_bank_ID':1, 'e_value':0.05,
              'strand':10, 'Ofas_scaffold':1, 'Ofas_start':1, 'Ofas_end':2}]
locus_list4=[]
summ_file=open('/Users/lab/Documents/Neta_Livne/no_double_hits_tfs3.txt', 'w')
for j in range(sheet.nrows):
    print j
    Dmel_gene = str(sheet.cell_value(j, 0))
    Dmel_isoform_gene_bank_ID1 = sheet.cell_value(j, 1)
    if type(Dmel_isoform_gene_bank_ID1) == float:
        Dmel_isoform_gene_bank_ID2 =int(Dmel_isoform_gene_bank_ID1)
        Dmel_isoform_gene_bank_ID =str(Dmel_isoform_gene_bank_ID2)
    else:
        Dmel_isoform_gene_bank_ID = str(Dmel_isoform_gene_bank_ID1)
    e_value = sheet.cell_value(j, 2)
    Ofas_scaffold = str(sheet.cell_value(j, 3))
    Ofas_start1 = sheet.cell_value(j, 4)
    if type(Ofas_start1)==float:
        Ofas_start2 = int(Ofas_start1)
        Ofas_start = str(Ofas_start2)
    else:
        Ofas_start = str(Ofas_start1)
    Ofas_end1 = sheet.cell_value(j, 5)
    if type(Ofas_end1)==float:
        Ofas_end2 = int(Ofas_end1)
        Ofas_end = str(Ofas_end2)
    else:
        Ofas_end = str(Ofas_end1)        
    strand = sheet.cell_value(j, 6)########
    gene_locus = {'j':j, 'Dmel_gene':Dmel_gene, 'Dmel_isoform_gene_bank_ID':Dmel_isoform_gene_bank_ID,
                  'e_value':e_value, 'strand':strand, 'Ofas_scaffold':Ofas_scaffold,
                  'Ofas_start':Ofas_start, 'Ofas_end':Ofas_end}
    if gene_locus not in locus_list:
            locus_list.append(gene_locus)
            locus_list2.append(gene_locus)

print "list comparssions"
   
for k in locus_list:
    print 'k', k['j']
    for n in locus_list3:     
##            print 'n', n['j']#
            if k['Ofas_scaffold']!=n['Ofas_scaffold']:
                flag= True
##                print 'flag scaffold',  flag#
            else:
                flag = False
                if k['strand']!=n['strand']:
                    flag = True
##                    print 'flag strand', flag#
                if k['strand']==n['strand']:
##                    print 'flag strand',  flag#
                    flag = False
                    if (eval(k['Ofas_start'])<= eval(n['Ofas_start'])) and (eval(k['Ofas_end'])>= eval(n['Ofas_end'])):
                            locus_list3.remove(n)
                            locus_list3.append(k)
##                            print "remove" , n['j'], 'append' , k['j']#
                    elif (eval(k['Ofas_start'])>= eval(n['Ofas_start'])) and (eval(k['Ofas_end'])<= eval(n['Ofas_end'])):
                        flag = False
##                        print "inner1 flag", flag#
                        break
                    else:
                        flag= True
##                    print "inner flag", flag#
    if flag==True:
##            print 'end of iter flag', flag#
            locus_list3.append(k)
##            print 'append k' , k['j']#

print 'del doubles'
for y in locus_list3:
    if y not in locus_list4:
        locus_list4.append(y)
        print y['j']
       
print "saving file"            
for l in locus_list4:
    print l['Ofas_scaffold']
    summ_file.write(str(l))
    summ_file.write('\n')

summ_file.close()
