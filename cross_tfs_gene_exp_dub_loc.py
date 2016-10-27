import xlrd
import numpy

#tfile= xlrd.open_workbook('/Users/lab/Documents/Neta_Livne/test8.xlsx')
tfile= xlrd.open_workbook('/Users/lab/Documents/Neta_Livne/tfs_min_overlaps_new.xlsx')
sheet2 = tfile.sheet_by_index(0)
print 'screen tf table'
tf_gene_list_pos = []
tf_gene_list_neg = []
for k in range(sheet2.nrows):
    print k
    Dmel_gene = str(sheet2.cell_value(k, 0))
    GB_iso = sheet2.cell_value(k, 1)
    tf_scaff = str(sheet2.cell_value(k, 3))
#    print tf_scaff
    tf_beg = sheet2.cell_value(k, 4)
    tf_end = sheet2.cell_value(k, 5)
    tf_dir = numpy.sign(sheet2.cell_value(k, 6))
    if tf_dir ==1:
        tf_gene={ 'Dmel_gene':Dmel_gene, 'GB_iso': int(GB_iso), 'tf_scaff':tf_scaff ,
                  'tf_beg': tf_beg ,'tf_end':tf_end ,'tf_dir':tf_dir, 'k':k}
        tf_gene_list_pos.append(tf_gene)
    if tf_dir ==-1:
        tf_gene={ 'Dmel_gene':Dmel_gene, 'GB_iso': int(GB_iso), 'tf_scaff':tf_scaff ,
                  'tf_beg': tf_beg ,'tf_end':tf_end ,'tf_dir':tf_dir, 'k':k}
        tf_gene_list_neg.append(tf_gene)
tfile.release_resources()
#expdata= xlrd.open_workbook()
rna_tables= xlrd.open_workbook('/Users/lab/Documents/Neta_Livne/Ofas_res_new_feb15(2)/gene_exp_OFAS_trial.xls')

from Bio import Entrez
Entrez.email = 'neta.livne88@gmail.com'
#rna_tables= xlrd.open_workbook('/Users/lab/Documents/Neta_Livne/test20.xlsx')
for r in rna_tables.sheet_names(): ####??? #iterate over exp data workbook
    sheet1 = rna_tables.sheet_by_name(r) #exp_data_sheet
    summ_file=open('/Users/lab/Documents/Neta_Livne/new_tf_files/exp_cross_tfs_'+r+'.txt', 'w')
    print 'screen expdata', r
    j=2
    exp_gene_list_pos = []
    exp_gene_list_neg = []
    exp_gene_list_non_direction = [] #locuses marked with ?
    while j <len(range(sheet1.nrows)): #iterate over lines in exp table sheet, create a list
        print j
        test_id= str(sheet1.cell_value(j, 0))
#        for v in corrected_loc_list: 
        locus1 = sheet1.cell_value(j, 3)
        if "?" in locus1: # gene should be add to a separate list
            print "?"
            locusdub=locus1.split('?')
            locusdu = locusdub[1]
            locus2 = locusdu.split(':')
            exp_scaff = locus2[0]
            locus3 = locus2[1].split('-')
            exp_be= int(locus3[0])
            exp_en = int(locus3[1])
            exp_gene_p = {'j':j, 'test_id': test_id, 'exp_scaff':exp_scaff, 'exp_beg':exp_be,
                        'exp_end':exp_en, 'exp_dir':'dub_plus'}
            exp_gene_list_pos.append(exp_gene_p)
            exp_gene_m = {'j':j, 'test_id': test_id, 'exp_scaff':exp_scaff, 'exp_beg':exp_be,
                        'exp_end':exp_en, 'exp_dir':'dub_minus'}
            exp_gene_list_neg.append(exp_gene_m)
        else:
            locus2 = locus1.split(':')
            exp_scaff = locus2[0]
            locus3 = locus2[1].split('-')
            exp_be= int(locus3[0])
            exp_en = int(locus3[1])
 #           print exp_be, exp_en
    #    true_loc = str(sheet.cell_value(j, 0)) ## change table index? ##if neccesary?
            if exp_be<exp_en: #check the direction of seq- farward/reverse
                exp_dir = 1
                exp_beg= exp_be
                exp_end= exp_en
                exp_gene = {'j':j, 'test_id': test_id, 'exp_scaff':exp_scaff, 'exp_beg':exp_beg,
                            'exp_end':exp_end, 'exp_dir':exp_dir}
 #               print exp_gene
                exp_gene_list_pos.append(exp_gene)
            if exp_en<exp_be: #fliping the reverse reading and tagging as
                exp_dir = -1
                exp_beg=exp_en
                exp_end=exp_be
 #               print exp_beg, exp_end
                exp_gene = {'j':j, 'test_id': test_id, 'exp_scaff':exp_scaff, 'exp_beg':exp_beg,
                            'exp_end':exp_end, 'exp_dir':exp_dir}
#                print exp_gene
                exp_gene_list_neg.append(exp_gene)
        j+=1
#    print "exp_gene_list_neg", exp_gene_list_neg
#    print "exp_gene_list_pos", exp_gene_list_pos
    xloc_list=[] #???
    cross_pos_list=[]
    cross_dub_list=[]
    cross_dun_ids = [] #for comparission and matching minus and plus entries
    united_dub = [] #for united genes with no orientation
    print 'cross positive lists', r
    for n in exp_gene_list_pos:
        print 'n', n['j']
        for p in tf_gene_list_pos:
            #print 'p', p['k']
            if n['exp_scaff']!= p['tf_scaff']:
    ##            print 'not same scaff'
                continue
            if n['exp_scaff']== p['tf_scaff']:
    ##            print 'same scaff'
                if (n['exp_beg']>p['tf_end']) or (n['exp_end']<p['tf_beg']): #no match
    #                print n['exp_beg'], n['exp_end'], p['tf_beg'], p['tf_end'] 
    ##                print 'no match'
                    continue
                else:
    ##                print n
                    min_loc=min(n['exp_beg'], n['exp_end'], p['tf_beg'], p['tf_end'])
                    #(tf_end,tf_beg, exp_end, exp_beg)
                    max_loc=max(n['exp_beg'], n['exp_end'], p['tf_beg'], p['tf_end'])
    ##                print n['exp_beg'], n['exp_end'], p['tf_beg'], p['tf_end']
    ##                print max_loc, min_loc
                    x=max_loc-min_loc
                    overlap_loc=n['exp_end']-n['exp_beg']-p['tf_beg']+p['tf_end']-x
                    over_ratio=overlap_loc/x
                    lexp=n['exp_end']-n['exp_beg']
                    ltf=p['tf_end']-p['tf_beg']
                    tf_cover=overlap_loc/ltf
                    exp_cover=overlap_loc/lexp
                    n['x']=x
                    n['overlap_length']=overlap_loc
                    n['overlap_ratio']=over_ratio
                    n['tf_cover']=tf_cover
                    n['exp_cover']=exp_cover
                    n['gb_id']=p['GB_iso']
                    handle = Entrez.efetch(db="protein", id=str(n['gb_id']), rettype="gp", retmode="xml")
                    record = Entrez.read(handle)
                    prot_details=[]
                    prot_details.append(record[0]['GBSeq_definition'])
                    for i in record[0]['GBSeq_feature-table']:
                        # print len(i.items())
                         for itdict in i.items():
                             for ydict in itdict:
                                 if isinstance(ydict, unicode):
                                     #print 'unicode'
                                     continue
                                 if str(type(ydict))=="<class 'Bio.Entrez.Parser.StringElement'>":
                ##                     print 'Bio.Entrez.Parser.StringElement'
                                     continue
                                 if str(type(ydict))=="<class 'Bio.Entrez.Parser.ListElement'>":
                ##                     print 'Bio.Entrez.Parser.ListElement'
                                     list_dict=str(ydict)
                                     if 'region_name' in list_dict:
                                         for k in ydict:
                                             prot_details.append(k.values())

                    prot_details1= str(prot_details)
                    prot_details2= prot_details1.translate(None, ',')
                    prot_details3= prot_details2.translate(None, ':')
                    n['prot_details2']=prot_details3
                    if n['exp_dir']==1:
                        cross_pos_list.append(n)
                    if n['exp_dir']=='dub_plus':
                        cross_dub_list.append(n)
                        cross_dun_ids.append(n['test_id'])
                    xloc_list.append(n['test_id']) #???
    ##                print 'appeding'
    ##                print n
    ##                print 'break'
                    break

    cross_neg_list=[]
    print 'cross negative lists'
 #   print exp_gene_list_neg
    for g in exp_gene_list_neg:
        print 'n', g['j'], g['test_id']
        for y in tf_gene_list_neg:
#            print 'y', y['k']
#            print "g['exp_scaff']",g['exp_scaff'], "y['tf_scaff']", y['tf_scaff']
            if g['exp_scaff']!= y['tf_scaff']:
#                print 'not same scaff'
                continue
            if g['exp_scaff']== y['tf_scaff']:
#                print 'same scaff'
 #               print g['exp_beg'], g['exp_end'], y['tf_beg'], y['tf_end'] 
                if (g['exp_beg']>y['tf_end']) or (g['exp_end']<y['tf_beg']): #previous- wrong ?
#                    print g['exp_beg'], g['exp_end'], y['tf_beg'], y['tf_end'] 
                    print 'no match'
                    continue
                else:
 #                   print g
                    min_loc=min(g['exp_beg'], g['exp_end'], y['tf_beg'], y['tf_end'])
                    max_loc=max(g['exp_beg'], g['exp_end'], y['tf_beg'], y['tf_end'])
                    x=max_loc-min_loc
                    overlap_loc=g['exp_end']-g['exp_beg']-y['tf_beg']+y['tf_end']-x
                    over_ratio=overlap_loc/x
                    lexp=g['exp_end']-g['exp_beg']
                    ltf=y['tf_end']-y['tf_beg']
                    tf_cover=overlap_loc/ltf
                    exp_cover=overlap_loc/lexp
                    g['x']=x
                    g['overlap_length']=overlap_loc
                    g['overlap_ratio']=over_ratio
                    g['tf_cover']=tf_cover
                    g['exp_cover']=exp_cover
                    g['gb_id']=y['GB_iso']
                    handle = Entrez.efetch(db="protein", id=str(g['gb_id']), rettype="gp", retmode="xml")
                    record = Entrez.read(handle)
                    prot_details=[]
                    prot_details.append(record[0]['GBSeq_definition'])
                    for i in record[0]['GBSeq_feature-table']:
                        # print len(i.items())
                         for itdict in i.items():
                             for ydict in itdict:
                                 if isinstance(ydict, unicode):
                                     continue
                                 if str(type(ydict))=="<class 'Bio.Entrez.Parser.StringElement'>":
                                     continue
                                 if str(type(ydict))=="<class 'Bio.Entrez.Parser.ListElement'>":
                                     list_dict=str(ydict)
                                     if 'region_name' in list_dict:
                                         for k in ydict:
                                             prot_details.append(k.values())
                    prot_details1= str(prot_details)
                    prot_details2= prot_details1.translate(None, ',')
                    prot_details3= prot_details2.translate(None, ':')
                    g['prot_details2']=prot_details3

                    if g['exp_dir']==-1:
                        cross_neg_list.append(g)
                        print "append neg_list",g['test_id']
                    if g['exp_dir']=='dub_minus':
                        cross_dub_list.append(g)
                        cross_dun_ids.append(g['test_id'])
                    xloc_list.append(g['test_id'])
                    print 'appeding', g['test_id']
 #                   print g
    ##                print 'break'
                    break
    
    print "going over genes with no orientation"
### all the genes that had a hit and a "?" in the scaffold were added to 2 lists:
####cross_dub_list=list dicts of each gene's dictionary ;
####and "cross_dun_ids" =contains only the gene id. both lists may contain the same gene twice,
####if it has 2 hits.
####this part goes over the dicts list. on each iteration it checks if the id of the item appears
####twice in the ids list, and if so, it deletes the first occourance of it from the id's list.
####then, it takes the index of the remaining id in the ids list, and 1 to it
####(since the dict list is in the same length as it was when we entered this iteration)
####and look for the item with that index(sec_index) in the dicts list.
####the entire dictionary of the second hit is add to the protein details of the first hit, and then
####the dict is add to united_dub list, which is the final list of doubles, and the second dict of the gene
####is deleted from the dicts list so it won't reocour in the file. if the gene occours only once in the ids
####list- no problem, it will be add to the final list.
 #   print cross_dub_list
    for q in cross_dub_list:
        print q['test_id']
        print "cross_dun_ids", cross_dun_ids
        if cross_dun_ids.count(q['test_id'])==2: #each orientaton has a different homolog
            print "HEY!"
 #           print "before del q cross_dun_ids", cross_dun_ids
            cross_dun_ids.remove(q['test_id'])
#            print "after del q cross_dun_ids", cross_dun_ids
            sec_index = cross_dun_ids.index(q['test_id'])+1
#            print "sec_index", sec_index
            seqq_index = cross_dub_list[sec_index]

            secc_index= "Equivalent String : %s" % str(seqq_index)
            seccc_index =secc_index.translate(None, ":")
            secccc_index=seccc_index.translate(None, ",")
#            print "seccc_index", secccc_index
#            print "str(q['prot_details2'])", str(q['prot_details2'])
            baaaa=  str(q['prot_details2'])+secccc_index
#            print "baaaa", baaaa
#            print type(baaaa)
            q['prot_details2']=baaaa
            #print 'q', q
            united_dub.append(q)
#            print "#appending === each orientaton has a different homolog", q['test_id']
            del cross_dub_list[sec_index]
 #           del cross_dun_ids[sec_index]
#            print "cross_dun_ids", cross_dun_ids
#            print "cross_dub_list", cross_dub_list
            continue
        if cross_dun_ids.count(q['test_id'])==1: #only one orientation has a match
            united_dub.append(q)
#            print "only one orientation has a match ---appending" , q['test_id']
            ###
    print "adding genes with no matches"
    no_match=[]
    for t in exp_gene_list_pos:
        if t['test_id'] not in xloc_list:
            t['overlap_length']='no tf match'
            t['x']=0
            t['overlap_ratio']=0
            t['tf_cover']=0
            t['exp_cover']=0
            t['gb_id']=0
            t['prot_details2']='none'
            print t['test_id']
            no_match.append(t)
            xloc_list.append(t['test_id'])
    for m in exp_gene_list_neg:
        if m['test_id'] not in xloc_list:
            m['overlap_length']='no tf match'
            m['x']=0
            m['overlap_ratio']=0
            m['tf_cover']=0
            m['exp_cover']=0
            m['gb_id']=0
            m['prot_details2']='none'
            print m['test_id']
            no_match.append(m)
            
    print "saving file"
    print "cross_pos_list"
    for l in cross_pos_list:####
        print l['test_id']
        summ_file.write(str(l))
        summ_file.write('\n')
    print "cross_neg_list"
    for l in cross_neg_list:####
        print l['test_id']
        summ_file.write(str(l))
        summ_file.write('\n')
    print "united"
    for l in united_dub:####
        print l['test_id']
        summ_file.write(str(l))
        summ_file.write('\n')
    print "no match"
    for l in no_match:####
        print l['test_id']
        summ_file.write(str(l))
        summ_file.write('\n')

    summ_file.close()
    rna_tables.unload_sheet(r)

