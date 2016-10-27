import xlrd
gfile="/Users/lab/Documents/Neta_Livne/genome_seqs/Ofas.scaffolds.fa"
genome=open(gfile, "r")
lines=genome.readlines()
#locus="Scaffold5:1352665-1355673"
#rna_tables= xlrd.open_workbook('/Users/lab/Documents/Neta_Livne/Ofas_res_updated/gene_exp_OFAS_some_false_locuses.xls')
rna_tables= xlrd.open_workbook('/Users/lab/Documents/Neta_Livne/Ofas_res_updated/Sig_10A_19A_for_fasta.xlsx')
for t in rna_tables.sheet_names(): ####???
    sheett = rna_tables.sheet_by_name(t)
    outfile10a=open('/Users/lab/Documents/Neta_Livne/SEQFILE_10a'+t+'.fasta', 'w') #upregulated an 10a
    outfile10p=open('/Users/lab/Documents/Neta_Livne/SEQFILE_19a'+t+'.fasta', 'w')# upregulated at 10p
    u=2
    while u <len(range(sheett.nrows)):
        print u
        test_idu= str(sheett.cell_value(u, 0))
        locus= str(sheett.cell_value(u, 3))
        updown= int(sheett.cell_value(u, 53))
        if locus.find('?')!= -1:
            print locus
            aa=locus.split(":")
            scaff='>'+aa[0]+'\n'
            ab=aa[1]
            ac=ab.split("-")
            ad=int(ac[0])
            ae=int(ac[1])

            if ad<ae:
                flag= True
                beg=ad
                ends=ae
            else:
                flag=False
                beg=ae
                ends=ad
            #lines=genome.readlines()
            for i in xrange(len(lines)):
                if lines[i].find(scaff)!=-1:
                    lib=i+1+(beg/50)
                    indb=beg%50
                    lie=i+1+(ends/50)
                    inde=ends%50
                    if lib<lie:
                        seq=lines[lib][(indb-1):]
            ##            print seq
                        k=lib+1
                        while k<lie:
                            seq=seq+lines[k]
            ##                print seq
                            k+=1
                        seq=seq+lines[lie][:(inde)]
                    if lib==lie: #for tzach's bug- short sequences
                        seq=lines[lib][(indb-1):(inde)]
                    if flag==False: #seq on the - strand
                        seqr= seq[::-1]
                        seqt= seqr.replace('T','a')
                        seqa= seqt.replace('A', 't')
                        seqg= seqa.replace('G', 'c')
                        seqc= seqg.replace('C', 'g')
                        seq= seqc.upper()
##                    print flag                       
##                    print seq
       #             outfile= open('/Users/lab/Documents/Neta_Livne/rnaseq_seqs/seq'+aa[0]+'@'+ab+'.fasta', 'w')
                    seq1=seq.translate(None,'\n')
                    lenseq1= len(seq1)
                    lenseq2= lenseq1/70
                    lenseq3= lenseq1%70
##                    header= '>'+locus
                    if updown==1: #upregulated at 10a
##                        print '#upregulated at 10a'
                        header= '>upregulated_at_10a|'+test_idu+'|'+locus
                        outfile10a.write(header)
                        outfile10a.write('\n')
                        n=0
                        while n<lenseq1:
##                            print "writing"
                            seqline=seq1[n:n+70]
                            outfile10a.write(seqline)
                            outfile10a.write('\n')
                            n=n+70
                    if updown==0: #upregulated at 10p
##                        print '#upregulated at 10p'
                        header= '>upregulated_at_19a|'+test_idu+'|'+locus
                        outfile10p.write(header)
                        outfile10p.write('\n')
                        n=0
                        while n<lenseq1:
##                            print "writing"
                            seqline=seq1[n:n+70]
                            outfile10p.write(seqline)
                            outfile10p.write('\n')
                            n=n+70                   

                                       
                    #print lines[i], lines[i+1][1:5]
                
        u+=1    
    outfile10a.close()
    outfile10p.close()
