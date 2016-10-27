
gfile="/Users/lab/Documents/Neta_Livne/genome_seqs/Ofas.scaffolds.fa" ## change directory
genome=open(gfile, "r")
lines=genome.readlines()
genome.close()
locus = raw_input("locus? in the format 'Scaffold5:1352665-1355673' , note the direction") ####if used as a function to be called in
###other script, del this line
def get_seq(locus):
##    flag1=True
##    while flag1:
        #locus="Scaffold5:1352665-1355673"
        aa=locus.split(":")
        scaff='>'+aa[0]+'\n'
        ab=aa[1]
        ac=ab.split("-")
        ad=int(ac[0])
        ae=int(ac[1])

        if ad<ae:
            beg=ad
            ends=ae
            flag=True #seq is on the positive strand
        else:
            beg=ae
            ends=ad
            flag=False #seq is on the negative strand

        for i in xrange(len(lines)): # i is the iterator- index!!! over the lines in the genome
        ##    if len(lines[i])=>51:
        ##    print i, len(lines[i])
            if lines[i].find(scaff)!=-1:
                lib=i+1+(beg/50) #number of line in the genome where the seq starts
                indb=beg%50 # location of the seq beginning within the line
                lie=i+1+(ends/50) 
                inde=ends%50 #location of the seq end within the line
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
        ##            print seq
                if flag==False: #seq on the - strand
                        seqr= seq[::-1]
                        seqt= seqr.replace('T','a')
                        seqa= seqt.replace('A', 't')
                        seqg= seqa.replace('G', 'c')
                        seqc= seqg.replace('C', 'g')
                        seq= seqc.upper()
        ##print "len(seq)" , len(seq)
        ##print "coordinates", ad , ae , "ad-ae", ad-ae
 #       print seq

        return seq
##ask=True
##while ask:
##    locus=raw_input("what's the locus? type it as 'ScaffoldX:Y-Z'")            
##        #print lines[i], lines[i+1][1:5]
print get_seq(locus)
##    
##    

