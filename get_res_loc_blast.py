import os
import numpy
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast.Applications import NcbitblastnCommandline
import math

#genome = "/home/barbara/data/genomes/Clectularius/Clec_Bbug02212013.genome.fa"
genome = '/Users/lab/Documents/Neta_Livne/genome_seqs/Ofas.scaffolds.fa'
#genome = "/home/barbara/data/genomes/Ofasciatus/Ofas.scaffolds.fa"
#species = "clec"
species = "ofas"
blastfolder = '/Users/lab/Documents/Neta_Livne/tfs_entrez_seqs/'
#blastfolder = "/home/barbara/Dropbox/oncopeltus/annotations/"

#f1 = open('/Users/lab/Documents/Neta_Livne/TFs_pool/test1.rtf', 'r')
f1 = open('/Users/lab/Documents/Neta_Livne/TFs_pool/file_names_list.rtf', 'r')
f2 = f1.read()
f3 = f2.split(' ')
f4 = f3[len(f3)-1].split('}')
filelist = f4[0].split(',')

summ_hits = open("/Users/lab/Documents/Neta_Livne/TFs_pool/local_summery_hits5.txt", "w")
#summ_hits.write('file_name \t max_e_val \t mean_e_val \t std \t median') order of parameters in the table
#filelist = ['FBtr0070075|FBgn0000137|ase', 'FBtr0070265|FBgn0000210|br', 'FBtr0070261|FBgn0000210|br', 'FBtr0070599|FBgn0021738|Crg-1'] #change to overlap_list
add = 1000
for i in filelist:
    outf = '/Users/lab/Documents/Neta_Livne/TFs_pool/local_blast_res_files/'+ i + "-" + species
    #tblastn_cline = NcbitblastnCommandline(query=blastfolder+i, db=genome, outfmt=5, out= outf+'.xml')
    #blast_record = NCBIXML.read(out)
    #stdout, stderr = tblastn_cline()
    #os.system(str(tblastn_cline))
    result_handle = open(outf+'.xml')
    blast_record = NCBIXML.read(result_handle)
    E_VALUE_THRESH = 0.04 ### change and fine-tune
    print i
    e_val_list = []
    expect_min=100
    expect_max=0.05 #maximum e-val for a hit
    if len(blast_record.alignments)> 0:
        for alignment in blast_record.alignments:
            #print 'alignment', alignment
            #expect_min=10
            for hsp in alignment.hsps:
                #print 'hsp', hsp
                if hsp.expect<expect_max:
                    #print 'hsp.expect<expect_max', hsp
                    if math.fabs(hsp.sbjct_start-hsp.sbjct_end)>100:
                        #print hsp.frame
                        hits_list= [i, hsp.expect, str(alignment.hit_def),'seqbeg',
                                    hsp.sbjct_start , 'seqend', hsp.sbjct_end, 'strand',
                                    hsp.frame[1]]
                        summ_hits.write(str(hits_list))
 #                       print str(hits_list)
                        summ_hits.write('\n')
                        #print '\n'
                        #print 'hits_list', hits_list
                        e_val_list.append(hsp.expect)
                        if hsp.expect<expect_min:
                            expect_min=hsp.expect
                            best_hit=[i, hsp.expect, str(alignment.hit_def), 'seqbeg',
                                      hsp.sbjct_start , 'seqend', hsp.sbjct_end, 'besthit']
                            #print best_hit
    if len(e_val_list) !=0:
#        print best_hit
##        hits_ave = numpy.mean(e_val_list)
##        hits_std = numpy.std(e_val_list)
##        hits_median = numpy.median(e_val_list)
        best_eval = min(e_val_list)
 #       print best_hit
        add_to_table = [best_hit]
        summ_hits.write(str(add_to_table))
 #      print str(add_to_table)
    else:
        summ_hits.write(i+"_no hits")
    summ_hits.write('\n')
summ_hits.close    
    
    
    
    
