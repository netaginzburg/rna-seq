

nt_file = open("/Users/lab/Documents/Neta_Livne/SEQFILE_upreg_at_10P_Sig_10A_10P.fasta", "r")
nt_source = nt_file.readlines()
nt_file.close()
header_list= []

print "go over nt fasta file"
nts_dict = {}
for line in nt_source:
    if ">" in line:
        print nt_source.index(line)
        header_list.append(nt_source.index(line))
nt_all_seqs =[]

print "go over nt seqs"
for header in header_list:
    print header_list.index(header)
    if header_list.index(header)==len(header_list)-1:
        nt_seq1 = nt_source[header+1:]
    else:
        nt_seq1 = nt_source[header+1:header_list[header_list.index(header)+1]]
    s1= ""
    for item in nt_seq1:
        s2 = item.replace('\n', '')
        s1 = s1+s2
    s4 = nt_source[header]
    s5 = s4[s4.find("|X")+1:s4.find("|len")]
    nts_dict[s5] = s1
print nts_dict
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from Bio import SeqIO

trans_raw= []

print "translating..."
for key, value in nts_dict.items():
    print key
    pos_strand = Seq(value)
    neg_strand = pos_strand.reverse_complement()
    
    trans_1 = pos_strand.translate()
    rec1= SeqRecord(trans_1, id= key+"_trans_1", description="_")
    trans_raw.append(rec1)

    trans_2 = pos_strand[1:].translate()
    rec2= SeqRecord(trans_2, id= key+"_trans_2", description="_")
    trans_raw.append(rec2)
    
    trans_3 = pos_strand[2:].translate()
    rec3= SeqRecord(trans_3, id= key+"_trans_3", description="_")
    trans_raw.append(rec3)
    
    trans_4 = neg_strand.translate()
    rec4= SeqRecord(trans_4, id= key+"_trans_minus_1", description="_")
    trans_raw.append(rec4)
    
    trans_5 = neg_strand[1:].translate()
    rec5= SeqRecord(trans_5, id= key+"_trans_minus_2", description="_")
    trans_raw.append(rec5)
    
    trans_6 = neg_strand[2:].translate()
    rec6= SeqRecord(trans_6, id= key+"_trans_minus_3", description="_")
    trans_raw.append(rec6)

print "saving into fasta file"
SeqIO.write(trans_raw, "/Users/lab/Documents/Neta_Livne/tf_domain_analysis/peptide_seq_upreg_at_10P_Sig_10A_10P.fasta", "fasta")

print "aaaaand we're done :-)"
