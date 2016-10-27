####### script: seqs_from_entrez1.py #######

from Bio import Entrez, SeqIO

old=open('/Users/lab/Documents/Neta_Livne/TFs_pool/Dmel_TFs_FBids.txt', 'r') #original list from GO
#old=open('/Users/lab/Documents/Neta_Livne/TFs_pool/Dmel_TFs_FBids_new1.txt', 'r') #original list from GO
old_tfs=old.readlines()

FBids = []
for i in xrange(len(old_tfs)): ##create a list of tfs fbid's from the GO list file
    z=old_tfs[i].split('\n')
    FBids.append(z[0])

Entrez.email = "neta.livne88@gmail.com"
files_list = []
for fbid in FBids: #iterate over the list of old_py tfs/whatever
    print fbid
    handle1 = Entrez.esearch(db="protein", term=fbid)
    record1 = Entrez.read(handle1)
    record1["IdList"] #Each of the IDs in IdList is a GenBank identifier. ###each ID is an isoform!!!
    #See section 9.6 for information on how to actually download these GenBank records.
    for gbid in record1["IdList"]: #get the information from genebank about the isoform
        reslist=[]
        handle2 = Entrez.efetch(db="protein", id=gbid, rettype="gp", retmode="xml")
        record2 = Entrez.read(handle2)
        handle2.close()
        record2[0]["GBSeq_definition"]
        prot_seq = record2[0]['GBSeq_sequence']
        path='/Users/lab/Documents/Neta_Livne/tfs_entrez_seqs/'
        file_name = path + fbid +'_GBprotein_' + gbid #add fb transcript id?
        seq_file=open(file_name, 'w')
        seq_file.write(prot_seq)
        seq_file.close()
        f_short_name = fbid +'_GBprotein_' + gbid
        files_list.append(f_short_name)
#        record2[0]['GBQualifier_value']
#        record2[0]["GBSeq_origin"]

print "\n", "files_list", files_list
