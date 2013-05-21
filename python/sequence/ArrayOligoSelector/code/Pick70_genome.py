import string, os,sys
import fasta

#return nonredundant seq dic to with all the ids
#[seq:[id,id...]]
def all_id_seq_dic(id_dic):
    dic ={}
    for id in id_dic.keys():
        seq = string.upper(id_dic[id])
        try:
            dic[seq].append(id)
        except KeyError:
            dic[seq] = [id]
    return dic

def all_id_seq_dic_update( seq_id_dic, newid, newseq):
    foundnewseq =0
    for seq in seq_id_dic.keys():
        if seq == string.upper(newseq):
            seq_id_dic[seq].append(newid)
            foundnewseq =1
            break
    if (not foundnewseq) :
         seq_id_dic[newseq] = [newid]
                                            
if len(sys.argv[:]) !=  4 and __name__ == "__main__":
    print 'USAGE: python Pick70_genome.py genome_fasta input_fasta out_samefile\n'

elif  __name__ == '__main__':
    print "Identifying the input sequences in the genome file ..."

    fdb = open(sys.argv[1],'r')
    finput = open(sys.argv[2],'r')
    out_same = open(sys.argv[3],'w')

    #process input seqfasta file { id:seq,.......}
    id_seq_dic = fasta.fasta_dic(finput)
    seq_id_dic = all_id_seq_dic(id_seq_dic)
    
    #only input ids
    id_dic ={}
    id_same_dic={}
    for id in id_seq_dic.keys():
        id_dic[id] = 0
        id_same_dic[id] =[]

    for seq in seq_id_dic.keys():
        gene_list = seq_id_dic[seq]
        #the values include inputid and genome id
        for id in gene_list:
            try:
                id_dic[id]
                id_same_dic[id].extend(gene_list)
                id_dic[id] = len(seq) #Update Info of the sequence length
            except KeyError:
                continue
            
    #add the genome sequence
    line = fdb.readline()
    seq_list =[]
    id =''
    while (line!=''):
        if line[0] =='>':
            if len(seq_list) == 0 and id =='':
                id = string.strip(line[1:])  #id is whole line
            else:
                for i in range(0,len(seq_list)):
                    seq_list[i] = string.strip(string.replace(seq_list[i],' ',''))
                    seq_list[i] = string.strip(string.replace(seq_list[i],'\n',''))
                    seq_list[i] = string.upper(seq_list[i])
                seq = string.join(seq_list,'')
                #converte to seq:[id,id,...]
                all_id_seq_dic_update( seq_id_dic, id, seq)
                gene_list = seq_id_dic[seq]
                #the values include inputid and genome id
                for id in gene_list:
                    try:
                        id_dic[id]
                        id_same_dic[id].extend(gene_list)
                        id_dic[id] = len(seq) #Update Info of the sequence length
                    except KeyError:
                        continue
                if len(seq_id_dic[seq]) <=1:
                    del seq_id_dic[seq]
                else:
                    seq_id_dic[seq].remove(id)
                id = string.strip(line[1:])  #id is whole line
                seq_list =[]
        else:
            seq_list.append(line)
        line = fdb.readline()
    if len(seq_list)==0 and id =='':
        pass
    else:
        seq = string.join(seq_list,'')
        seq = string.upper(string.strip(string.replace(seq,' ','')))
        seq = string.upper(string.strip(string.replace(seq,'\n','')))
        #converte to seq:[id,id,...]
        all_id_seq_dic_update( seq_id_dic, id, seq)
        gene_list = seq_id_dic[seq]
        #the values include inputid and genome id
        for id in gene_list:
            try:
                id_dic[id]
                id_same_dic[id].extend(gene_list)
                id_dic[id] = len(seq) #Update Info of the sequence length
            except KeyError:
                continue
        if len(seq_id_dic[seq]) <=1:
            del seq_id_dic[seq]
        else:
            seq_id_dic[seq].remove(id)
    fdb.close()




    #output id_same_dic

#here is the groupfile format 
# [id,...]
# id
# strand hit_start hit_end hit_start hit_end .....
# id
# strand hit_start hit_end hit_start hit_end .....
# ....
#
    for id in id_same_dic.keys():  
        out_same.write('[%s]\n' % (id))
        seq_len = id_dic[id]
        genes = id_same_dic[id]
        templist =[]
        for gene in genes:  #the first one is the query
            if gene not in templist:
                out_same.write("%s\n" %(gene))
                out_same.write("%s %s %s %s %s\n" %('+', '1', str(seq_len), '1', str(seq_len)))
                templist.append(gene)
        out_same.write('\n')
    out_same.close()

    print "Identifying the input sequences in the genome file ... DONE"
    sys.exit(0)






