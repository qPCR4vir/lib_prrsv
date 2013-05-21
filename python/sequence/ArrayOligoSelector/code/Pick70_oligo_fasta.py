import string,sys, os
import nonredundant

def outfasta(id, id_dic, out):
    out.write('>'+id+'\n')
    oligo_seq = id_dic[id]
    out.write(oligo_seq+'\n')

if len(sys.argv[:]) !=3:
    print 'USAGE: python Pick70_oligo_fasta.py oligo_result fout\n'
else:
    fin = open(sys.argv[1],'r')
    out = open(sys.argv[2],'w')

    #{[oligoid:seq],} dictionary
    id_dic ={}
    for seg in string.split(fin.read(),'>'):
        if seg not in ['','\n']:
            id = string.split(seg,'\n')[0]
            for line in string.split(seg,'\n')[1:]:
                if line in ['\n','']:
                    continue
                oligo_seq =  string.split(string.strip(line))[-1]
                oligo_pos = string.split(line)[0]
                oligo_id = id+'_'+str(oligo_pos)
                id_dic[oligo_id]= oligo_seq

    
    dic_id = nonredundant.non_redu_seq_dic(id_dic)
    index = dic_id.keys()
    index.sort()
    dup_id =[]
    for id in index:
        if len(dic_id[id]) != 1:
            list = dic_id[id]
            list.sort()
            if list not in dup_id:
                dup_id.append(list)
        else:
            outfasta(id, id_dic, out)

    for id_list in dup_id:
        outfasta(id_list[0],id_dic,out)

    out.close()
    






