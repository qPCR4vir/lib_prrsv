import string

#groupfile format: id are group member id, as long as the union of all the ids define a set for teh group
# [id]
# id
# strand q_start q_end hit_start hit_end hit_start hit_end .....
# id
# strand q_start q_end hit_start hit_end hit_start hit_end .....
# ....
#
# [id,...]
# id
# strand q_start q_end hit_start hit_end hit_start hit_end .....
# id
# strand q_start q_end hit_start hit_end hit_start hit_end .....
# ....
#
#{gene_id:[[hit_id, hit_strand, qstart, qend, hitstart, hitend],[hit_id, hit_strand, qstart, qend, hit_start, hit_end]]..., ...}

def All_Group_dic(groupfile): # IN parameter file pointer
    Group_dic = {}
    for contig in string.split(string.strip(groupfile.read()), '//\n'):
        contig = string.strip(contig)
        if contig[0] != "[":
            return 0
        if contig != '':
            #each contig
            #process grouping info
            for lines in string.split('\n'+contig,'\n[')[1:]:
                #each group
                gene_list = []
                entries = string.split(string.strip(lines),'\n')
                queryid = entries[0][:-1]
                queryid = string.replace(queryid,' ','')
                queryid = string.replace(queryid,'\t','')
                for i in range (1, len(entries), 2):
                    #each gene
                    gene_id = string.strip(entries[i])
                    gene_id = string.replace(gene_id,' ','')
                    gene_id = string.replace(gene_id,'\t','')
                    strand, qstart, qend, hstart, hend = string.split(string.strip(entries[i+1]))
                    qstart = int(qstart)
                    qend = int(qend)
                    hstart = int(hstart)
                    hend = int(hend)
                    each_gene_info = [gene_id,strand, qstart, qend, hstart, hend]
                    gene_list.append(each_gene_info)  #no-sorting, keep order
                try:
                    Group_dic[queryid].extend(gene_list)
                except KeyError:
                    Group_dic[queryid] = gene_list
                    
    groupfile.close()
    return Group_dic

#def Find_Group (id, groups):
#    for group in groups:
#        if id in group:
#            return group
#    return []  #[]:error

def Generate_Group(fastafile, outfile): #fastafile: ptr or filename
    import string,sys,os
    try:
        fastafile +'a'
        fastafile = open(fastafile,'r')
    except TypeError:
        pass
        
    import fasta
    id_seq_dic = fasta.fasta_dic(fastafile)
    ids = id_seq_dic.keys()
    for id in ids:
        outfile.write('[%s]\n' % (id))
        outfile.write(id+'\n')
        outfile.write('%s %s %s %s %s\n' %('+','1','1','1','1'))
    outfile.close()





