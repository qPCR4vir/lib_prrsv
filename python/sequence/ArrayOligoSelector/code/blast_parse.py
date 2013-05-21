import string, os,sys,re

#blastp_output no alignment
def blast_parse_blastp (blast_outfile):
    parse_dic = {}
    content = []
    line = " "
    while (line !=""):
        line = blast_outfile.readline()
        if line[:6] == "BLASTP":
            parse(content, parse_dic)
            content =[line]
        else:
            content.append(line)
    parse(content, parse_dic)
    blast_outfile.close()
    return parse_dic  # { Query: [Query_len, [hit_id, hit_score, hit_E, hit_len, hit_identical_bp, hit_cover_bp],[hit_id, hit_score, hit_E, hit_len, hit_identical_bp, hit_cover_bp]],..... }
  # no-hit:  { Query: [Query_len ],..... } 

#blastn output no alignment
# backward compatable to blast_parse
def blast_parse_blastn (blast_outfile):
    return blast_parse(blast_outfile)

###################################################################
#blastn output  with alignment
def blast_parse (blast_outfile): #blast_output_ptr
    ALIGN = 1
    parse_dic = {}
    content = []
    line = " "
    while (line !=""):
        line = blast_outfile.readline()
        if line[:6] == "BLASTN":
            parse_alignment(content, parse_dic,ALIGN)
            content =[line]
        else:
            content.append(line)
    parse_alignment(content, parse_dic, ALIGN)
    blast_outfile.close()
    return parse_dic  # { Query: [Query_len, [hit_id, hit_score, hit_E, hit_len, hit_identical_bp, hit_cover_bp,record],[hit_id, hit_score, hit_E, hit_len, hit_identical_bp, hit_cover_bp,record]],..... }
  # no-hit:  { Query: [Query_len ],..... } 
#####################################################################

#parse one complete output as one string or a list of strings
#no alignment information being included
def parse(content_list, parse_dic):
    ALIGN =0
    return parse_alignment(content_list, parse_dic, ALIGN)


#blastp_output with alignmentg
def  blast_parse_blastp_ali (blast_outfile,ALIGN):
    parse_dic = {}
    content = []
    line = " "
    while (line !=""):
        line = blast_outfile.readline()
        if line[:6] == "BLASTP":
            parse_alignment(content, parse_dic,ALIGN)
            content =[line]
        else:
            content.append(line)
    parse_alignment(content, parse_dic, ALIGN)
    blast_outfile.close()
    return parse_dic  # { Query: [Query_len, [hit_id, hit_score, hit_E, hit_len, hit_identical_bp, hit_cover_bp],[hit_id, hit_score, hit_E, hit_len, hit_identical_bp, hit_cover_bp]],..... }
  # no-hit:  { Query: [Query_len ],..... } 

#blastn output with alignment
# backward compatable to blast_parse
def blast_parse_blastn_ali (blast_outfile,ALIGN):
    return blast_parse_ali(blast_outfile,ALIGN)

#blastn output with alignment
def blastn_coordinate (blast_outfile): #blast_output_ptr
    parse_dic = {}
    content = []
    line = " "
    while (line !=""):
        line = blast_outfile.readline()
        if line[:6] == "BLASTN":
            content = string.join(content,'')
            temp_dic = parse_query_alignment (content) #            parse_alignment(content, parse_dic,ALIGN)
            parse_dic.update(temp_dic)
            content =[line]
        else:
            content.append(line)
    content = string.join(content,'')
    temp_dic = parse_query_alignment (content) 
    parse_dic.update(temp_dic)
    blast_outfile.close()
    return parse_dic  # { Query: [querylen,{hit_id:[strand, pstart, pend, hstart,hend, identity]}]} 

#blastn output with alignment
def blast_parse_ali (blast_outfile,ALIGN): #blast_output_ptr
    parse_dic = {}
    content = []
    line = " "
    while (line !=""):
        line = blast_outfile.readline()
        if line[:6] == "BLASTN":
            parse_alignment(content, parse_dic,ALIGN)
            content =[line]
        else:
            content.append(line)
    parse_alignment(content, parse_dic,ALIGN)
    blast_outfile.close()
    return parse_dic  # { Query: [Query_len, [hit_id, hit_score, hit_E, hit_len, hit_identical_bp, hit_cover_bp],[hit_id, hit_score, hit_E, hit_len, hit_identical_bp, hit_cover_bp], [alignment, alignment,...],..... }
  # no-hit:  { Query: [Query_len ],..... } 

##############################################################################3
#parse one complete output as one string or a list of strings
#alignment information depends
#record is the complete information for the one query such as in a list each entry is for one hit
#gi|9625671|AD169_whole_genome
#            Length = 229354
#            
# Score = 1855 bits (936), Expect = 0.0
# Identities = 936/936 (100%)
# Strand = Plus / Plus
#
#                                                                        
#Query: 1    atgccagccacagacacaaacagcacccacaccacgccgcttcacccagacgcccaacac 60
#            ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#Sbjct: 970  atgccagccacagacacaaacagcacccacaccacgccgcttcacccagacgcccaacac 1029
#
#                                                                        
#Query: 61   acgttacccttacaccacagcaacacacaaccgcatgtccaaacctcggacaaacacgcc 120
#            ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#Sbjct: 1030 acgttacccttacaccacagcaacacacaaccgcatgtccaaacctcggacaaacacgcc 1089
#
#
# Score = 1855 bits (936), Expect = 0.0
# Identities = 936/936 (100%)
# Strand = Plus / Minus
#
#                                                                          
#Query: 1      atgccagccacagacacaaacagcacccacaccacgccgcttcacccagacgcccaacac 60
#              ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#Sbjct: 188497 atgccagccacagacacaaacagcacccacaccacgccgcttcacccagacgcccaacac 188438
#
#                                                                          
#Query: 61     acgttacccttacaccacagcaacacacaaccgcatgtccaaacctcggacaaacacgcc 120
#              ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#Sbjct: 188437 acgttacccttacaccacagcaacacacaaccgcatgtccaaacctcggacaaacacgcc 188378
#
#
#>gi|1167917|gb|U33331.1|HCU33331_Toledo_18KB
#           Length = 18535
#           
# Score = 26.3 bits (13), Expect = 2.8
# Identities = 13/13 (100%)
# Strand = Plus / Minus
#
#                        
#Query: 589 ggccacggccgcc 601
#           |||||||||||||
#Sbjct: 782 ggccacggccgcc 770

def parse_alignment(content_list, parse_dic, ALIGN):
    try:
        seg = content_list + ''
    except TypeError:
        seg = string.join(content_list,'')
    if seg not in ['','\n']:
        part = string.split(seg,'\n\n')
        #Query
        part[2] = string.strip(part[2])
        Query = string.strip(string.join(string.split(string.strip(part[2]),'\n')[0:-1],'') [7:])
        Query_len =  int(float(string.replace(string.replace(string.split(string.strip(part[2]))[-2],'(',''), ',','')))
        parse_dic[Query] = [Query_len]
        Hit_list = []
        if string.find (seg, 'No hits found') == -1:
            record = string.split(string.split(seg,'\n  Database')[0],'\n>')[1:]
            for i in range (0, len(record)):
                data = string.split( string.split(record[i],'Score =')[0],  'Length =')
                hit_id = string.strip(string.join(string.split(data[0],'\n'),''))
                hit_len = int(string.strip(data[1] ))
                Hit_list.append([hit_id, hit_len]) 

        parse_dic[Query].extend (Hit_list)
        if ALIGN and Hit_list:
            parse_dic[Query].append (record)
    return parse_dic


def parse_alignment_list(content_list,parse_dic):
    try:
        seg = content_list + ''
    except TypeError:
        seg = string.join(content_list,'')
    if seg not in ['','\n']:
        part = string.split(seg,'\n\n')
        #Query
        part[2] = string.strip(part[2])
        Query = string.strip(string.join(string.split(string.strip(part[2]),'\n')[0:-1],'') [7:])
#        Query_len =  int(float(string.replace(string.replace(string.split(string.strip(part[2]))[-2],'(',''), ',','')))
        parse_dic[Query] = []
        Hit_list = []

        if string.find (seg, 'No hits found') == -1:
            record = string.split(string.strip(part[6]),'\n')
            for i in range (0, len(record)):
                hitid, hitscore, hitE = string.split(record[i])
                Hit_list.append([hitid, hitscore, hitE]) 

        parse_dic[Query].extend (Hit_list)
    return parse_dic #single entry dic {query_id: [[hitid, hitscore, hitE], [hitid, hitscore, hitE],....]


#####################################################################################
#parse the the entire query alignment to a dictionary as
#{hitid:[],.............}
#the entire query alignment is as with multiple hits and each hit  with multiple segments

#gi|9625671|AD169_whole_genome
#            Length = 229354
#            
# Score = 1855 bits (936), Expect = 0.0
# Identities = 936/936 (100%)
# Strand = Plus / Plus
#
#                                                                        
#Query: 1    atgccagccacagacacaaacagcacccacaccacgccgcttcacccagacgcccaacac 60
#            ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#Sbjct: 970  atgccagccacagacacaaacagcacccacaccacgccgcttcacccagacgcccaacac 1029
#
#                                                                        
#Query: 61   acgttacccttacaccacagcaacacacaaccgcatgtccaaacctcggacaaacacgcc 120
#            ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#Sbjct: 1030 acgttacccttacaccacagcaacacacaaccgcatgtccaaacctcggacaaacacgcc 1089
#
#
# Score = 1855 bits (936), Expect = 0.0
# Identities = 936/936 (100%)
# Strand = Plus / Minus
#
#                                                                          
#Query: 1      atgccagccacagacacaaacagcacccacaccacgccgcttcacccagacgcccaacac 60
#              ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#Sbjct: 188497 atgccagccacagacacaaacagcacccacaccacgccgcttcacccagacgcccaacac 188438
#
#                                                                          
#Query: 61     acgttacccttacaccacagcaacacacaaccgcatgtccaaacctcggacaaacacgcc 120
#              ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#Sbjct: 188437 acgttacccttacaccacagcaacacacaaccgcatgtccaaacctcggacaaacacgcc 188378
#
#
#>gi|1167917|gb|U33331.1|HCU33331_Toledo_18KB
#           Length = 18535
#           
# Score = 26.3 bits (13), Expect = 2.8
# Identities = 13/13 (100%)
# Strand = Plus / Minus
#
#                        
#Query: 589 ggccacggccgcc 601
#           |||||||||||||
#Sbjct: 782 ggccacggccgcc 770

#{ hit_id : [[hit_strand(+/-), qstart, qend, hit_start, hit_end, percentage ], [..],..],.....}
# one query at a time
def parse_query_alignment (alignment):
    hit_alignment_dic = {}
    if  alignment =='':
        return hit_alignment_dic
    
    header = string.split(alignment,'Query= ')[1]
    header = string.split(header,"letters)")[0]
    queryid = string.split(header,'\n')[:-1]
    queryid= string.strip(string.replace(string.join(queryid,''),'\n',''))
    querylen = int(string.strip(string.split(string.replace(header,',',''),'\n')[-1])[1:])
        
    if string.find(alignment,'No hits found')!=-1 :
        return {queryid:{}}
    
    alignment = string.split(alignment,'\n  Database')[0]
    hits = string.split(alignment,'\n>')[1:]
    for hit in hits:
        #to get the id
        header = string.strip(string.split(hit,'Score =')[0])
        lines = string.split(header,'\n')
        hitid = string.join(lines[:-1],'')
        #
        hit_list=[]
        segments = string.split(hit,'Score =')[1:]
        for seg in segments:
            entry_segs  = string.split(string.strip(seg),'\n\n')

            header_seg =string.strip(entry_segs[0])
            start_seg = string.strip(entry_segs[1])
            end_seg = string.strip(entry_segs[-1])
            hit_end = int(string.split(string.split(end_seg,'\n')[-1])[-1])
            hit_start = int(string.split(string.split(start_seg,'\n')[-1])[1])
            query_end = int(string.split(string.split(end_seg,'\n')[0])[-1])
            query_start = int(string.split(string.split(start_seg,'\n')[0])[1])

            if (hit_end > hit_start) and (query_end > query_start):
                strand ='+'
            else:
                strand ='-'
                if (query_start > query_end):
                    temppos = hit_start
                    hit_start = hit_end
                    hit_end = temppos
                    temppos = query_end
                    query_end = query_start
                    query_start = temppos
                
            hit_identity, hit_cover = string.split(string.split(string.split(header_seg,'\n')[1])[2] ,'/')
            hit_identity = string.strip(hit_identity)
            hit_cover= string.strip(hit_cover)
            percentage = 100* float(hit_identity)/float(hit_cover)
            seg_list = [strand, query_start, query_end, hit_start, hit_end, percentage]
            hit_list.append(seg_list)
        hit_alignment_dic [hitid]= hit_list
    return {queryid: [querylen,hit_alignment_dic]}
    # { Query: [querylen, {hit_id:[strand, pstart, pend, hstart,hend, identity]}]} 

###########################################################################################
#parse alignmnet into plusline, minus line, plus_start, minus_start, plus_end, minus_end :
#the alignment has exactly one hit and one segment as the following:

# Score =  140 bits (70), Expect = 3e-37
# Identities = 70/70 (100%)
# Strand = Plus / Plus
#
#                                                                       
#Query: 1   attaataaagaaatgaaaaaccaaaatgaaaatgtacccgaacatgtacaacataatgct 60
#           ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#Sbjct: 676 attaataaagaaatgaaaaaccaaaatgaaaatgtacccgaacatgtacaacataatgct 735
#
#                     
#Query: 61  gaagcaaatg 70
#           ||||||||||
#Sbjct: 736 gaagcaaatg 745

def parse_align(alignment):
    seg = string.split(alignment,'Query:')[1:]
    minusline =''
    plusline =''
    for i in range(0,len(seg)):
        line = string.split(string.strip(seg[i]),'\n')
        if len(line) != 3:
            print seg[i]
            return []
        else:
            plusline = plusline+string.split(line[0])[1]
            minusline = minusline+string.split(line[-1])[2]
            
    plus_start = int( string.split(string.split(string.strip(seg[0]),'\n')[0])[0])
    minus_start = int( string.split(string.split(string.strip(seg[0]),'\n')[-1])[1])
    plus_end = int( string.split(string.split(string.strip(seg[-1]),'\n')[0])[-1])
    minus_end = int( string.split(string.split(string.strip(seg[-1]),'\n')[-1])[-1])

    if plus_start > plus_end: #for parsing blat-out put
        temp = plus_end
        plus_end = plus_start
        plus_start = temp 

        temp = minus_end
        minus_end = minus_start
        minus_start = temp

        plusline = reverse(plusline)
        minusline = reverse(minusline)
    return [plusline,  minusline, plus_start, minus_start, plus_end, minus_end]

compen ={'A':'T', 'T':'A', 'G':'C', 'C':'G', 'U':'A',
	 'M':'K', 'R':'Y', 'W':'W', 'S':'S', 'Y':'R', 'K':'M',
	 'V':'B', 'H':'D', 'D':'H', 'B':'V',
	 'N':'N', 'X':'X',
         'a':'t', 't':'a', 'g':'c', 'c':'g','u':'a',
         'm':'k', 'r':'y', 'w':'w', 's':'s', 'y':'r', 'k':'m',
         'v':'b', 'h':'d', 'd':'h', 'b':'v',
         'n':'n', 'x':'x'}

def reverse (seq):
    from para import compen 
    l = len (seq)
    out_seq = []
    for i in range (0,l):
        out_seq.append(compen[seq[i]])
    out_seq.reverse()
    out_seq = string.join(out_seq,'')
    return out_seq

#output table view of blast result as:
#  >query
#  hit score hit score
#>query
#  hit score hit score

def tableview (parse_dic, hit_cutoff, out):
    keys= parse_dic.keys()
    keys.sort()
    for key in keys:
        Query = key
        out.write('>'+Query+'\n')
        hit_list = parse_dic[key][1:]
        for hit in hit_list:
            try:
                hit_id, hit_score, hit_E, hit_len, hit_identical, hit_cover = hit
            except ValueError:
                print "error:", Query, hit
                break
#LIAM modification
            if hit_E < hit_cutoff:
                out.write(hit_id+' '+str(hit_identical)+' '+ str(hit_E)+' '+ str(hit_cover)+'\t')
#LIAM modification                
                
#            if hit_identical > hit_cutoff:
#                out.write(hit_id+' '+str(hit_identical)+'\t')
        out.write('\n')
    out.close()

def table_to_dic (parse_dic, ftable):
    for seg in string.split(ftable.read(),'>')[1:]:
        data = string.split(seg,'\n')
        id = data[0]
        parse_dic[id]= []
        if len(data) > 1:
            hits = string.split(data[1])
            for i in range (0, len(hits),2):
                parse_dic[id].append([hits[i], int(hits[i+1])])
    ftable.close()















