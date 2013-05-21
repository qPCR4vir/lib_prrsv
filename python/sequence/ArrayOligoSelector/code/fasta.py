import string, sys,os

# IN parameter:
#      fastafile: file pointer
# RETURN:
#      [[id,seq],...]

def fasta(fastafile):
    id_seq_list =[]
    for record in string.split('\n'+string.strip(fastafile.read()) ,'\n>'):
        if record != '':
            data = string.split(record,'\n')
            id = data[0] #id is line'>"
            seq = string.strip(string.join(data[1:],''))
            id_seq_list.append([id,seq])
    fastafile.close()
    return id_seq_list

def fasta_dic(fastafile):
    id_seq_dic ={}

    line = fastafile.readline()
    seq_list =[]
    id =''
    while (string.strip(line) ==""):
        if line =="":
            return id_seq_dic
        line = fastafile.readline()
        
    while (line!=''):
        if line[0] =='>':
            if len(seq_list) == 0 and id =='':
                id = string.strip(line[1:])  #id is whole line
            else:
                for i in range(0,len(seq_list)):
                    seq_list[i] = string.strip(string.replace(seq_list[i],' ',''))
                    seq_list[i] = string.strip(string.replace(seq_list[i],'\n',''))
#                    seq_list[i] = string.upper(seq_list[i])
                seq = string.join(seq_list,'')
                id_seq_dic[id] = seq
                id = string.strip(line[1:])  #id is whole line
                seq_list =[]
        else:
            seq_list.append(line)
        line = fastafile.readline()

    if len(seq_list)==0 and id =='':
        pass
    else:
        seq = string.join(seq_list,'')
        seq = string.strip(string.replace(seq,' ',''))
        seq = string.strip(string.replace(seq,'\n',''))
#        seq = string.upper(string.strip(string.replace(seq,' ','')))
#        seq = string.upper(string.strip(string.replace(seq,'\n','')))
        id_seq_dic[id] = seq

    fastafile.close()
    return id_seq_dic


def generate_fasta (fasta_dic, fout):  #fout : ptr to output file
    for key in fasta_dic.keys():
        fout.write('>'+key+'\n')
        fout.write(fasta_dic[key]+'\n')
    fout.close()
    return

#all files are in fasta format
def find_fasta(fseq, fid, fout):
    id_dic = fasta_dic(fid)
    seq_dic = fasta_dic(fseq)

    temp_dic ={}
    keys = id_dic.keys()
    for id in keys:
        try:
            temp_dic[id] = seq_dic[id]
        except KeyError:
            print "no seq in fseq file"

    generate_fasta(temp_dic, fout)

#check if file fseq is in the correct fasta format
#return 0 as fail
#return 1 as success
def check_format_fasta(fseq, seq_size): #seq_size is the dictionary of id:length of squence
    WRONGFORMAT = 0
    SUCCESS = 1
    WRONGSIZE= 2
    IDWITHSPACE =3
    FAIL = 0

    id =""
    while (1):
        line = fseq.readline()
        if line == '':
            return WRONGFORMAT
        elif string.strip(line) == '':
            continue
        elif string.strip(line[0]) != '>':
            return FAIL
        else:
            id = string.strip(line)[1:]
            if string.find(id, ' ')!=-1:
                return IDWITHSPACE
            break

    idline = 0
    sequence= 0
    while (1):
        line = fseq.readline()
        if line == '':
            break
        if (idline ==0):
            if string.strip(line) == '':
#                if sequence != 0:
#                    idline = not (idline)
#                    try:
#                        if seq_size[id] != sequence:
#                            return WRONGSIZE
#                    except KeyError:
#                        seq_size[id] = sequence
#                        id = string.strip(line)[1:]  #
#                        sequence =0 #
                continue
            elif string.find(line, '>') != -1:
                idline = not(idline)
                try:
                    if seq_size[id] != sequence:
                        return WRONGSIZE
                except KeyError:
                    seq_size[id] = sequence
                id = string.strip(line)[1:]
                if string.find(id, ' ')!=-1:
                    return IDWITHSPACE
                sequence = 0
            else:
                sequence = sequence + len(string.strip(line))
        if (idline == 1):
            if string.strip(line) == '':
                continue
            elif string.find(line, '>') != 0:
                return FAIL
            else:
                idline = not(idline)
    fseq.close()

    if ( idline !=1 and sequence != 0 ) or (idline ==1 and sequence ==0):
        return SUCCESS

    return FAIL

#converting all sequences into UPPER case
def ToUpper(seq_dic):
    for id in seq_dic.keys():
        seq = seq_dic[id]
        seq_dic[id] = string.upper(seq)

