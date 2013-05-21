import string, fasta
import os
#GOAL: return the longest length of subsequece that only have the symbols in symbol list with a
#      tolerant level of non-symbol sequence.
#      * tolerant level : [0,1] interval float number

MASKAT = 15

def mask(sequence, symbol_list, TOLERANT):
    if symbol_list == []:
        return 0
    else:
        TOLERANT = float(TOLERANT)
        total_len =0  #current total scanned length
        non_mask_len  =0  #current sequence length DOES not have  symbol
        mask_len  =0  #current BEST sequence length only  have  the symbols with the tolerant level
        for i in range(0, len(sequence)):
            total_len = total_len +1
            if sequence[i] not in symbol_list:
                non_mask_len = non_mask_len +1
            if non_mask_len > total_len* TOLERANT:
                mask_len = max (mask_len, total_len-1)
                total_len =0
                non_mask_len =0
                continue
        mask_len = max (mask_len, total_len)
        return mask_len  
                
#the ad-hoc filter is stretch of AT only region of >= MASKAT bp  in length
def mask_longAT(seq, type): #type "lower", "N"
    retseq = seq

    start = 0
    AT_length=0
    cutoff = MASKAT
    
    mask_list=[]
    for i in range(0,len(retseq)):
        letter = retseq[i]
        if letter not in [ 'G','C', 'g','c']:
            if AT_length ==0:
                start = i
            AT_length= AT_length +1
        else:
            if AT_length >= cutoff:
                mask_list.append([start,i])
            AT_length =0
            
    seq_list=[]
    for i in range(0,len(retseq)):
        seq_list.append(retseq[i])
        
    for start,end in mask_list:
        for i in range(start,end+1):
            if type == 'N':
                seq_list[i]='N' #ad-hoc masking
            if type == "lower":
                seq_list[i]=string.lower(seq_list[i])  #ad-hoc masking
    masked_seq = string.join(seq_list,'')
    return masked_seq

def mask_dust(seq,case):
    fseq = open("~tempdust",'w')
    fseq.write(">temp\n")
    fseq.write(seq+"\n")
    fseq.close()
    os.system('dust ~tempdust > ~tempdustout')
    fseq = open("~tempdustout",'r')
    retseq = string.join(fseq.readlines()[1:],'')
    fseq.close()
    tempseq=[]
    if case == "lower":
        for i in range (0, len(seq)):
            if retseq[i] == 'N':
                tempseq.append(string.lower(seq[i]))
            else:
                tempseq.append(seq[i])
    tempseq = string.join(tempseq,'')
    os.system("rm ~tempdust")
    os.system("rm ~tempdustout")
    return tempseq

def mask_file(infile, outfile,type,case):
    fin = open(infile,'r')
    seq_dic = fasta.fasta_dic(fin)
    index = seq_dic.keys()
    for id in index:
        seq =seq_dic[id]
        if type == "longAT":
            ret_seq = mask_longAT(seq,case)
        if type == "dust":
            ret_seq = mask_dust(seq, case)
        seq_dic[id] = ret_seq
    fout = open(outfile,'w')
    fasta.generate_fasta(seq_dic, fout)
