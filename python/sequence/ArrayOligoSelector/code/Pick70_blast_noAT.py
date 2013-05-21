import string, os, sys, re
import fasta
import Pick70_mask

TRUE =1
FALSE =0
ADD =0
SUB =1
WORD = 10
MAX = 50
E_VALUE = 10
FLAGED  = 'F'
MASKALL = 0.5

#energy_dic = {pos:[min_energy,[hit_id, energy],[hit_id, energy],...], ....}
def Blast(SEQ, id, Group_dic_same, DB,  OLIGOLEN,  TRACEFLAG, STRAND): 
    seqfile = Write_SeqFile(SEQ,id)
    blast_output = Excute_Blast(seqfile, DB,  TRACEFLAG, STRAND)
    if blast_output != 0: #success
        energy_dic = Parse_Blast (id, blast_output, Group_dic_same,  OLIGOLEN)
        Mask_energy_dic(energy_dic, seqfile, len(SEQ), OLIGOLEN)
    else:
        energy_dic =0 #fail
    Del_Temp_Seqfile (seqfile)
    return energy_dic

def Del_Temp_Seqfile (seqfile):
    os.environ['file'] = seqfile
    os.system('rm "$file"')
    return

#write sequencefile for one sequence and perform the ad-hoc filtering to compensate for the dust program because I have noticed that dust filter will miss long at rich region
#the ad-hoc filter is stretch of AT only region of >=MASKAT  in length
def  Write_SeqFile(seq,id) :
    masked_seq= Pick70_mask.mask_longAT(seq,"lower")
    outfile = '~temp' + string.replace(id[:100],'/','')
    out = open(outfile, 'w')
    out.write('>'+id+'\n')
    out.write(masked_seq +'\n') 
    out.close()
    return outfile

def Mask_energy_dic(energy_dic, seqfile, seqlen, OLIGOLEN):
    os.environ['file'] = seqfile
    f = os.popen('./code/dust "$file"')  #mask with dust
    result = f.read()
    f.close()
    masked_seq = string.join(string.split(result,'\n')[1:],'')

    temp_seq=[]
    for i in range (0,len(masked_seq)):
        if masked_seq[i] == string.lower(masked_seq[i]):
            temp_seq.append('N')
        else:
            temp_seq.append(masked_seq[i])
    masked_seq = string.join(temp_seq,'')
    
    cutoff = OLIGOLEN* MASKALL
    for i in range(0,seqlen - OLIGOLEN+1):
        oligoseq = masked_seq[i:i+OLIGOLEN+1]
        if string.count(oligoseq,'N') >= cutoff:
            energy_dic[i] = [FLAGED]

def Excute_Blast(seqfile, DB,  TRACEFLAG, STRAND):
    os.environ['DB']= DB
    os.environ['seqfile'] = seqfile
    os.environ['E_value'] = str(E_VALUE)
    os.environ['STRAND'] = STRAND
    os.environ['MAX'] = str(MAX)
    fr,fw,fe = os.popen3('./code/blastall -p blastn -d $DB -i "$seqfile" -e $E_value -F "m D" -U -v $MAX -b $MAX -S $STRAND -a $cpunum -G 5 -E 0 -o "~out.blastout"')
    try:
        os.wait()
    except OSError:
        pass
    
    errormessage = fe.read()
    fr.close()
    fw.close()
    fe.close()
    if errormessage != "":
        print "Fail on blastn, program terminated."
        print "Error message from blastn:", errormessage
        #sys.exit(1)
        return 0
    
    fin = open("~out.blastout", 'r')
    result = fin.read()
    fin.close()

    if TRACEFLAG:
        print result
    return result


#the same alignment has been covered
def  WithInRange(coverArea, p_start, p_end):
    for s,e in coverArea:
        if s<= p_start and e>= p_end:
            return TRUE
    coverArea.insert(0,[p_start,p_end])
    return FALSE

#def Parse_Blast (output, Group, seqlen, OLIGOLEN, start, end, SCAN_INTERVAL):
#energy_dic = {pos:[min_energy,[hit_id, energy],[hit_id, energy],...], ....}
#Group_dic_same is as [[hit_id, hit_strand, qstart, qend,hstart,hend,...],[hit_id, hit_strand, qstart,qend, hstart, hend,...]]
def Parse_Blast (id , output, Group_dic_same,  OLIGOLEN):
    import blast_parse
    parse_dic = {}
    blast_parse.parse_alignment(output, parse_dic, 1)
    #now we need the information to check that not only that we hit the same id but in the corresponding region too
    new_id = id[:]
    new_id = string.replace(new_id,' ','')
    new_id = string.replace(new_id,'\t','')

    try:
        Group = Group_dic_same[new_id] # [[hit_id, hit_strand, qstart, qend,hstart,hend,...],[hit_id, hit_strand, qstart,qend, hstart, hend,...]]
        print Group
    except KeyError:
        Group = []
    GroupIds=[] # the ids
    for groupHit in Group:
        hitid = groupHit[0]
        GroupIds.append(hitid)
    
    ret = parse_dic.values()[0] #only one entry there for the entire result
    Query_len = ret[0]
    gene_end = Query_len - OLIGOLEN 
    hit_list = ret[1:-1]
    alignment = ret[-1]

    energy_dic={}
    empty_dic ={}
    for pos in range (0, gene_end+1):
        empty_dic[pos] =[ 0,-1] #energy, location

    if len(hit_list) == 0:
        return {}
    else:
        hitGroup =[]
        for j in range (0, len(hit_list)):
            hit_id = hit_list[j][0] #hit_id , hit_len 
            hit_id = string.replace(hit_id,' ','')
            hit_id = string.replace(hit_id,'\t','')
            
            aligns = string.split(alignment [j], ' Score =')[1:]
            coverArea = []
            temp_energy_dic = empty_dic.copy()
            for i in range(0, len(aligns)): #individula alignment of individula hit
                plusline, minusline, p_start, m_start, p_end, m_end  =  blast_parse.parse_align(aligns[i])

                #we need to check if the hit can be ignored because it is the query itself
                #according to the Group_Same_Dic
                ignore =0
                if hit_id in GroupIds: # found id , possible
                    for hitInfo  in Group: #hitInfo is [hitid, hitstrand(+/-), qstart, qend, hitstart, hitend]
                        HITid, HITstrand, Qstart, Qend, Hstart, Hend = hitInfo
                        #check if id match
                        if HITid != hit_id:
                            continue

                        #check if strand information match
                        tolerant_diff = 5
                        if m_start < m_end: #+ strand in  the alignment
                            if HITstrand =='-':
                                continue
                        else:
                            if HITstrand =='+':
                                continue
                        # a key point is if the hit is itself it must be detected as a complete unit, which is almost the same as the hit in the GroupDicSame # + strand in the alignment
                        #if abs(Qstart - p_start)<=5  and abs(Qend- p_end)<=5 and abs(Hstart-m_start)<=5 and abs(Hend -m_end)<=5:
                        # Qstart, Qend : in group file
                        # p_end,hit_id: in blast
                        

                        if ((Qstart <=  p_start)  and ( Qend >=  p_end ))  :
                            if (( HITstrand =='+') and (Hstart <= m_start) and (Hend >= m_end)):
                                ignore = 1
                            elif (( HITstrand =='-') and (Hstart >= m_start) and (Hend <= m_end)):
                                ignore =1
                            if ignore:
                                print "ignored segment: query:",p_start,p_end,hit_id,":", Hstart,Hend
                                break

                        if abs(Qstart - p_start)<=5  and abs(Qend- p_end)<=5 and abs(Hstart-m_start)<=5 and abs(Hend -m_end)<=5:
                            ignore = 1
                            print "ignored segment: query:",p_start,p_end,hit_id,":", Hstart,Hend
                            break

                if ignore :
                    continue

                if WithInRange(coverArea, p_start, p_end):
                    continue
                a_length = len(plusline)
                
                p_start = p_start -1
                p_end = p_end -1
                if (m_start< m_end):
                    m_start = m_start-1
                    m_end = m_end  - 1
                else:
                    m_start = m_start +1
                    m_end = m_end +1
                
                energy_list = compute_energy2 (plusline, minusline)
                #start energy
                start_energy = 3.4
                for i in range (0, min(OLIGOLEN, a_length)):
                    start_energy = start_energy + energy_list[i][ADD]

                #position = p_start
                if (p_start <= gene_end):
                    if start_energy < temp_energy_dic[p_start][0]:
                        temp_energy_dic[p_start] = [ start_energy, m_start]

                poffset = 0 #check gap
                moffset = 0
                #the alignment of after the alignment
                energy = start_energy
                for pos in range(1,a_length - WORD):
                    if pos != 0 and plusline[pos-1] == '-': #check gap
                        poffset = poffset - 1
                    if pos!=0 and minusline[pos-1] =='-':
                        moffset = moffset - 1
                        
                    position = p_start+pos
                    preal_position = position + poffset #compensate for gap

                    if (m_start< m_end):
                        position = m_start+pos
                        mreal_position = position +moffset
                    else:
                        position = m_start -pos
                        mreal_position = position - moffset
                    
                    if (preal_position > gene_end) :#or (real_position < 0): #check bound
                        break

                    else:
                        energy = energy + energy_list[pos-1][SUB]
                        end = pos+OLIGOLEN -1
                        if end < a_length:
                            energy = energy + energy_list[end][ADD]

                        if energy < temp_energy_dic[preal_position][0]:
                            temp_energy_dic[preal_position] = [energy,mreal_position]

                poffset = 0  #no gap
                moffset =0
                end = min(OLIGOLEN, a_length)

                energy = start_energy
                #the alignemnet is completely inside of the oligo selection
                if (end == a_length):
                    start = p_start+end-OLIGOLEN  
                    for i in range(max(0,-(start) ), OLIGOLEN-end):
                        preal_position = start + i
                        if preal_position >gene_end:
                            break

                        if energy < temp_energy_dic[preal_position][0]:
                            if m_start < m_end:
                                m_location = m_start + preal_position - p_start
                            else:
                                m_location = m_start - (preal_position - p_start)
                            temp_energy_dic[preal_position] = [energy,m_location]

                #before the alignment 
                for pos in range(end-1, max(WORD, OLIGOLEN - p_start)-1 , -1): #version 2.7 correction
                    preal_position = p_start+ pos -OLIGOLEN
                    if preal_position > gene_end:
                        break
                    energy = energy - energy_list[pos][ADD]
                    if energy < temp_energy_dic[preal_position][0]:
                        if m_start < m_end:
                            m_location = m_start + preal_position - p_start
                        else:
                            m_location = m_start - (preal_position - p_start)
                        temp_energy_dic[preal_position] = [energy,m_location]
            
            for key in temp_energy_dic.keys():
                if temp_energy_dic[key] != [0,-1]:
                    try:
                        energy_dic[key][0] = min(energy_dic[key][0], temp_energy_dic[key][0])
                    except KeyError:
                        energy_dic[key] = [temp_energy_dic[key][0]]
                    energy_dic[key].append([hit_id, temp_energy_dic[key][0], temp_energy_dic[key][1]])                   
    return energy_dic #{pos:[min_energy,[hit_id, energy],[hit_id, energy],...], ....}

def compute_energy2(p_line, m_line): 
    os.environ['plusline']= p_line
    os.environ['minusline']= m_line
    f = os.popen('./code/energy \"$plusline\" \"$minusline\" 1')
    retData =f.read()
    f.close()
    energy_list = []
    for line in string.split(retData,'\n')[:-1]:
        data = string.split(line)
        add = float(data[0])
        sub = float(data[1])
        energy_list.append([add, sub])
    return energy_list

def compute_energy(p_line, m_line): 
    os.environ['plusline']= p_line
    os.environ['minusline']= m_line
    f= os.popen('./code/energy \"$plusline\" \"$minusline\" 1')
    retData = f.read()
    f.close()
    energy =  float(string.split(string.strip(retData),'\n') [-1])
    return energy





