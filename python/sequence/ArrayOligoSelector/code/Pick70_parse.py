import string, os, sys
import Pick70_mask

FALSE = 0
TRUE = 1
SW_INCREMENT = 10
REPEAT_INCREMENT = 1
GC_INCREMENT = 1

#format for data_dic : data_dic[pos] = [bl, gc, rp, sw, seq, hit_info] 
def filter_unqiue (data_dic, TOP_FIVE, BEST, USER_UNIQUE):
    print BEST, TOP_FIVE
    new_dic = {}
    for key in data_dic.keys():
        score = data_dic[key][0] #eg. -17.34
        #if ( (score < USER_UNIQUE) and ( ( score < TOP_FIVE) or  ( score < BEST- 5 ))) :  #the oligo is not satisfy uniqueness criteria  # OR gate
        if ( (score < USER_UNIQUE) or ( ( score < TOP_FIVE) or  ( score < BEST- 5 ))) :  #the oligo is not satisfy uniqueness criteria # AND gate
            continue
        else:
            new_dic[key] = data_dic[key][:]
    print 'uniqueness:',len(data_dic.keys()), len(new_dic.keys())
    return new_dic

def filter_SW (data_dic, SW):
    new_dic = {}
    for key in data_dic.keys():
        score = data_dic[key][3]
        if score <= SW:
            new_dic[key] = data_dic[key][:]
    print 'self-annealing:',len(data_dic.keys()), len(new_dic.keys())
    return new_dic


def filter_REPEAT (data_dic, REPEAT):
    new_dic = {}
    for key in data_dic.keys():
        score = data_dic[key][2]
        if score <= REPEAT:
            new_dic[key] = data_dic[key][:]
    print 'internal repeat:',len(data_dic.keys()), len(new_dic.keys())
    return new_dic

def filter_mask (data_dic, SYMBOL_LIST, TOLERANT, MASK):
    new_dic = {}
    for key in data_dic.keys():
        seq = data_dic[key][4]
        os.environ['seq'] = seq
        mask = Pick70_mask.mask(seq,SYMBOL_LIST, TOLERANT)
        if mask <= MASK:
            new_dic[key] = data_dic[key][:]
    print 'mask:',len(data_dic.keys()), len(new_dic.keys())
    return new_dic

def filter_GC (data_dic, GC, GC_extend):
    new_dic = {}
    if len(data_dic.keys()) ==0 :
        return new_dic

    for key in data_dic.keys():
        score = float(data_dic[key][1])
        if (score <= (GC +GC_extend)) and ( score >= (GC -GC_extend)):
            new_dic[key] = data_dic[key][:]

    print 'GC',str(GC-GC_extend),str(GC+GC_extend),':',len(data_dic.keys()), len(new_dic.keys())
    return new_dic

#def filter_pos (data_dic, pos_threshold, seq_len):
def filter_pos (data_dic, candidate_list ): #pos_threshold):
    index = data_dic.keys()[:]
    for pos in index:
        candidate_list.append(pos)  
    candidate_list.sort() 
    candidate_list.reverse()  #3' perference
    return candidate_list
    #new_dic = {}
    #if len(data_dic.keys()) ==0 :
    #    return new_dic

    #for key in data_dic.keys():
    #    score = seq_len - float(key) #distance from end of gene
    #    if score <= pos_threshold:
    #        new_dic[key] = data_dic[key][:]

    #print "Distance To 3' end", pos_threshold,":", len(data_dic.keys()), len(new_dic.keys())
    #return new_dic

def find_unique(data_dic):
    index = data_dic.keys()[:]
    list =[]
    for pos in index:
        bl = data_dic[pos][0]  #eg. -15.0
        list.append(bl)
    list.sort()
    list.reverse() #[-15,-16,....]
    top_five_pos = len(list)/20
    return [list[top_five_pos], list[0]] #top five, best uniqueness score
    
def find_SW(data_dic):
    index = data_dic.keys()[:]
    list =[]
    for pos in index:
        sw = data_dic[pos][3] 
        list.append(sw)
    list.sort()
    return [ list[len(list)/3], list[-1]] #top 33% and the highest
    
def find_REPEAT(data_dic):
    index = data_dic.keys()[:]
    list =[]
    for pos in index:
        rp = data_dic[pos][2]
        list.append(rp)
    list.sort()
    return [ list[len(list)/3], list[-1]] #top 33% and the highest

def oligo_dup(index,data_dic,geneid, out):
    for i in range(0, len(index)):
        key = index[i]
        gc = data_dic[key][1]
        rp =data_dic[key][2]
        sw =data_dic[key][3]
        seq =data_dic[key][4]
        hit_info= data_dic[key][5]
        hit_info = parse_hit_info(hit_info)
        out.write('>'+geneid+'_'+str(key)+'\t'+str(gc)+'\t'+str(rp)+'\t'+str(sw)+'\n')
        os.environ['seq'] = seq
        energy =  float(string.split(os.popen('./code/energy \"$seq\" \"$seq\" 1').read(),'\n')[-1])
        hit_info.insert(0,[geneid,energy,key])

        for i in range(0, len(hit_info)):
            id , energy, pos = hit_info[i]
            out.write(id+'\t'+str(energy)+'\t'+str(pos)+'\t')
        out.write('\n')
            
def Oligo_fasta (index, data_dic,id, out):
    for i in range(0, len(index)):
        key = index[i]
        out.write('>%s_%s\n' % (id, key))
        out.write(data_dic[key][4])
        out.write("\n")
 
def process(info, no_design_list,SYMBOL_LIST, TOLERANT, MASK, USER_UNIQUE, NUMBER_OLIGO, oligo_dup_file,out):
    if info == '':
        return
    else:
        id , data_dic ,seq_len= parse_entry (info)
        if data_dic == {}:
            no_design_list.append(id)
            del data_dic
            return
        
        oligo_pos_list =[] #final selected list

        #filter for uniqueness
        TOP_FIVE, BEST = find_unique(data_dic)
        candidate_dic = filter_unqiue (data_dic, TOP_FIVE, BEST, USER_UNIQUE)

        #filter for mask sequences
        if (MASK >= 0):
            candidate_dic = filter_mask (candidate_dic, SYMBOL_LIST, TOLERANT, MASK)

        #if no oligos in the new_dic set, we have to only use the uniqueness filter
        if len(candidate_dic) == 0:
            candidate_dic = filter_unqiue (data_dic, TOP_FIVE, BEST,USER_UNIQUE)

        #self-binding and internal repeat relaxation procedure
        SW_START, SW_limit = find_SW (data_dic)
        REPEAT_START, REPEAT_limit = find_REPEAT(data_dic)
        SW = SW_START -SW_INCREMENT
        REPEAT = REPEAT_START -REPEAT_INCREMENT
        pass_dic ={}
        anchor_pos = -1
        
        while len(oligo_pos_list) < NUMBER_OLIGO:
            if SW < SW_limit:
                SW = SW+SW_INCREMENT
            if REPEAT < REPEAT_limit:
                REPEAT = REPEAT +REPEAT_INCREMENT
            print SW, REPEAT,SW_limit, REPEAT_limit
            pass_dic = filter_SW (candidate_dic, SW)
            pass_dic = filter_REPEAT(pass_dic, REPEAT)

            #find oligos from pass_dic:
            oligo_pos_list = find_oligos (anchor_pos, pass_dic, NUMBER_OLIGO, GC, GC_extention, seq_len)
            if oligo_pos_list:
                anchor_pos = oligo_pos_list[0]
            
            if (SW >= SW_limit ) and (REPEAT >= REPEAT_limit):
                break
        
        if oligo_pos_list:
            add_collection( oligo_pos_list, data_dic) 
            Oligo_fasta (oligo_pos_list, data_dic,id, out)
            oligo_dup(oligo_pos_list,data_dic,id, oligo_dup_file)
        else:
            no_design_list.append(id)
        return

def add_collection (oligo_list, data_dic):
    for pos in oligo_list:
        print '*',pos, data_dic[pos][:-1],
        for entry in data_dic[pos][-1]:
            print entry,
        print
    return

def no_homolog_pass(data_dic, SW, REPEAT,SYMBOL_LIST, TOLERANT):
    new_dic = filter_SW (data_dic,SW)
    new_dic = filter_REPEAT (new_dic,REPEAT)
    if SYMBOL_LIST != []:
        new_dic = filter_mask (new_dic, SYMBOL_LIST, TOLERANT, MASK)
    return new_dic

def filter_pass(data_dic,SYMBOL_LIST, TOLERANT):
    if SYMBOL_LIST != []:
        new_dic = filter_mask (data_dic, SYMBOL_LIST, TOLERANT, MASK)
    else:
        new_dic = data_dic.copy()
            
    new_list =[]
    for pos in new_dic.keys():
        bl, gc, rp, sw, seq, hit_info = new_dic[pos] 
        new_list.append([-bl,rp,-pos,gc])  # -pos is in the list to be sorted from 3' to 5' #energy reversed
    new_list.sort()
    if len(new_list) == 0:
        return []
    else:
        ret =  [ - new_list[0][2]]
        return ret

def  parse_hit_info(data_list):

    list =[]
    for i in range (0,len(data_list),3):
	id = data_list[i]
        value = data_list[i+1]
        pos =  data_list[i+2]
        list.append([id,value,pos])
    return list

def parse_entry (information):
    if information == '':
        return []
    #each entry
    data_dic ={}
    lines = string.split(information,'\n')
    title_line = lines[0]
    id,seq_len =string.split( string.strip(title_line),'\t')
    seq_len= int(seq_len)
    print '>'+id, seq_len
    for line in lines[1:]:
        if line not in ['','\n']:
            data =string.split(line,'\t')
            gc = data[2]
            if gc == 'F':
                continue
            else:
                gc = float(gc)
            rp = int(data[3])
            sw = int(data[4])  
            bl = data[1]
            if bl =='F':
                continue
            else:
                bl = float(bl)
            pos = int(data[0])
            seq = data[5]
            hit_info = data[6:]
            data_dic[pos] = [bl, gc, rp, sw, seq, hit_info] #,hit_list]
    return [id, data_dic,seq_len]

#anchor_pos is used to ensure the optimal pick for unique, self-binding and internal repeat and gc
#anchor_pos is always the first in the pos_list
def find_oligos (anchor_pos , candidate_dic, NUMBER_OLIGO, GC, GC_extention,  seq_len):
    if candidate_dic == {}:
        return []
    
    pos_list =[]
    pass_dic={}
    interval = min (OLIGO_LEN*0.7, seq_len/ float( NUMBER_OLIGO))
    not_done = 1
    GC_extention = GC_extention - GC_INCREMENT
#    pos_threshold = int(seq_len*0.2)
#    pos_increment = int(seq_len*0.1)
    
    while not_done:
        candidate_list =[]
        GC_extention = GC_extention + GC_INCREMENT
#        pos_threshold = pos_threshold + pos_increment
        pass_dic = filter_GC (candidate_dic,GC, GC_extention)
        
        if len(pass_dic)!=0:         #found at least one oligo in the GC range
            candidate_list = filter_pos (pass_dic, candidate_list)

            #index = pass_dic.keys()[:]
            #for pos in index:
            #    bl, gc, rp,sw,seq, hit_info = pass_dic[pos] 
            #    candidate_list.append(pos)  
            #candidate_list.sort() 
            #candidate_list.reverse()  #3' perference
                
        else: #no oligo fits increase GC
            continue

        if anchor_pos == -1:
            pos_list = [ candidate_list[0]]  #pos being reverted,
        else:
            pos_list = [anchor_pos]
        count = 1
        
        if  count < NUMBER_OLIGO: #if more than one oligo desired
            while 1:
                print "interval:", interval
                for j in range (1, len(candidate_list)):
                    pos = candidate_list[j]
                    if match_interval (pos, pos_list, interval):
                        pos_list.append(pos)
                        count = count +1
                        if count >= NUMBER_OLIGO:
                            not_done = 0
                            return pos_list
                if interval >=10:
                    interval = interval -5
                else:
                    break
        else:
            return pos_list
        
        if len(candidate_list) == len(candidate_dic):
            return pos_list

def match_interval (pos, pos_list, interval):
    for p in pos_list:
        if abs(p - pos) <  interval:
            return FALSE
    return TRUE
    

if len(sys.argv[:]) not in [7, 8, 11] and __name__ == '__main__':
    print 'USAGE: python Pick79_parse.py filein fileout foligo_dup GC_percent(eg: 35.5)   OLIGO_LEN(bp)  HowManyOligos [optional binding energy cutoff] [mask_length(interger)  mask_symbol(eg, AT)  mask_tolerance(any number between 0 and 1, 0 minimium tolerance)]\n'

elif __name__ == '__main__':
    filename  = sys.argv[1]
    filecount =0
    
    out = open(sys.argv[2],'w')
    oligo_dup_file = open(sys.argv[3],'w')

    GC = float(sys.argv[4])
    OLIGO_LEN = int(sys.argv[5])
    NUMBER_OLIGO = int(sys.argv[6])
    GC_extention = 0.0
    SYMBOL_LIST =[]
    TOLERANT =0
    MASK = -1

    if len(sys.argv[:]) ==8:
        USER_UNIQUE = float(sys.argv[7])
    else:
#        USER_UNIQUE = 0 #OR gate
        USER_UNIQUE = -1E308 #AND gate
        
    if len(sys.argv[:]) ==11:
        SYMBOL_LIST = list(sys.argv[9])
        TOLERANT = float(sys.argv[10])
        MASK = int(sys.argv[8])

    no_design_list =[]

    currentfile = filename+str(filecount)
    while (os.path.exists(currentfile)):
        fin = open(currentfile,'r')
        fin.readline()
        line = fin.readline()
        entry = []
        while (line):
            if string.find(line,'SUMMARY OLIGO DATA:') == -1:
                entry.append(line)
            else:
                content = string.join(entry,'')
                process(content, no_design_list,SYMBOL_LIST, TOLERANT, MASK, USER_UNIQUE, NUMBER_OLIGO, oligo_dup_file, out)
                entry =[]
                content =""
            line = fin.readline()
        fin.close()
        content = string.join(entry,'')
        process(content, no_design_list,SYMBOL_LIST, TOLERANT,MASK, USER_UNIQUE, NUMBER_OLIGO, oligo_dup_file,out)
        entry =[]
        content =""

        filecount = filecount +1
        currentfile = filename+str(filecount)
        
    out.close()
    oligo_dup_file.close()

    #no-design-output
    print "\nthe following genes were not designed for oligos because their sequences have too much low complexity regions.\n"
    print no_design_list
    out = open ('nodesign','w')
    for gene in no_design_list:
        out.write('>'+gene+'\n')

    out.close()
    
    print "\nThe genes with no design have been written to file \"nodesign\"\n"
