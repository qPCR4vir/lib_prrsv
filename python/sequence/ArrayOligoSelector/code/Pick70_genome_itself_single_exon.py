import string, os, sys
import math
import blast_parse
import fasta
import Pick70_mask

length_diff_tolerant = 20
bigSize = 50

#generate the rule dictionary as id:okid_list using the position of the input list as ids
#the rule is the query start and end position of the any two segment can not include each other
#the covereage_dic is a dictionary with id of hiti and values are a vector of vector, each vector represents the alignments ids the cover that position. 
def rules(hit_candidate,querylen):
    # each entry is the list : [hitid, strand,qstart,qend, hstart,hend]
    rule_dic={}
    coverage_dic= {} 
    bigPieceInterest =[] 
    for i in range(0, len(hit_candidate)):
        rule_dic[i] =[]
        hiti, strandi, qstarti,qendi, hstarti, hendi  = hit_candidate[i][:]

        #build a dictionary with id of hiti and values are a vector of vector, each vector represents the alignments ids the cover that position. 
        try:
            for n in range(qstarti,qendi):
                if len(coverage_dic[hiti][n]) <2:
                    coverage_dic[hiti][n].append(i)
                                
        except KeyError:
            coverage_dic[hiti] = []
            for n in range(0,querylen):
                coverage_dic[hiti].append([])
                
        if qendi - qstarti > bigSize:
            if i not in bigPieceInterest:
                bigPieceInterest.append(i)
        if  i not in bigPieceInterest:
            continue
        
        for j in range(i+1, len(hit_candidate)):
            hitj, strandj, qstartj,qendj,hstartj, hendj = hit_candidate[j][:]
            #the two hits have to be from the same contig
            if hiti != hitj:
                break
            #the two hits have to be from the same strand
            if strandi != strandj:
                continue

            #test if two query segments include each other
            if (qstarti >= qstartj) and (qendi <= qendj): #i is inside j
                continue
            if (qstarti <= qstartj) and (qendi >= qendj): #j is inside i
                continue

            #test if two hits segments include each other
            if (hstarti >= hstartj) and (hendi <= hendj): #i is inside j
                continue
            if (hstarti <= hstartj) and (hendi >= hendj): #j is inside i
                continue

            minimum = min ((qendj-qstartj), (qendi-qstarti))

            if ( (qstarti<qstartj) and (qstartj<qendi) and ( (qendi-qstartj) > 0.7 * minimum)): #overlap more than 80% of i
                continue
            if ( (qstartj<qstarti) and (qstarti<qendj) and ( (qendj-qstarti) > 0.7 * minimum)): #overlap more than 80% of j
                continue
            
            if abs (hstarti - hstartj) > 3000000:
                continue

            rule_dic[i].append(j)
            if j not in bigPieceInterest:
                bigPieceInterest.append(j)

    must_coverage =[]
    index = coverage_dic.keys()

    for id in index:
        for n in range(0, len(coverage_dic[id])):
            if len(coverage_dic[id][n]) == 1:
                alignid = coverage_dic[id][n][0]
                if alignid not in must_coverage:
                    must_coverage.append(alignid)

    return rule_dic, must_coverage
 
#generate a list of all the combinations that do not controdict with the rule_dic:
def generate_combination(rule_dic, must_cover_list, hit_candidate, seqlen): #id:okids
    max = len(rule_dic.keys())
    queue =[] #[[[selected_ids], [left_ids]],.....]
    ret_list =[]

    #add the rule_dic into the queue first
    for key in rule_dic.keys():
        hiti, strandi, qstarti,qendi, hstarti, hendi  = hit_candidate[key]
        if abs(qendi - qstarti) > bigSize:
            queue.append([[key], rule_dic[key]])

    while (len(queue)!=0):
        temp_queue=[]
        while (len(queue)!=0):
            #pop the last entry , add the entry to return,
            selected_ids, left_ids = queue.pop()
            ret_list.append(selected_ids)            
            #if the size of selected_ids not reach the max, add more the selected_ids from left_ids
            for pick_id in left_ids:
                new_selected = selected_ids[:]
                new_selected.append(pick_id)
                new_left=[]
                #eliminate from left_ids the ones that do not match pick_id
                for left in left_ids:
                    if ( left  in rule_dic[pick_id]):
                        new_left.append(left)
                total_length =0
                count = 0
                numberExon =100  #assume of number of exons
                if seqlen <200:
                    numberExon = 5

                for id in new_selected:
                    count = count +1
                    if count > numberExon:
                        break
                    hiti, strandi, qstarti,qendi, hstarti, hendi  = hit_candidate[id]
                    total_length = total_length + qendi - qstarti

                for id in new_left:
                    count = count +1
                    if count > numberExon:
                        break
                    hiti, strandi, qstarti,qendi, hstarti, hendi  = hit_candidate[id]
                    total_length = total_length + qendi - qstarti

                #the alignment ids must have at least one must_cover_list
                found =0
                for id in new_selected + new_left:
                    if id in must_cover_list:
                        found = 1
                        break

                if (found ==0):
                    continue

                if total_length < seqlen - length_diff_tolerant:
                    continue
                
                temp_queue.append([new_selected,new_left])

                hiti, strandi, qstarti,qendi, hstarti, hendi  = hit_candidate[pick_id]
                if (qendi - qstarti) > 150:
                    break
        queue = temp_queue
    
    return ret_list
                
    
#[[hitid, strand, pstart, pend, hstart, hend],...]        
#constrains are:
# 1. the query segments add up together to the query length, but due to blast program , the following situatin might happen: 1,758,757,1000, but the length of the combined segment will always bigger than the query length. However, if the exon segments are too small(<20), they might be missed by the blastn program, therefore, the requirement is relaxed as the sum of all segment is bigger than the size of the query -20bp. 

# 2. the query segments has structure of 1,300,301,456,457,700, but due to blast program , the following situatin might happen: 1,758,757,1000, allow overlapping, but not overshadow, that is do not allow 40,100, 50,60 , also not allowed in missed connectin: that is 10,50,60,100.  

# 3. hit segments are strand compatible  (do not need to do, because function rules has done that)

# 4. hit segments has no overlapping structure

# 5. same contig (do not need to do anything, because function rules has garanteed that)
def constrains( query_len, candidate_list):
    #test constrain no. 1
    query_seg_list =[]
    size =0
    for list in candidate_list:
        strand, qstart,qend, hstart, hend = list[1:]
        size = size + qend - qstart +1
        query_seg_list.append([qstart,qend, hstart,hend,strand])
    if size < query_len - length_diff_tolerant:
        return 0

    #test constrain no. 2 no. 4
    query_seg_list.sort()
    for i in range(0,len(query_seg_list)-1):
        if query_seg_list[i][1] <  query_seg_list[i+1][1] and query_seg_list[i][0] < query_seg_list[i+1][0] :
            pass
        else:
            return 0
        if strand=='+':
            if query_seg_list[i][2] <  query_seg_list[i+1][2] and query_seg_list[i][3] < query_seg_list[i+1][3] :
                pass
            else:
                return 0
        if strand =='-':
            if query_seg_list[i][2] >  query_seg_list[i+1][2] and query_seg_list[i][3] > query_seg_list[i+1][3] :
                pass
            else:
                return 0
            
    #test constrain no.3 (no need)
    return 1

#def adhocfilter(input): #, output):#filter seq with AT rich region
#    fin = open(input,'r')
#    seq_dic= fasta.fasta_dic(fin)
#    new_seq={}
    
#    for id in seq_dic.keys():
#        seq = string.upper(seq_dic[id])
#        start = 0
#        AT_length=0
#        cutoff = 20
#        mask_list=[]
#        for i in range(0,len(seq)):
#            letter = seq[i]
#            if letter not in [ 'G','C']:
#                if AT_length ==0:
#                    start = i
#                AT_length= AT_length +1
#            else:
#                if AT_length >= cutoff:
#                    mask_list.append([start,i])
#                AT_length =0

#        seq_list=[]
#        for i in range(0,len(seq)):
#            seq_list.append(seq[i])
#        for start,end in mask_list:
#            for i in range(start,end+1):
#                seq_list[i]=string.lower(seq[i]) #ad-hoc masking
#        masked_seq = string.join(seq_list,'')
#        new_seq[id] = masked_seq
#    return new_seq

def parsing_each_query(content,fout):
    #1. have to be 100 similarity
    #2. different segment add up together to be the complete length
    #3. check strand constrains
    #4. check no overlapping constrains
    content = string.strip(string.join(content,''))
    if content =='':
        return
    
    parse_dic = blast_parse.parse_query_alignment(content)
    
    query = parse_dic.keys()[0]
    print query,
    hit_candidate_segment =[] #[[hitid, strand, pstart, pend, hstart, hend],...]
    group_list =[]
    query_len= parse_dic[query][0]
    hit_dic = parse_dic[query][1] #only two
    hit_list = hit_dic.keys()

    for hitid in hit_list:
        #to get the id
        hit_align_list = hit_dic[hitid]

        first_100 =0
        for seg_list in hit_align_list:
            strand, qstart, qend, hstart, hend, percentage = seg_list
            if percentage < 100 :#99.5:
                continue
            if abs(qend - qstart) +1 <query_len:
                continue
            hit_candidate_segment.append([hitid, strand,qstart,qend, hstart,hend])
         
    #write to file of the group list
    fout.write('['+query+']\n')
    #write the query itself
    fout.write(query+'\n')
    fout.write("%s %s %s %s %s\n" % ( '+','1',str(query_len),'1',str(query_len)))
    hit_dic={}

    if gfdir !="":
        string_start = len(gfdir)+1
    else:
        string_start = 0
        
    for entry in hit_candidate_segment:
        hitid, hitstrand, qstart,qend, hstart,hend = entry
        hitid = hitid[string_start:]
        if gfdir !="":
            string_end = string.find(hitid, ".nib:")
            hitid = hitid[:string_end]
        fout.write(hitid+'\n')
        fout.write("%s %s %s %s %s\n" % (hitstrand, str(qstart), str(qend), str(hstart), str(hend)))
    fout.write('\n')
    print "done"

if len(sys.argv[:]) != 6:
    print "USAGE: python Pick70_genome_itself.py inputfasta genomefasta groupfile(fout) version(BLAST/BLAT/GFCLIENT) masklower(yes, no)"
else:
    print "Identifying the input sequences' genomic target ... "

    input = sys.argv[1]
    genome = sys.argv[2]
    fout =open( sys.argv[3],'w')
    VERSION = string.upper(sys.argv[4])
    os.environ['genome'] = genome
    os.environ['input'] = input
    gfdir =""
        
    #convert to upper case if no masking user lower case
    finput = open(input,'r')
    seq_dic= fasta.fasta_dic(finput)
    index= seq_dic.keys()
    if sys.argv[5] !="yes":
        for id in index:
            seq_dic[id] = string.upper(seq_dic[id])
    tempinput = '~ftempinput'
    ftempinput =open(tempinput,'w')
    for id in index:
        ftempinput.write('>'+id+'\n')
        ftempinput.write(seq_dic[id]+'\n')
    ftempinput.close()

    #blast version
    if (VERSION =="BLAST"):
        #generate the blast db if the db files do not exist
        if (os.access(genome+".nsq", os.R_OK)) and (os.access(genome+".nin", os.R_OK)) and (os.access(genome+".nhr", os.R_OK)):
            pass
        else:
            fr, fw, fe = os.popen3('formatdb -i $genome -p F')
            errormessage = fe.read()
            fr.close()
            fw.close()
            fe.close()
            if errormessage!="":
                print errormessage
                print "program terminated"
                sys.exit(1)

        maskinput =open('~maskedinput','w')
        for id in index:
            masked_seq = Pick70_mask.mask_longAT(seq_dic[id],"lower")
            maskinput.write('>'+id+'\n')
            maskinput.write(masked_seq+'\n')
        maskinput.close()

        os.environ['input'] = '~maskedinput'
        fr, fw, fe = os.popen3('blastall -i "$input" -d $genome -a $cpunum -p blastn -v 50 -b 50 -e 1  -o ~tempblatout  -F "m D" -U')

    elif (VERSION =="BLAT"):
        Pick70_mask.mask_file(tempinput,input+".mask","dust", "lower")
        #blat version
        fr,fw,fe  = os.popen3('blat "$genome" "$input".mask -out=blast -qMask=lower ~tempblatout')
        
    elif (VERSION =="GFCLIENT"):
        Pick70_mask.mask_file(tempinput,input+".mask","dust", "lower")
        #gfCleint version
        fgfinfo = open(".gfServer_information",'r')
        gfServer= string.strip(fgfinfo.read())
        fgfinfo.close()
        fieldnum =  string.split(string.strip(gfServer))
        if len(fieldnum) !=3:
            print "Wrong Format of .gfServer_information file. It has only",fieldnum,"fields"
            print "program terminated"
            sys.exit(1)
        os.environ['gfServer'] = gfServer
        fr,fw,fe  = os.popen3('./code/gfClient $gfServer "$input"  ~tempblatout  -qMask=lower -out=blast')
        gfdir = string.split(string.strip(gfServer))[-1]
    else:
        print "Unknown method for identifying genomic targets. Program terminated."
        sys.exit(1)

    os.system("rm ~ftempinput")
    
    try:
        os.wait()
    except OSError:
        pass
    
    errormessage = fe.read()
    fr.close()
    fw.close()
    fe.close()
    if errormessage != "":
        print "Error message from:", errormessage
        if  (VERSION =="GFCLIENT"):
            print "If you choose to use gfclient, you need first set up the gfServer, and enter the server information in the \".gfServer_information\" file.  It is VERY IMPORTANT that you set up the gfServer from the directory where the .nib files locate. "
            print "Program terminated."
            sys.exit(1)
        if (VERSION == "BLAT"):
            print "Program terminated."
            sys.exit(1)
            
    fparse = open('~tempblatout', 'r')
    parse_dic = {}
    content = []
    line = " "
    ALIGN =1
    while (line !=""):
        line = fparse.readline()
        if line[:6] == "BLASTN":
            parse_dic={}
            parsing_each_query(content,fout)
            content =[line]
        else:
            content.append(line)
    parse_dic={}
    parsing_each_query(content,fout)
    fparse.close()
    os.system('rm ~tempblatout')
    print "Identifying the input sequences' genomic target ... DONE"
    sys.exit(0)
    fout.close()
