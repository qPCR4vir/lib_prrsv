#!/usr/bin/python
import string,sys,os,math
import fasta
import  Pick70_group
import Pick70_blast

TRUE = 1
FALSE = 0
UNDEFINED = 0
FLAGED = 'F'
BLAST =0
GC =1
REPEAT =2
S_COM =3
NULL =0
GB = math.pow(10,9)

#return nonredundant seq dic to with all the ids
#[seq:[id,id...]]
def all_id_seq_dic(id_dic):
    dic ={}
    for id in id_dic.keys():
        seq = id_dic[id]
        try:
            dic[seq].append(id)
        except KeyError:
            dic[seq] = [id]
    return dic
                                            
#input fastafile
#database already build
#file named temp will be changed during excution
def Output_Data_Dic (id,seq,data_Dic, dic_hit, finalout,filecount ):
    out = open ('~tempout','w')
    out.write('SUMMARY OLIGO DATA:\n')
    out.write(id+'\t'+ str(len(seq))+'\n')
    index = data_Dic.keys()
    index.sort()
    for key in index:
        out.write( str(key) + '\t'+ str(data_Dic[key][BLAST])+'\t'+ str(data_Dic[key][GC])+'\t'+ str(data_Dic[key][REPEAT])+'\t'+ str(data_Dic[key][S_COM])+'\t' + seq[key:key+OLIGOLEN] )
        #write genome target
        for entry in dic_hit[key]:
            hit_id, energy,pos = entry
            out.write( '\t'+ hit_id+ '\t' + str(energy)+'\t'+str(pos) )
        out.write('\n')
    out.close()

    if os.path.getsize(finalout+str(filecount)) > GB:
        filecount = filecount+1
        os.environ['finalout'] = finalout+str(filecount)
        os.system('touch $finalout')
        
    os.environ['finalout'] = finalout+str(filecount)
    os.system('cat ~tempout >> $finalout ; rm ~tempout')
    return filecount 
        
def Print_Data_Dic ( id, seq, data_Dic):
    print 'SUMMARY OLIGO DATA:'
    print id, len(seq)
    index = data_Dic.keys()
    index.sort()
    for key in index:
        print key, data_Dic[key], seq[key:key+OLIGOLEN]
    return
#gc
def run_gc(data_Dic, seq):
    import GC_compute

    keys = data_Dic.keys()
    keys.sort()
    if len(keys) == 0:
        return
    total, gc = GC_compute.Compute_GC(seq[keys[0]:keys[0]+OLIGOLEN])
    if total < 0.7*OLIGOLEN:
        data_Dic[keys[0]][GC] = FLAGED
    else:    
        data_Dic[keys[0]][GC] = gc*100 / float(total)
    for key in keys[1:]:
        if string.upper(seq[key-1]) in ['A','T','G','C']:
            total = total -1
            if string.upper(seq[key-1]) in ['G','C']:
                gc = gc -1 
        if string.upper(seq[key+OLIGOLEN-1]) in ['A','T','G','C']:
            total = total +1
            if string.upper(seq[key+OLIGOLEN-1]) in ['G','C']:
                gc = gc +1
        if total < 0.7*OLIGOLEN:
            data_Dic[key][GC] = FLAGED
        else:    
            data_Dic[key][GC] = gc *100 / float(total)


#self-complementation
def run_self_com (data_Dic,seq):
    os.environ['len'] = str(OLIGOLEN)
    fout = open("~temp","w")
    fout.write(seq)
    fout.close()
    os.environ['file'] = '~temp'
    f = os.popen('./code/SW $file 0 $len')
    retDat = f.read()
    f.close()
    os.system('rm $file');
    lines = string.split(retDat,'\n')[:-1]
    for i in range (0, len(lines)):
        data_Dic[i][S_COM] = int(lines[i])

#repeat
#check repeat within subsequence consist of SYMBOL_LIST, if omit, check whole sequence
def run_repeat(data_Dic,seq):
    os.environ['len'] = str(OLIGOLEN)
    fout = open("~temp","w")
    fout.write(seq)
    fout.close()
    os.environ['file'] = '~temp'
    f = os.popen('./code/LZW $file $len')
    retDat = f.read()
    f.close()
    os.system('rm $file');
    lines = string.split(retDat,'\n')[:-1]
    for i in range (0, len(lines)):
        data_Dic[i][REPEAT] = int(lines[i])
        
#blast 
def run_blast(data_Dic,id,  seq, start,end,dic_hit_info,Group_dic_same, DB):
    energy_dic  =  Pick70_blast.Blast(seq, id, Group_dic_same, DB,  OLIGOLEN, TRACEFLAG, STRAND)
    if energy_dic ==0 : #fail
        #error log
        ferror = open("error.log",'a')
        ferror.write(id+" failed on blast\n")
        ferror.close();
        return 0 #fail

    for pos in energy_dic.keys():
        data_Dic[pos][BLAST] = energy_dic[pos][0] #minimum energy
        if energy_dic[pos][0] != FLAGED:
            if energy_dic[pos][0] >=  0 :
                data_Dic[pos][BLAST] = UNDEFINED
            else:
                dic_hit_info[pos] = energy_dic[pos][1:]  #hit_info : [[id,energy,pos],[id,energy,pos],...] #energy is in original format
    return 1 #sucess

def compute (id,  seq, Data_Dic, start, end, Group_dic_same, dic_hit_info,DB):
    if start > end: #sequence length smaller than the oligo length
        return
    
    for i in range(start, end+1, 1):
        Data_Dic[i] =[UNDEFINED, UNDEFINED, UNDEFINED, UNDEFINED]
        dic_hit_info[i]=[]
        
    print id,len(seq)
    print "start gc counting ...."
    run_gc (Data_Dic,seq)
    print "start internal self binding calculation ...."
    run_self_com(Data_Dic,seq)
    print "start blasting ....."
    ret = run_blast(Data_Dic,id, seq,start, end,dic_hit_info, Group_dic_same, DB)
    if ret == 0:
        return 0
    print "start internal repeat calculation .... "
    run_repeat(Data_Dic,seq)

    return 1

if ( len(sys.argv[:]) < 8  or len(sys.argv[:]) > 9 )  and __name__ == '__main__':
     print 'USAGE python Pick70.py fastafile outfile DBfile TRACEFLAG(0:off,1:on) STRAND(1:plus 3:both) oligo_len masklower(yes, no) [groupfile_same]\n'
    
elif  __name__ == '__main__':
    print "OLIGO program in progress ... "

    OLIGOLEN = int(sys.argv[6])
    DB = sys.argv[3]
    
    TRACEFLAG = int(sys.argv[4]) #print to stdio
    STRAND = sys.argv[5]

    inputfile_name = sys.argv[1]
    infile = open(inputfile_name,'r')

    id_seq_dic={}
    #process input seqfasta file { id:seq,.......}
    id_seq_dic = fasta.fasta_dic(infile)
    ids = id_seq_dic.keys() #only input ids
    #does this need to mask lower case in the furthur blast
    if (sys.argv[7] != "yes"):
        for id in ids:
            id_seq_dic[id] = string.upper(id_seq_dic[id])

    #clean up redundant sequence in input file
    seq_id_dic = all_id_seq_dic(id_seq_dic)
    input_dup_list = []
    print "the following sequences are duplicated in the input file:",
    #error log
    ferror = open("error.log",'a')
    for values in seq_id_dic.values():
        if len(values) >1:
            for id in values[1:]:
                ids.remove(id)
                input_dup_list.append(id)
                print "WARNING:",id," is duplicated in sequence", values[0]
                ferror.write("WARNING:"+id+" is duplicated in sequence"+values[0]+"\n")
    print "done"
    ferror.close();
    
    #generate the blast db if the db files do not exist
    dbfile = sys.argv[3]
    if (os.access(dbfile+".nsq", os.R_OK)) and (os.access(dbfile+".nin", os.R_OK)) and (os.access(dbfile+".nhr", os.R_OK)): 
        pass
    else:
        os.environ['db'] = sys.argv[3]
        fr, fw,fe = os.popen3('formatdb -i $db  -p F')
        try:
            os.wait()
        except OSError:
            pass
        error = fe.read()
        if error!="":
            print error
            print "Error occured. Program terminated."
            sys.exit(1)
        fr.close()
        fw.close()
        fe.close()
    
    #reset output file
    outfile_name = sys.argv[2]
    filecount =0
    os.environ['finalout'] = outfile_name+str(filecount)
    os.system('touch $finalout')

    #process same_group filea to detect same sequence in input and genome files
    if len(sys.argv[:]) >= 9:
        groupfile = open(sys.argv[8],'r')
    else:
        groupfile = open('temp_group','w')
        Pick70_group.Generate_Group(inputfile_name, groupfile)
        groupfile = open('temp_group','r')
    Group_dic_same  = Pick70_group.All_Group_dic (groupfile) #{gene_id:[[hit_id, hit_strand, start, end...],[hit_id, hit_strand, start, end...]..., ...}

    ids.sort()
    
    for id in ids:
        seq = id_seq_dic[id]

        dic_hit_info ={}
        Data_Dic ={}
        
        ret = compute (id, seq, Data_Dic, 0, len(seq)- OLIGOLEN, Group_dic_same, dic_hit_info,DB)
        if ret == 1 :# success
            Print_Data_Dic (id, seq, Data_Dic)
            filecount  = Output_Data_Dic (id,seq,Data_Dic, dic_hit_info, outfile_name,filecount)
        else:
            #error log
            ferror = open("error.log",'a')
            ferror.write("WARNING: "+id+" is "+ str(len(seq))+" bp.\n")
            ferror.close();
                
    if len(input_dup_list)>0:
        print "duplicated sequences in input file:", input_dup_list

    os.system("rm ~out.blastout")
    print "OLIGO program in progress ... DONE"
    sys.exit(0)

