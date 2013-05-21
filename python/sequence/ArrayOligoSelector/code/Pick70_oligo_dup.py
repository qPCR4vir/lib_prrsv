import string, sys,os
import blast_parse
import fasta
import GC_compute


def oligo_dup_output(out,parse_dic):
    keys= parse_dic.keys()
    keys.sort()
    for key in keys:
        Query = key
        if len(parse_dic[key]) > 1:
            temp_list = []
            hit_list = parse_dic[key][1:-1]
            alignment = parse_dic[key][-1]
            for i in range (0, len(hit_list)):
#                hit_id, hit_score, hit_E, hit_len, hit_identical, hit_cover = hit_list[i]
                hit_id, hit_score = hit_list[i]
                align = string.split(alignment [i], ' Score =')[1] #1st HSP
                plusline, minusline, p_start, m_start, p_end, m_end  =  blast_parse.parse_align(align)
                os.environ['plusline']= plusline
                os.environ['minusline']= minusline

                f = os.popen('./code/energy \"$plusline\" \"$minusline\" 1')
                retData = f.read()
                f.close()
                energy =  float(string.split(string.strip(retData),'\n') [-1])
                if energy >= -10:
                    break
#                out.write(hit_id+' '+str(hit_identical)+'\t')
#                out.write(hit_id+' '+str(hit_E)+'\t')
#                out.write(hit_id+' '+str(energy)+'\t')

                temp_list.append([energy,hit_id])
            temp_list.sort()
            for i in range (0, len(temp_list)):
                energy, hit_id = temp_list[i]
#                out.write(hit_id+' '+str(energy)+'\t'+str(hit_identical)+'\t')
                out.write(hit_id+'\t'+str(energy)+'\t')
        out.write('\n')


if len (sys.argv[:]) != 5:
    print 'python Pick70_oligo_dup.py oligofastafile(inputfasta) genome_file(inputfasta)  oligodupfile(output) genome/gene(3,1)\n'
    
else:
#    OLIGOLEN = 70
    os.environ['inputf'] = sys.argv[1]
    os.environ['genomef'] =sys.argv[2]
    genomefile = sys.argv[2]
    out = open(sys.argv[3],'w')
    os.environ['STRAND'] = sys.argv[4]

    #test if the program can procede
    if os.access(genomefile+'.nhr', os.R_OK):
        if not os.access(genomefile+'.nhr', os.W_OK):
            print "The same program is running.  the current program can not procede"
            sys.exit()


    #redo the sw_score ,TM part
    inputf = open(sys.argv[1],'r')
    input_dic = fasta.fasta_dic(inputf)

    sw_dic ={}
    gc_dic ={}
    repeat_dic={}
    
    for key in input_dic.keys():
        seq =input_dic[key]

        os.environ['len'] = str(len(seq))
        fout = open("~temp","w")
        fout.write(seq)
        fout.close()
        os.environ['file'] = '~temp'

        f = os.popen('./code/SW $file 0 $len')
        retData =f.read()
        f.close()
        lines = string.split(retData,'\n')[:-1]
        sw_dic[key] = int(string.strip(lines[0]))

        f = os.popen('./code/LZW $file $len')
        retData = f.read()
        f.close()
        lines = string.split(retData,'\n')[:-1]
        repeat_dic[key] = int(lines[0])

        os.system('rm $file');

        total, gc = GC_compute.Compute_GC(seq)
        gc_dic[key] = gc * 100 / float(total)

        
    MAX = 20
    os.environ['MAX'] = str(MAX)

    #redo the genome match

    os.system('formatdb -p F -i $genomef')
    os.system('chmod a-wx $genomef'+'.*')
    os.system('blastall -p blastn -a 2 -i $inputf -d $genomef -S $STRAND  -F F  -v $MAX -b $MAX -G 5 -E 0 -o oligo_blastout')
    os.system('chmod u+w $genomef'+'.*')
    os.system('rm $genomef'+'.*')

    resultf = open('oligo_blastout','r')

    parse_dic={}
    content = []
    line = " "
    while (line !=""):
        line = resultf.readline()
        if line[:6] == "BLASTN":
            if content ==[]:
            	continue 
            blast_parse.parse_alignment(content, parse_dic,1)
            Query = parse_dic.keys()[0]
            print Query
            out.write('>'+Query+'\t'+str(gc_dic[Query])+'\t'+str(repeat_dic[Query])+'\t'+str(sw_dic[Query])+'\n')
            oligo_dup_output(out, parse_dic)
            content =[line]
            parse_dic={}
        else:
            content.append(line)

    blast_parse.parse_alignment(content, parse_dic,1)
    if parse_dic != {}:
        Query = parse_dic.keys()[0]
        out.write('>'+Query+'\t'+str(gc_dic[Query])+'\t'+str(repeat_dic[Query])+'\t'+str(sw_dic[Query])+'\n')
        oligo_dup_output(out, parse_dic)
    
    resultf.close()
    out.close()
