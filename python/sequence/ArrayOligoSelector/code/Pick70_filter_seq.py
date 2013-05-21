import fasta , sys, os
import string

if len(sys.argv[:]) != 5:
    print "USAGE: python Pick70_filter.py fin_input fin_db oligolen"
    print "       This is a pre-run module to varify the sequence format and the oligolength" 
else:
    print "Checking file format and all sequences ..."
    seq_size ={} #id: length of sequence
    if os.path.isfile(sys.argv[1]):
    	fin = open(sys.argv[1],'r')
    else:
	sys.exit(1)

    result = fasta.check_format_fasta(fin, seq_size)
    fin.close()
    if result !=1 :
        if result == 0: #wrong format
            sys.exit(1)
        if result == 2: #conflict sequence size
            sys.exit(4)
        if result == 3: #idwith space
            sys.exit(5)

    if os.path.isfile(sys.argv[2]):
	fin = open(sys.argv[2],'r')
    else:
	sys.exit(2)

    result = fasta.check_format_fasta(fin,seq_size)
    fin.close()
    if result!=1:
        if result == 0:
            sys.exit(2)
        if result == 2:
            sys.exit(4)
        if result == 3: #idwith space
            sys.exit(5)
        
    try:
        oligo_len = string.atoi(sys.argv[3])
        if oligo_len <= 0:
            sys.exit(3)
    except ValueError:
	sys.exit(3)

    if ( sys.argv[4] not in ["yes","no"] ):
        sys.exit(6)
         
    print "Checking file format and all sequences ... DONE"
    sys.exit(0)







