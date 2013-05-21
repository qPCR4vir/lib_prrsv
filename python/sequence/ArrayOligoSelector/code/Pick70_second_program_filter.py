import string, sys, os

if len(sys.argv[:]) not in [4,5,8]:
    print "usage: python Pick70_second_program_filer.py GCpercentage(float) oligo_length(integer) numberOLigo(integer) [cutoff (flaot) ] [masklength(integer) mask_symbol(combination of ATGC) mask_tolerance(0-1,float)]"
else:
    #parse no 1. GC percentage
    try:
        value = string.atof(sys.argv[1])
        if value >100.0 or value <0:
            print "please enter an valid GC percentage number between 0 and 100"
            sys.exit(1)
    except ValueError:
        print "please enter a integer or floating point number for GC percentage"
        sys.exit(1)
        
    #parse no 2, oligo length
    try:
        value = string.atoi(sys.argv[2])
        if value <=0:
            print "please enter a valid oligo length"
            sys.exit(1)
    except ValueError:
        print "please enter an integer for oligo length"
        sys.exit(1)

    #parse no 3, oligo number
    try:
        value = string.atoi(sys.argv[3])
        if value <=0:
            print "please enter a valid oligo number"
            sys.exit(1)
    except ValueError:
        print "please enter an integer for oligo number"
        sys.exit(1)

    if len(sys.argv[:]) in [5,8]:
        #parse no 4, energy
        try:
            value = string.atof(sys.argv[4])
            if value > 0:
                print "please enter a valid binding energy cutoff"
                sys.exit(1)
        except ValueError:
            print "please enter a floating number for binding energy cutoff"
            sys.exit(1)

        
    if len(sys.argv[:]) == 8:
        #parse no 5, maks length
        try:
            value = string.atoi(sys.argv[5])
            if value <0:
                print "please enter a valid mask length"
                sys.exit(1)
        except ValueError:
            print "please enter an integer for mask length"
            sys.exit(1)
        
        #parse no 6, maks symbol
        for i in range (0, len(sys.argv[6])):
            letter = sys.argv[6][i]
            combo =[]
            if string.upper(letter) not in ['A','T','G','C','N']:
                print "please enter a valid combination of bases"
                sys.exit(1)
            if letter in combo:
                print "please do not enter duplication in the base combination"
                sys.exit(1)
            else:
                combo.append(letter)
            
        #parse no 7, maks tolerance
        try:
            value = string.atof(sys.argv[7])
            if value <0 or value >=1:
                print "please enter a valid mask tolerance number between 0 and 1 representing 0% to 100% tolerance"
                sys.exit(1)
        except ValueError:
            print "please enter a floating point number for the mask tolerance"
            sys.exit(1)
        
    sys.exit(0)
