codon ={ 'TTT':  'F',      'TCT':  'S' ,      'TAT':  'Y',      'TGT':  'C',
         'TTC':  'F',      'TCC':  'S' ,      'TAC':  'Y',      'TGC':  'C',
         'TTA':  'L',      'TCA':  'S' ,      'TAA':  '*',      'TGA':  '*',
         'TTG':  'L',      'TCG':  'S',       'TAG':  '*',      'TGG':  'W',
         
         'CTT':  'L',      'CCT':  'P' ,      'CAT':  'H',      'CGT':  'R',
         'CTC':  'L',      'CCC':  'P' ,      'CAC':  'H',      'CGC':  'R',
         'CTA':  'L',      'CCA':  'P' ,      'CAA':  'Q',      'CGA':  'R',
         'CTG':  'L',      'CCG':  'P' ,      'CAG':  'Q',      'CGG':  'R',

         'ATT':  'I',      'ACT':  'T' ,      'AAT':  'N',      'AGT':  'S',
         'ATC':  'I',      'ACC':  'T' ,      'AAC':  'N',      'AGC':  'S',
         'ATA':  'I',      'ACA':  'T' ,      'AAA':  'K',      'AGA':  'R',
         'ATG':  'M',      'ACG':  'T' ,      'AAG':  'K',      'AGG':  'R',

         'GTT':  'V',      'GCT':  'A' ,      'GAT':  'D',      'GGT':  'G',
         'GTC':  'V',      'GCC':  'A' ,      'GAC':  'D',      'GGC':  'G',
         'GTA':  'V',      'GCA':  'A' ,      'GAA':  'E',      'GGA':  'G',
         'GTG':  'V',      'GCG':  'A' ,      'GAG':  'E',      'GGG':  'G',

         'ttt':  'F',      'tct':  'S' ,      'tat':  'Y',      'tgt':  'C',
         'ttc':  'F',      'tcc':  'S' ,      'tac':  'Y',      'tgc':  'C',
         'tta':  'L',      'tca':  'S' ,      'taa':  '*',      'tga':  '*',
         'ttg':  'L',      'tcg':  'S',       'tag':  '*',      'tgg':  'W',
         
         'ctt':  'L',      'cct':  'P' ,      'cat':  'H',      'cgt':  'R',
         'ctc':  'L',      'ccc':  'P' ,      'cac':  'H',      'cgc':  'R',
         'cta':  'L',      'cca':  'P' ,      'caa':  'Q',      'cga':  'R',
         'ctg':  'L',      'ccg':  'P' ,      'cag':  'Q',      'cgg':  'R',

         'att':  'I',      'act':  'T' ,      'aat':  'N',      'agt':  'S',
         'atc':  'I',      'acc':  'T' ,      'aac':  'N',      'agc':  'S',
         'ata':  'I',      'aca':  'T' ,      'aaa':  'K',      'aga':  'R',
         'atg':  'M',      'acg':  'T' ,      'aag':  'K',      'agg':  'R',

         'gtt':  'V',      'gct':  'A' ,      'gat':  'D',      'ggt':  'G',
         'gtc':  'V',      'gcc':  'A' ,      'gac':  'D',      'ggc':  'G',
         'gta':  'V',      'gca':  'A' ,      'gaa':  'E',      'gga':  'G',
         'gtg':  'V',      'gcg':  'A' ,      'gag':  'E',      'ggg':  'G'}

       
start_codon ={'TTG':  'M' , 'CTG': 'M',   'ATG':  'M', 'ttg':  'M' , 'ctg': 'M',   'atg':  'M'}
              
compen ={'A':'T', 'T':'A', 'G':'C', 'C':'G', 'U':'A',
	 'M':'K', 'R':'Y', 'W':'W', 'S':'S', 'Y':'R', 'K':'M',
	 'V':'B', 'H':'D', 'D':'H', 'B':'V',
	 'N':'N', 'X':'X',
         'a':'t', 't':'a', 'g':'c', 'c':'g','u':'a',
         'm':'k', 'r':'y', 'w':'w', 's':'s', 'y':'r', 'k':'m',
         'v':'b', 'h':'d', 'd':'h', 'b':'v',
         'n':'n', 'x':'x'}

code= {'A':['A'], 'T':['T'], 'G':['G'], 'C':['C'],
       'M':['A','C'],'R':['A','G'],'W':['A','T'], 'S':['C','G'], 'Y':['C','T'], 'K':['G','T'],
       'V':['A','C','G'], 'H':['A','C','T'], 'D':['A','G','T'], 'B':['C','G','T'],
       'N':['A','T','G','C'],'X':['A','T','G','C'],
       'a':['a'], 't':['t'], 'g':['g'], 'c':['c'],
       'm':['a','c'],'r':['a','g'],'w':['a','t'], 's':['c','g'], 'y':['c','t'], 'k':['g','t'],
       'v':['a','c','g'], 'h':['a','c','t'], 'd':['a','g','t'], 'b':['c','g','t'],
       'n':['a','t','g','c'],'x':['a','t','g','c']}
      










