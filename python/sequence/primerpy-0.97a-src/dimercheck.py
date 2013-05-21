#!/usr/bin/env python
"""dimercheck.py
examine dimer formation of primer pairs 
return structures in plain text or html
support degenerate primers
support batch mode
non-base characters will be skipped

dimer structures are defined like:
[  (set([1, 2, 3, 4]), set([16, 17, 18, 15])), 
(set([12, 13, 14, 15]), set([4, 5, 6, 7])), 
(set([13, 14, 15]), set([0, 1, 2]))  ]

last update: 11/19/2007
"""

wcdict={'A':'T',
        'T':'A',
        'C':'G',
        'G':'C'}

IUPAC={'R': ['A', 'G'],
        'Y': ['C', 'T'],
        'S': ['G', 'C'],
        'W': ['A', 'T'],
        'K': ['G', 'T'],
        'M': ['A', 'C'],
        'B': ['C', 'G', 'T'],
        'D': ['A', 'G', 'T'],
        'H': ['A', 'C', 'T'],
        'V': ['A', 'C', 'G'],
        'N': ['A', 'T', 'C', 'G']
        }


class SequenceUtil:
    def __init__(self):
        self.note = 'the smaller total_score the better; <br>for QPCR applications, a total_score under 6 is advised.'
    def basecheck(self, seq):
        """check if seq contain illegal bases"""
        for item in seq:
            if item not in ['A', 'T', 'C', 'G']:
                raise ValueError
            
    def deambiguity(self, seq):
        """translate ambigious seq to a batch of all possible combinations; non-base characters will be skipped"""
        global IUPAC
        iupac_keys = IUPAC.keys()
        batch = ['', ]
        for base in seq:
            N = len(batch)
            if base in ['A', 'T', 'C', 'G']:
                for i in range(N):
                    batch[i] += base
            elif base in iupac_keys:
                baselist = IUPAC[base]
                forks = len(baselist)
                batch *= forks
                for i in range(N):
                    for j in range(forks):
                        batch[i*forks+j] += baselist[j]
        return batch
    def MakeBatch(self, seq1, seq2):
        """return a list of [(seq1, seq2), ...]"""
        batch = []
        list_sense = self.deambiguity(seq1)
        list_antisense = self.deambiguity(seq2)
        for a in list_sense:
            for b in list_antisense:
                batch.append((a, b))
        return batch
    def BatchCheck(self, inlist):
        summary = ''
        s = ''
        for apair in inlist:
            pp = PrimerPair(apair[0], apair[1])
            pp.compute()
            pp.display()
            summary += pp.header + '<br>'
            s += '<h4>' + pp.header + '</h4>'
            s += '<pre>' + pp.structures + '</pre>'
        html_text = '<html><body><h1>BatchCheck Result</h1>' +self.note + '<p>' + summary + '</p><hr>' + s + '</body></html>'
        return html_text
    def BatchCheck_single(self, inlist):
        """Same as BatchCheck, but a list of unpaired single primers"""
        summary = ''
        s = ''
        for a in inlist:
            sp = SinglePrimer(a)
            sp.compute()
            sp.display()
            summary += sp.header + '<br>'
            s += '<h4>' + sp.header + '</h4>'
            s += '<pre>' + sp.structures + '</pre>'
        html_text = '<html><body><h1>BatchCheck Result</h1>' +self.note + '<p>' + summary + '</p><hr>' + s + '</body></html>'
        return html_text

    def html_write(self, html_text, html_out):
        out = open(html_out, 'w')
        out.write(html_text)
        out.close()
    def degenerate_check(self, seq1, seq2):
        batch = self.MakeBatch(seq1, seq2)
        html_text = self.BatchCheck(batch)
        return html_text
    def degenerate_check_single(self, seq):
        batch = self.deambiguity(seq)
        html_text = self.BatchCheck_single(batch)
        return html_text
    
    def reformat(self, text):
        """reformat the wrapped primer3 output to html, with dimerchecked"""
        plist, size = [], []
        s = '<h2>Dimer check report</h2>'
        x = text.find('PRIMER PICKING RESULTS FOR')
        html_text = '<html><body>'
        w = text[x:].splitlines()
        if len(w) > 17: #check for failed runs
            summary =  '<h2>PrimerPy RUN SUMMARY</h2><pre>' + text[:x].rstrip() + '</pre>' + self.note
            for line in w:
                if 'LEFT PRIMER' in line:
                    x = line.split()
                    a, posa = x[-1], int(x[-7])
                elif 'RIGHT PRIMER' in line:
                    x = line.split()
                    b, posb = x[-1], int(x[-7])
                    plist.append((a, b))
                    size.append(posb-posa+1)
            i = 0
            for apair in plist:
                i += 1
                pp = PrimerPair(apair[0], apair[1])
                pp.compute()
                pp.display()
                summary += '<p>primer_pair_' + str(i) + '<br>   <b>' + apair[0] + '<br>   ' + apair[1] + '</b><br>' + 'amplicon size = ' + str(size[i-1]) + '<br>total_score = ' + str(pp.total_score) + '</p>'
                s += '<h4>' + pp.header + '</h4>'
                s += '<pre>' + pp.structures + '</pre>'
        else:
            summary = '<h2>Failed run.</h2>'
        html_text += summary + '<hr><h2>Original primer3 report</h2><pre>' + text.replace('>', '&gt;').replace('<', '&lt;') + '</pre>' + '<hr>' + s + '</body></html>'
        return html_text




class PrimerFunctions:
    """superclass for primer functions"""
    def rc(self, instr):
        """get reverse complementary seq"""
        global wcdict
        N = len(instr)
        s = ''
        for i in range(N):
            s += wcdict[instr[-i-1]]
        return s
    def reverse(self, instr):
        s = ''
        for i in range(len(instr)):
            s += instr[-i-1]
        return s
    def findall(self, a, seq):
        """find all occurance of a in seq, return list of positions"""
        pos = []
        i = 0
        p = seq.find(a)
        while p != -1:
            pos.append(p+i)
            i += p+1
            p = seq[i:].find(a)
        return pos
    def triplet_seed(self, str1, str2):
        """both seq are 5'-3', check complementariness by triplets in str1"""
        N = len(str1)
        rc2 = self.rc(str2)
        triplethit = []
        for i in range(N-3):
            hit = self.findall(str1[i:i+3], rc2)
            for p in hit:
                triplethit.append((i, p))
        #print triplethit
        return triplethit
    #expand the seeds to maximum, then remove the redundant
    def expand(self, triplethit, str1, str2):
        #output as [([1, 2, 3, 4], [15, 16, 17, 18]), ..., ()]
        raw_result, result = [], []
        rc2 = self.rc(str2)
        if len(triplethit) > 0:
            for item in triplethit:
                #e.g. in [(1, 15), (2, 16), (12, 4), (13, 0), (13, 5)]
                a, b = item
                pairs = (set([a, a+1, a+2]), set([b, b+1, b+2]))
                downstr1, downrc2 = a-1, b-1
                upstr1, uprc2 = a+3, b+3
                #match up/down stream along the shorter ends, end at mismatch
                if a > b:
                    while upstr1 < len(str1) and str1[upstr1] == rc2[uprc2]:
                        pairs[0].add(upstr1)
                        pairs[1].add(uprc2)
                        upstr1 += 1
                        uprc2 += 1
                    while downrc2 > -1 and str1[downstr1] == rc2[downrc2]:
                        pairs[0].add(downstr1)
                        pairs[1].add(downrc2)
                        downstr1 -= 1
                        downrc2 -= 1
                else:
                    while uprc2 < len(rc2) and str1[upstr1] == rc2[uprc2]:
                        pairs[0].add(upstr1)
                        pairs[1].add(uprc2)
                        upstr1 += 1
                        uprc2 += 1
                    while downstr1 > -1 and str1[downstr1] == rc2[downrc2]:
                        pairs[0].add(downstr1)
                        pairs[1].add(downrc2)
                        downstr1 -= 1
                        downrc2 -= 1
                raw_result.append(pairs)
            for item in raw_result:
                if item not in result:
                    result.append(item)
        return result

    def list_structure(self, str1, str2):
        trihit = self.triplet_seed(str1, str2)
        return self.expand(trihit, str1, str2)

    def score(self, pairedlist, str1_len):
        """score = longest stretch of base-pairing + 3'-complement penalty"""
        score3end, score_longest, worst_list = 0, 0, []
        if len(pairedlist) > 0:
            score3end, worst_list = self.complement3end(pairedlist, str1_len)
            score_longest = max([len(item[0]) for item in pairedlist])
            for item in pairedlist:
                if len(item[0]) == score_longest and item not in worst_list:
                    worst_list.append(item)
            #r = score3end + score_longest
        return score3end, score_longest, worst_list
    def complement3end(self, pairedlist, str1_len):
        #str1_len = len(str1), 3'end of rc(str2) is 0
        score3end, worst_list = 0, []
        for item in pairedlist:
            if max(item[0]) == str1_len-1 or min(item[1]) == 0:
                score3end += len(item[0])
                worst_list.append(item)
        return score3end, worst_list
    
    def draw_structure(self, pairedlist, str1, str2):
        """plain text representation of structures"""
        #e.g. [(set([1, 2, 3, 4]), set([16, 17, 18, 15])), ..., ]
        s = ''
        if len(pairedlist)>0:
            for apair in pairedlist:
                sense, antisense = list(apair[0]), list(apair[1])
                sense.sort()
                antisense.sort()
                if sense[0] > antisense[0]:
                    s += "5'- " + str1 + " -3'" + '\n'
                    s += ' '*(sense[0]+4) + '|'*len(sense) + '\n'
                    s += ' '*(sense[0]-antisense[0]) + "3'- " + self.reverse(str2) + " -5'" + '\n\n'
                else:
                    s += ' '*(antisense[0]-sense[0]) +"5'- " + str1 + " -3'" + '\n'
                    s += ' '*(antisense[0]+4) + '|'*len(sense) + '\n'
                    s += "3'- " + self.reverse(str2) + " -5'" + '\n\n'
        else:
            s += 'No dimer!\n'
        return s
    
    def hairpin(self, seq):
        pairedlist = self.list_structure(seq, seq)
        seq_len = len(seq)
        hairpinlist = []
        #at least 3 bases separation!
        for item in pairedlist:
            a, b = item[0], set([seq_len-1-x for x in item[1]])
            if max(b)-min(a) > 7:
                while min(b)-max(a) < 4:
                    b.remove(min(b))
                    a.remove(max(a))
                hairpinlist.append((a, b))
        return hairpinlist
    def score_hairpin(self, hlist):
        #because self-dimers already count part of hairpin effect, half penalty here
        score_longest = 0
        if len(hlist) > 0:
            score_longest += max([len(item[0]) for item in hlist])
        return score_longest/2
    def draw_hairpin(self, hlist, seq):
        s = ''
        seq_len = len(seq)
        if len(hlist)>0:
            for apair in hlist:
                a, b = apair
                linker = (min(b)-max(a))-1
                linkera = linker/2
                linkerb = linker - linkera
                head, tail = min(a), seq_len-max(b)-1
                breaker = max(a)+linkera+1
                if head > tail:
                    s += "5'- " + seq[:breaker] +' '*(linkerb-linkera) + '+\n'
                    s += ' '*(4+head) + '|'*len(a) + ' '*linkerb + '+\n'
                    s += "3'- " + ' '*(head-tail) + self.reverse(seq[breaker:]) + '+\n\n'
                else:
                    s += "5'- " + ' '*(tail-head) + seq[:breaker] +' '*(linkerb-linkera) + '+\n'
                    s += ' '*(4+tail) + '|'*len(a) + ' '*linkerb + '+\n'
                    s += "3'- " + self.reverse(seq[breaker:]) + '+\n\n'
        else:
            s += 'No hairpin!\n'
        return s



class SinglePrimer(PrimerFunctions):
    def __init__(self, seq):
        self.seq = seq
    def compute(self):
        self.hairpin = self.hairpin(self.seq)
        self.selfdimer = self.list_structure(self.seq, self.seq)
        score3end1, sc, self.worst = self.score(self.selfdimer, len(self.seq))
        self.total_score = score3end1 + sc + self.score_hairpin(self.hairpin)
    def display(self):
        self.header = 'input_sequence: ' + self.seq + ' total_score: ' + str(self.total_score)
        self.structures = self.draw_hairpin(self.hairpin, self.seq) + self.draw_structure(self.worst, self.seq, self.seq)

class PrimerPair(PrimerFunctions):
    def __init__(self, str1, str2):
        self.seq1, self.seq2 = str1, str2
        #self.compute()
    def compute(self):
        self.hairpin1 = self.hairpin(self.seq1)
        self.selfdimer1 = self.list_structure(self.seq1, self.seq1)
        score3end1, sc1, self.worst1 = self.score(self.selfdimer1, len(self.seq1))
        self.hairpin2 = self.hairpin(self.seq2)
        self.selfdimer2 = self.list_structure(self.seq2, self.seq2)
        score3end2, sc2, self.worst2 = self.score(self.selfdimer2, len(self.seq2))
        self.crossdimer = self.list_structure(self.seq1, self.seq2)
        score3endc, sc, self.worstc = self.score(self.crossdimer, len(self.seq1))
        self.total_score = score3end1+score3end2+score3endc + max([sc1, sc2, sc]) + self.score_hairpin(self.hairpin1) + self.score_hairpin(self.hairpin2)
    def display(self):
        self.header = 'input_sequences: ' + self.seq1 +', '+ self.seq2 + ' total_score: ' + str(self.total_score)
        self.structures = 'sense primer:\n' + self.draw_hairpin(self.hairpin1, self.seq1) +  self.draw_structure(self.worst1, self.seq1, self.seq1) + 'antisense primer:\n' + self.draw_hairpin(self.hairpin2, self.seq2) + self.draw_structure(self.worst2, self.seq2, self.seq2) + 'Both primers:\n' + self.draw_structure(self.worstc, self.seq1, self.seq2)



if __name__ == '__main__':
    #test
    b = SequenceUtil()
    """
    t = b.degenerate_check('TTAAGCCATGCAAGTGTAAGTACA', 'RTTCCCAGAGCCAGTABCCA')
    b.html_write(t, 'test1.html')
    """
    t = b.degenerate_check_single('TTAAGCCATGCAAGTGTAAGTACA')
    print t