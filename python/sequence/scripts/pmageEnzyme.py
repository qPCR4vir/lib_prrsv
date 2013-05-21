#!/usr/local/bin/python

import sys
import re
from optparse import OptionParser

from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg

from ncbi import genbank, giInfo
from sequence import NAambig

def insert_in_dict_with_list(dictionary, key, value):
    if not dictionary.has_key(key):
        dictionary[key] = list()
    dictionary[key].append(value)

class Analizer(object):
    def __init__(self, enzyme, record, tag_length, originalEnzyme=None):
        self.record = record
        self._original_tag_length = tag_length
        if originalEnzyme is None:
            self.originalEnzyme = enzyme
        else:
            self.originalEnzyme = originalEnzyme
        self.enzyme = enzyme
        self.tag_length = tag_length - len(self.originalEnzyme)

        self.tags = dict()
        self._genes_c = 0

    def _get_feature_id(self, feature):
        if feature.qualifiers.has_key('gene'):
            return feature.qualifiers['gene']
        return feature.qualifiers['protein_id']

    def process(self):
        for feature in self.record.features():
            if feature.type != 'CDS':
                continue
            if len(feature.location.regions) != 1:
                continue

            region = feature.location.regions[0]
            start = region.start - 1
            end = region.end
            feat_id = self._get_feature_id(feature)
            sequence = self.record.sequence[start:end].upper()
            if region.complement:
                sequence = genbank.reverseComplement(sequence)
            if not sequence.startswith('ATG'):
                msg = ('Feature %s don\'t start with ATG. It starts '
                       'with %s') % (str(feat_id), sequence[0:3])
                print >> sys.stderr, msg
            if type(self.enzyme) is str:
                parts = sequence.split(self.enzyme)
            else:
                parts = self.enzyme.split(sequence)
                
            if len(parts) == 1:
                insert_in_dict_with_list(self.tags, '', feat_id)
            else:
                part = parts[-1]
                key = self.originalEnzyme + part[0:self.tag_length]
                insert_in_dict_with_list(self.tags, key, feat_id)
            self._genes_c += 1

    def getUniqueGenes(self):
        result = dict()
        for tag, genes in self.tags.iteritems():
            if len(genes) == 1:
                insert_in_dict_with_list(result, genes[0], tag)
        return result

    def saveToFile(self):
        results = list()
        if self.tags.has_key(''):
            results.append(('', self.tags['']))
        for tag, genes in self.tags.iteritems():
            if tag == '':
                continue
            results.append((tag, genes))
        
        fout = file('enzyme-%s-length-%i.txt' % (self.originalEnzyme,
                                                 self._original_tag_length),
                    'w')
        for tag, genes in results:
            genes_output = '\t'.join(genes)
            print >> fout, '%s\t%s' % (tag, genes_output)
        fout.close()
        
    def generateChart(self):
        u_genes = self.getUniqueGenes()
        data = dict()
        for gene, tags in u_genes.iteritems():
            if not data.has_key(len(tags)):
                data[len(tags)] = 0
            data[len(tags)] += 1
        data[0] = self._genes_c - len(u_genes.keys())

        xs = list()
        ys = list()
        # Convert the values to %
        for k, v in data.iteritems():
            xs.append(k)
            ys.append((float(v) / self._genes_c))

        fig = Figure()
        ax = fig.add_subplot(111)
        ax.bar(xs, ys, width=0.5, align='center')
        fig.get_axes()[0].set_ylabel('% of genes')
        fig.get_axes()[0].set_xlabel('# of unique tags')
        #fig.get_axes()[0].set_yscale('log')
        canvas = FigureCanvasAgg(fig)
        canvas.print_figure('enzyme-%s-length-%i.png' % \
                            (self.enzyme, self._original_tag_length),
                            dpi=96)
        return data


class EnzymeFileReader(object):
    def __init__(self, fileobject):
        self.input = fileobject
        self._enzymes_read = set()

    def _process_enzyme(self, enzyme):
        have_parenthesis = enzyme.find('(') != -1 or \
                           enzyme.find(')') != -1
        if have_parenthesis:
            return None
        elif enzyme.find('^') != -1:
            enzyme = enzyme.replace('^', '')
        enzyme_re = enzyme
        for na, replacements in NAambig.iteritems():
            repl_re = '(?:%s)' % '|'.join(replacements)
            enzyme_re = re.sub(na, repl_re, enzyme_re)
        if enzyme_re != enzyme:
            enzyme_re = '(?:%s)' % (enzyme_re,)
            return (re.compile(enzyme_re), enzyme)
        return (enzyme, None)

    def enzymes(self):
        for line in self.input:
            stripped_line = line.strip()
            if len(stripped_line) == 0 or \
               stripped_line.startswith('REBASE') or \
               stripped_line.startswith('=-=') or \
               stripped_line.startswith('Copyright') or \
               stripped_line.startswith('Rich'):
                continue

            columns = line.split('\t')
            if len(columns[5]) == 0:
                continue
            if len(columns[1]) != 0:
                continue
            p = self._process_enzyme(columns[2])
            if p is None:
                continue
            enzyme, enzyme_original = p
            if enzyme in self._enzymes_read:
                continue
            self._enzymes_read.add(enzyme)
            name = columns[0]
            yield name, enzyme, enzyme_original


def main(args):
    parser = OptionParser('usage: %prog [options]')
    parser.add_option('-g', '--gi', dest='gi', type='int',
                      help='gi of the genbank record to analyze.')
    parser.add_option('-r', '--record', dest='record', type='string',
                      help='read genbank record from FILE',
                      metavar='FILE')
    parser.add_option('-e', '--enzyme', dest='enzyme', type='string',
                      help='enzyme to process')
    parser.add_option('-l', '--list', dest='enzymefile',
                      type='string', help='file of enzymes to process',
                      metavar='FILE')
    parser.add_option('-s', '--size', dest='size', type='int',
                      help='size of the tag')
    (options, args) = parser.parse_args(args)

    if options.size:
        length = options.size
    else:
        print >> sys.stderr, 'Usage error: must specify a length'
        print >> sys.stderr, parser.format_help()
        sys.exit(1)

    if options.record:
        record = genbank.Record(file(options.record))
    elif options.gi:
        record = giInfo.GiRecord(options.gi, True)
    else:
        print >> sys.stderr, 'Usage error: a record file (-r) or ' + \
              'a gi (-g) is required'
        print >> sys.stderr, parser.format_help()
        sys.exit(1)

    if options.enzymefile:
        enzyme_reader = EnzymeFileReader(file(options.enzymefile))
        for name, enzyme, originalEnzyme in enzyme_reader.enzymes():
            if len(originalEnzyme or enzyme) >= length:
                continue
            analizer = Analizer(enzyme, record, length, originalEnzyme=originalEnzyme)
            analizer.process()
            analizer.saveToFile()
    elif options.enzyme:
        if length <= len(options.enzyme):
            msg = 'Length of the tag cannot be smaller than the enzyme.'
            print >> sys.stderr, msg
            sys.exit(1)
        analizer = Analizer(options.enzyme, record, length)
        analizer.process()
        analizer.saveToFile()
    else:
        print >> sys.stderr, 'Usage error: an enzyme (-e) or ' + \
              'a file with enzymes (-l) is required'
        print >> sys.stderr, parser.format_help()
        sys.exit(1)


if __name__ == '__main__':
    main(sys.argv)
