#!/usr/local/bin/python
import hybes


#Define 2 groups here

groups=[
    {'VR2332':[19,
               13,
               3,
                115,
               118,
               ]},


    {'SDSU73':[9,
               53,
               98,
               108,
               110,
               107,
               113,
               ]}
    ]

columnHeaders = ('Oligo',
                 'A_mean', 'A_std', 'A_n',
                 'B_mean', 'B_std', 'B_n',
                 't','p')

a=hybes.IntensityArray('cluster532.txt')
a.maskByIntenisty(0)
a.maskNonViral()
a.sumNormalize()
stats=a.arrayGroupCompare(groups[0].values()[0],
                          groups[1].values()[0],
                          maxP=0.1)

print a.clusterPath
print ("Group A (%s) Arrays: %s" % (groups[0].keys()[0],
       ','.join([str(x) for x in sorted(groups[0].values()[0])])))
print ("Group B (%s) Arrays: %s" % (groups[1].keys()[0],
       ','.join([str(x) for x in sorted(groups[1].values()[0])])))

print "Data Manipulations:"
for op in a.operations:
    print "    %s" % op

print
print "\t".join(columnHeaders)

for row in stats:
    print "\t".join([str(x) for x in row])
