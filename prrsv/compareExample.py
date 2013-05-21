#!/usr/local/bin/python
import hybes

outputFile=sys.stdout
# send the data to a file like so:
# outputFile=file("myexample.txt","w")


#Define 2 groups here
# a somewhat complicated data structure, but it has the advantage
# of keeping the group labels and the microarray lists in one place.
#
# the "groups" variable is a list (python data type).
# Lists are ordered collections of data.  Groups has 2 items, 
# one for each group of microarrays you want to compare.
# 
# Each of the items in groups is a dictionary (python data type).
# They are enclosed by {}.  So the overall structure of groups is:
# [ {}, {} ], a list with 2 dictionaries inside it. 
#
# Dictionaries are sets of values that are each associated with a key.  
# In this groups variable each dictionary has a single key:value pair.
# The key is the name of the set of arrays, and the value is a list 
# containing integers that specify the set of arrays
#
# The sets and names can be changed as desired.  Pay special attention 
# to making sure that you keep commas seperating the integers and the 
# dictionaries, and keep the keys surrounded by quotes. 
#
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
# this reads the data in the same format that culster uses for input.
# the varible that holds the data is 'a'. 
a=hybes.IntensityArray('cluster532.txt')  # IntensityArray is defined in
                                          # the hybes.py file, so it is 
                                          # referneced hybes.Intensity.
                                          # n.b. hybes was imported above 

#
# Examples of data transformations.
# the variable 'a' represents an 'instance' of the 'class' IntensityArray.
# class instances are a nice place to put data and also pirces of code
# to operated on the type of data in the class.
#
# IntensityArray instances (like 'a') hold intensities for 1 or more arrays 
# in a matrix.  And they know how to mask and normaize in various ways.
# e.g.
a.maskByIntenisty(0)
a.maskNonViral()
a.sumNormalize()

# This next part compares the oligo responses in the 2 groups
stats=a.arrayGroupCompare(groups[0].values()[0],
                          groups[1].values()[0],
                          maxP=0.1)

# print out the history of what we have done, the input file etc.
print >> outputFile,  a.clusterPath
print >> outputFile,  ("Group A (%s) Arrays: %s" % (groups[0].keys()[0],
       ','.join([str(x) for x in sorted(groups[0].values()[0])])))
print >> outputFile,  ("Group B (%s) Arrays: %s" % (groups[1].keys()[0],
       ','.join([str(x) for x in sorted(groups[1].values()[0])])))

print >> outputFile,  "Data Manipulations:"
for op in a.operations:
    print >> outputFile,  "    %s" % op

print >> outputFile, 

# useful headers for each column
columnHeaders = ('Oligo',
                 'A_mean', 'A_std', 'A_n',
                 'B_mean', 'B_std', 'B_n',
                 't','p')
                 
# how to print a tab seperated line of data
print >> outputFile,  "\t".join(columnHeaders)

# note trick below to convert to strings, if fields include numbers

for row in stats:
    # this is the trick - it uses "list comprehension"
    print >> outputFile,  "\t".join([str(x) for x in row])

# done
