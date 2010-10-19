#log# Automatic Logger file. *** THIS MUST BE THE FIRST LINE ***
#log# DO NOT CHANGE THIS LINE OR THE TWO BELOW
#log# opts = Struct({'__allownew': True, 'logfile': 'ipython_log.py'})
#log# args = []
#log# It is safe to make manual edits below here.
#log#-----------------------------------------------------------------------
_ip.system("ls -F ")
more ecoli_data.txt
_ip.system("cat ecoli_data.txt")
ifile = open('ecoli_data.txt','r')
lines = ifile.readlines()
lines[0]
lines[0].split('|')
for line in lines:
    line = line.split('|')
    operon_name = line[0]
    locus_tag = line[1] + '|' + line[2]
line
operon_name
abc = ['a','b','c']
xyz = ['x','y','z']
something = []
something.append(abc)
something
something.append(xyz)
something
for some in something:
    print some

group_name_a = 'JEM_7707'
output_contents = """SELECT c.sequence_id, c.rank, c.left_end, c.right_end, g.group_id
    FROM
    (SELECT group_id, sequence_id
    FROM groups
    WHERE group_id = '%s') g
    INNER JOIN coding_sequences c
    ON c.sequence_id = g.sequence_id
    """%(group_name_a)
output_contents
astring = ''
if not astring:
    print "yes"
import pylab
nos = [0,0,0,1,2,333,4]
hist(nos,5)
from pylab import *
hist(nos,5)
_ip.magic("logstart ")

from scipy import *
a = zeros(1000)
a[:100] = 1
b = fft(a)
from pylab import *
plot(abs(b))
show()
quit()
