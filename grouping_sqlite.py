#! usr/bin/python


import sqlite3


   
conn = sqlite3.connect('orthomcl.db')
cursor = conn.cursor()

table_name = 'groups'
contents = """ CREATE TABLE %s
(
sequence_id VARCHAR(30),
group_id VARCHAR(8)
)
"""%(table_name)

#cursor.execute(contents)

ifile = open('groups.txt','r')

lines = ifile.readlines()
for line in lines:
    split_line = line.split(' ')	#Space is all that separates the different sequenceids
    group_name = (split_line.pop(0))[:-1]	#It has a ':' in the end
    for sequences in split_line:
	sequences = sequences.replace('\n','')	#Some seq id has this annoying carriage return
	
	#Inserting this information into the table
	insert_table = """ INSERT INTO groups (sequence_id, group_id) VALUES ('%s','%s')
	"""%(sequences,group_name)
	cursor.execute(insert_table)
	print cursor
	print contents

ifile.close()
#insert_table = """  INSERT INTO groups (sequence_id, group_id) VALUES ('NC_011740|EFER_4020','JEM_54831') """

cursor.execute(insert_table)



cursor.close()
conn.commit()
conn.close()
