#! /usr/bin/python
# This is a trial python script to check if blastall is working fine and how long does it take to execute one blastall

import os

def blastall_two():
    
        
    
    command_line_1 = 'blastall -p blastp -d /project/marcotte/jagannath/projectfiles/bacterialoperons/parsed_gbk_files_genomes/NC_007925/NC_007925.fasta -i /project/marcotte/jagannath/projectfiles/bacterialoperons/parsed_gbk_files_genomes/NC_007925/NC_007925.fasta -e 0.001 -m 8  -o /project/marcotte/jagannath/projectfiles/bacterialoperons/blast_results2/NC_007925vsNC_007925.blast'
    
    command_line_2 = 'echo /project/marcotte/jagannath/projectfiles/bacterialoperons/blast_results/NC_007925vsNC_007925.blast'
    
    print command_line_1
    os.system(command_line_1)
    
    os.system(command_line_2)
    
    return

if __name__ == '__main__':
    #blast_sh_file = os.getcwd() + '/blast_sh/blastall_1.sh'

    blastall_two()
    from timeit import Timer
    t = Timer(stmt = "blastall_two()",setup = "from __main__ import blastall_two",)
    print t.timeit(0)
#    print strftime("%d %b %Y %H:%M:%S",localtime())
    