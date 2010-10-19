#! usr/bin/env python
# Script file to pipeline the orthomcl files. I will run two arguments (1) for stage 1 which can be done in either of the servers. Preferably in einstein or lovelace as I will need to run blastp on this. (2) for stage 2 - it needs to be run on noether as i have mysql database access created in there. 
# I am creating a database orthomcl100 to say i am doing it for 100 orgs. 
# Also I need to run Step 4 - orthomclInstallSchema on noether. This is important as it sets up the database to which we will need to access. So this must be done prior to running this script. This requires a database to be setup which will be orthomcl100 in this case. 
# The taxon_ids are the same as that given in the parsinggbk.py. Note that the numbering for taxon_number which will become taxon_id is the same as that of the list file that is passed to parsinggbk.  

import sys
import os
import subprocess



def open_file(file_name,open_status = 'r'):
    """ This function just opens the file for reading and returns back the file handle for the file. If the file cannot be opened it just exits! It can open and read any type of files. It can also return back the file handle for writing. The default status for opening is to read the file. Note that the file handle is and must be closed in the place from where this function was called
    """
    
    #Opening and reading/writing the passed file_name """
    try:
	file = open(file_name,open_status)
    except IOError:	# If the file cannot be opened i am exiting
	print "File %s cannot be opened. Exiting "%file_name
	sys.exit(0)
	
    return file
	

def stage_one():
    """ This pipes the orthomcl programme from Stage 1 to stage 9 which includes doing an allvsall blast and parsing the blast results. This results is then taken to the next stage. Note here - Formatdb must be run only after the faa files are copied to the faafiles directory which will be used for sequence_id and blast etc """
    
    
    # Initializing
    orthomcl_bin_directory = '/home/jaggu/documents/orthomcl/bin'
    orthomcl_results_directory = '/home/jaggu/documents/orthomcl/results'
    orthomcl_faafiles_directory = '/home/jaggu/documents/orthomcl/faafiles'
    orthomcl_complaint_files = '/home/jaggu/documents/orthomcl/complaint_files'
    cwd = os.getcwd()
    status = 'Incomplete'
    
    #Step 1 - 4 is inconsequential. It is installation of mcl, mysql etc. Note as said before Step 4 must be run prior to executing the programme.
    
    #Step 5 - orthomclAdjustFasta
    # The main crux of this script will be to identify the locations of the different .faa files parsed by parsinggbk.py. The perl programme actually does the fasta file in proper format but which I have already done. It seems important that the faa file locations are identified and according to the list file another parsed_list file be made which has two columns or delimiters (1) taxon_number (0..) and (2) destination of the parsed file. The ordering must match between the initial listfile and the destination list file. For this sake i created the destination file in the parsinggbk script. Note that I want to keep ecoli ref MG1355 as the first line or taxon_id 0. So in my listfile i will always change the order manually to keep the ecoli-ACCNO: NC_000913
    # So I will just have to open the destination_list.lst file which contains the taxon_id|destination_path. Thus there is no requirement to copy all these files to the faafiles directory. 
    # The perl script is in this format - perl orthomclAdjustFasta taxon_code fasta_file id_field.
    # The fasta_file and the taxon_code can be obtained just by opening the destination_list.lst file. I found out that long accno is also ok for a taxon code. For all the trouble i took to make. But nevertheless the destination_list.lst is important. It is needed and I can get the accno for the taxon_code as well. The id_field that will be used throughout will be the locus_tag which is in the 7th position. 
    
    ifile = open_file('destination_list.lst','r')
    lines = ifile.readlines()
    ifile.close()
    os.chdir(orthomcl_complaint_files)
    
    for line in lines:
	split_line = line.split('|')
	taxon_id = (split_line[1].split('/'))[7]
	commandline = 'perl '+orthomcl_bin_directory + '/orthomclAdjustFasta ' + taxon_id + ' ' + split_line[1][:-1] + ' 7'
	os.system(commandline)
    
    os.chdir(cwd)
    print "Step 5 - orthomclAdjustFasta : Completed"
	
    
    #Step 6 - orthomclFilterFasta
    # This perl script will handle the Fastafiles from step 5 created in orthomcl_complaint_fasta_files directory. This step filters all the fasta files in the complaint directory and makes a result file - containing all the good proteins and all the bad proteins (in separate file). The command is just orthomclFilterFasta input_dir min_length max_percent_stops. Generally i dont get any sequences in the bad proteins. 
    
    os.chdir(orthomcl_results_directory)
    commandline = 'perl '+ orthomcl_bin_directory + '/orthomclFilterFasta ' + orthomcl_complaint_files + ' 10 20 '
    print commandline
    os.system(commandline)
    
    os.chdir(cwd)
    
    print "Step 6 - orthomclFilterFasta : Completed"
    
    
    #Step 7 - allvsall blast. This is the big blast which may take quite some time and so mandatory run on lovelace and einstein. Note that you need to do a formatdb and then a blastp and in the form -m8. This format is critical and take pains to ensure that it is -m8 or it is a lot of time wasted. 
    
    os.chdir(orthomcl_results_directory)
    
    inputfile = 'goodProteins.fasta'
    logfile = 'formatdb.log'
    commandline = 'formatdb -i ' + inputfile + ' -oF ' + '-pT -l ' + logfile
    os.system(commandline)
    
    evalue = '0.00001' # evalue (set at) = E-5
    result_file_name = 'allvsall.bst'
    
    database = 'goodProteins.fasta'
    query = 'goodProteins.fasta'
    commandline = 'blastall -p blastp -d ' + database + ' -i ' + query + ' -e ' + evalue + ' -m 8 ' + ' -o ' + result_file_name
    print "Blasting...."
    os.system(commandline)
    
        
    print "Step 7 - AllvsAll Blast : Completed"
    
    #Step 8 - orthomclBlastParser. This is the step where the blast file result is parsed. It also calculates the percent match of everyhit (not so much needed at this stage i think). This is the last step that can be run on the faster computers. 
    
    commandline = 'perl '+ orthomcl_bin_directory + '/orthomclBlastParser ' + result_file_name + ' ' + orthomcl_complaint_files + ' >> similarSequences.txt'
    os.system(commandline)
    
    print "Step 8 - orthomclBlastParser : Completed"
    
    status = 'Successful'
           
    return status
    

def stage_two():
    """ This function carries out the stage_two of the entire orthomcl operation wherein access to mysql database is needed. Thus for the time being atleast it must be done in noether only. Probably it will take time, but then no option. 
    """
    
    # Initializing
    orthomcl_bin_directory = '/home/jaggu/documents/orthomcl/bin'
    orthomcl_results_directory = '/home/jaggu/documents/orthomcl/results'
    orthomcl_faafiles_directory = '/home/jaggu/documents/orthomcl/faafiles'
    orthomcl_complaint_files = '/home/jaggu/documents/orthomcl/complaint_files'
    cwd = os.getcwd()
    status = 'Incomplete'
    
    
    #Step 9 - orthomclLoadBlast. This loads the results onto the mysql database orthomcl100
    # EXAMPLE: orthomclSoftware/bin/orthomclLoadBlast my_orthomcl_dir/orthomcl.config my_orthomcl_dir/similarSequences.txt
    
    os.chdir(orthomcl_results_directory)
    commandline = 'perl '+ orthomcl_bin_directory + '/orthomclLoadBlast ' + 'orthomcl.config ' + 'similarSequences.txt'
    os.system(commandline)
    print "Step 9 - orthomclLoadBlast : Completed"

    #Step 10 - orthomclPairs. Well this is computationally intensive step and it makes the kind of pairs. There are many options like cleanup options or options to start over once it has stalled. But for now i think it is clean and will do without those options, just the default
    #EXAMPLE: orthomclSoftware/bin/orthomclPairs my_orthomcl_dir/orthomcl.config my_orthomcl_dir/orthomcl_pairs.log cleanup=no
 
    commandline = 'perl '+ orthomcl_bin_directory + '/orthomclPairs ' + 'orthomcl.config  orthomcl_pairs.log cleanup=no'
    os.system(commandline)
    print commandline
    print "Step 10 - orthomclPairs : Completed"
    
    # Step 11 - orthomclDumpPairsFile. This step makes a directory pairs and prints the files of orthologs, paralogs etc as txt files. 
    # EXAMPLE: orthomclSoftware/bin/orthomclDumpPairsFile my_orthomcl_dir/orthomcl.config 
    
    commandline = 'perl '+ orthomcl_bin_directory + '/orthomclDumpPairsFiles ' + 'orthomcl.config'
    os.system(commandline)
    print commandline
    print "Step 11 - orthomclDumpPairsFile : Completed"
    
    #Step 11 - The major programme - mcl. This too needs to be run in noether as none of the other machines have it installed. 
    # EXAMPLE: mcl my_orthomcl_dir/mclInput --abc -I 1.5 -o my_orthomcl_dir/mclOutput

    commandline = 'mcl mclInput --abc -I 1.5 -o mclOutput'
    os.system(commandline)
    print "Step 12 - mcl programme execution : Completed"
    
    #Step 12 - OrthoMclToGroups. This is the final step and makes a groups.txt file, the most important file for this project !
    # EXAMPLE: mclOutput2groupsFile prefix starting_id_num
    #  orthomclMclToGroups my_prefix 1000 < mclOutput > groups.txt
    # I am going to keep the prefix as EM_ and the starting_id_num from 1. 
    
    commandline = orthomcl_bin_directory + '/orthomclMclToGroups EM_ 1 < mclOutput > groups.txt'
    os.system(commandline)
    print "Step 13 - orthomcltogroups : Completed"
    
    status = 'Successful'
    
    return status
        
    

try:
    argument = sys.argv[1]
except IndexError:
    print "No argument passed. Exiting.."
    sys.exit(0)
    
if argument == '1': 
    status = stage_one()
    print status
    
elif argument == '2':
    status = stage_two()
    print status
    
else: 
    print "Invalid argument passed. Exiting.."
    sys.exit(0)
	
