# ---------------------------------------------------
# --                                               --
# --    Hairpin Analyzer v0.3                      --
# --                                               --
# --    --> Main file                              --
# --                                               --
# --                                               --
# --  Developed during 2011-Aug-01 to 2011-Aug-23  --
# --  by Mathias Bader (mail@mathiasbader.de)      --
# --  Universitaet des Saarlandes                  --
# ---------------------------------------------------
#

import os, imp, sys
from glob import glob
from time import localtime, strftime

# import classes, configuration and specification
HairpinClasses = imp.load_source('HairpinData', 'classes/HairpinClasses.py')
Configuration = imp.load_source('Configuration', 'configuration/configuration.py')
Specification = imp.load_source('Specification', 'classes/HairpinSpec.py')

# output some general information about the software
line_software_outline = "------------------------------------------------------------"
line_software_empty   = "--                                                        --"
line_software_inf     = "--              Hairpin Analyzer v" + Specification.software_version
for i in range(24 - len(Specification.software_version)):
    line_software_inf += " "
line_software_inf += "--"
print
print line_software_outline
print line_software_inf
print line_software_outline

# function to print a line to stdout with boxed style
def print_line(line):
    sys.stdout.write("--  ")
    sys.stdout.write(line)
    for i in range(52-len(line)):
        sys.stdout.write(" ")
    sys.stdout.write("  --\n")
    sys.stdout.flush()

# print the current time
print_line("                                            " + strftime("%H:%M:%S", localtime()))

# import configuration
import imp
Configuration = imp.load_source('Configuration', 'configuration/configuration.py')
filename_datafile    = Configuration.filename_results
filename_summaryfile = Configuration.filename_results_summary
main_folder = Configuration.data_path_main_folder



# search for files that contain data for hairpin
path_to_data = main_folder + Configuration.data_path_cpg
filelist = glob(path_to_data + "*/*/" + filename_datafile)


# if no files were found in the specified folder
if (len(filelist) == 0):
    print line_software_outline
    print_line("Couldn't find any data in the specified data folder.")
    print_line("Please open the file configuration/configuration.py")
    print_line("to specify the data path.")
    print line_software_empty
    print_line("The folder you specified is as follows:")
    print_line(path_to_data)
    print line_software_empty
    print line_software_outline
    print
else:
    # calculate and output progress bar length
    prog_bar_len = 44
    if len(filelist) == 1:
        print_line("One file found in ")
        print_line(path_to_data)
    else:
        print_line(str(len(filelist)) + " files found in")
        print_line(path_to_data)
    print_line("")
    sys.stdout.write("--  0%|")
    sys.stdout.flush()
    prog_bar_piece_len = prog_bar_len / len(filelist)
    prog_bar_piece = ""
    for i in range(prog_bar_piece_len):
        prog_bar_piece += "-"
    
    # figure out whether slashes or backslashes are used (Windows <-> Unix)
    seperator = "\\"	 # windows system
    if (filelist[0].find(seperator) == -1):
        seperator = "/"  # unix system

    # sort the filelist alphabetically
    filelist.sort()

    # figure out longest filename for alignment of output to stdout
    max_filename_length = 0
    for filename in filelist:
        if len(filename) > max_filename_length:
            max_filename_length = len(filename)

    # dictionary for collecting summary information
    # seperatly for each amplicon type
    summary_info = {}
    
    cur_file_count = 0
    processed_file_list = []
    for filepath in filelist:
        cur_file_count += 1
        folders = filepath[len(path_to_data):len(filepath)-len(filename_datafile)-1]
        folders = folders.split(seperator)
        sub_path = folders[0] + "/" + folders[1] + "/"

        # ----------------------------------------------------------
        # -- This command does the main work by calling the       --
        # -- according methods from the Hairpin class: reading    --
        # -- file, cutting specified positions, mapping data,     --
        # -- cutting again, reading linker data, sort data and    --
        # -- store all information in appropriate data variables. --
        # ----------------------------------------------------------
        hp = HairpinClasses.HairpinData(sub_path)


        # figure out whether there is additional data for this amplicon
        add_file = {}
        add_file["linker"]  = hp.get_additional_data_available("linker")
        add_file["non_cpg"] = hp.get_additional_data_available("non_cpg")
        add_file["snp"]     = hp.get_additional_data_available("snp")

        # remember processed file for later output to stdout
        new_line = " "
        if (add_file["linker"]):
            new_line += "l "
        else:
            new_line += "  "
        if (add_file["non_cpg"]):
            new_line += "n "
        else:
            new_line += "  "
        if (add_file["snp"]):
            new_line += "s "
        else:
            new_line += "  "
        new_line += "  "
        new_line += sub_path + " "
        # for aligned output
        for i in range(max_filename_length - len(filepath)):
            new_line += " "

        if (hp.result_row_count == 0):
            # there are no lines in this file
            new_line += " [no matching ids]"
        else:

            if add_file["linker"]:
                # get the statistics for matching with the linker file
                linker_match = hp.get_linker_match_information()
                linker_match_sum = linker_match[0] + linker_match[1] + linker_match[2] + linker_match[3]
                new_line += " "
                if (linker_match[0] < 10):
                    new_line += " "
                if (linker_match[0] < 100):
                    new_line += " "
                if (linker_match[0] < 1000):
                    new_line += " "
                new_line += str(linker_match[0]) + " "
                if (linker_match[1] < 10):
                    new_line += " "
                new_line += str(linker_match[1]) + " "
                if (linker_match[2] < 10):
                    new_line += " "
                new_line += str(linker_match[2]) + " "
                if (linker_match[3] < 10):
                    new_line += " "
                if (linker_match[3] < 100):
                    new_line += " "
                if (linker_match[3] < 1000):
                    new_line += " "
                new_line += str(linker_match[3])
            
            # print the results to std out (out-comment the following line if necessary)
            #hp.print_results()
            
            # write results to textfile
            hp.write_results_to_file('results/' + folders[0] + '_' + folders[1] + '_results.txt')

            # create heatmap as png-image
            hp.create_image('results/' + folders[0] + '_' + folders[1] + '_results.png')

            # collect summary information about this amplicon
            amplicon_type = hp.get_amplicon_type()
            if (not amplicon_type in summary_info):
                # if this is the first data line, include meta data line before
                summary_info[amplicon_type] = [hp.get_amplicon_meta_information(),
                                               hp.get_amplicon_information()]
            else:
                # add data line
                summary_info[amplicon_type].append(hp.get_amplicon_information())

        # delete the hairpin-object since we don't need
        # it anymore and therefore free the memory
        del(hp)

        # print the progressbar
        sys.stdout.write(prog_bar_piece)
        if (cur_file_count <= prog_bar_len % len(filelist)):
            sys.stdout.write("-")
        sys.stdout.flush()
        processed_file_list.append(new_line)
        
    # write a summary file for each amplicon type
    summary_filename       = Configuration.results_summary_filename
    summary_file_extension = Configuration.results_summary_file_extension
    for amplicon in summary_info:
        results_summary_filename = "results/" + summary_filename + amplicon + summary_file_extension
        try:
            file_summary = open(results_summary_filename, 'w')
        except IOError as e:
            print 
            print_line("Error: Can't write to file '" + results_summary_filename + "")
            assert(False)

        # write summary information to file
        for line in summary_info[amplicon]:
            file_summary.write(line)
        file_summary.close()
    

    # output final information to stdout
    print "|100%  --"
    print_line("")
    print_line("The file " + filename_datafile)
    print_line("has been analyzed in the following folders:")
    print_line("")
    print_line("            (data <-> linker) matched lines")
    print_line(" linker file found            |  under threshold")
    print_line(" | file with non CpGs found   |  |  only data line")
    print_line(" | | file with SNPs found     |  |  |    only linker")
    print_line(" | | |                        |  |  |    |")
    for line in processed_file_list:
        print_line(line)
    print_line("")
    print_line("Results have been written to the folder 'results/'")
    print_line("")
    print line_software_outline
    print
