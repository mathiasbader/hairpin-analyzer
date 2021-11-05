# ---------------------------------------------------
# --                                               --
# --    Hairpin Analyzer v0.3                      --
# --                                               --
# --    --> class HairpinData                      --
# --                                               --
# --                                               --
# --  2011-Aug-01 to 2011-Aug-23                   --
# --  by Mathias Bader (mail@mathiasbader.de)      --
# --  at Saarland University                       --
# --  at Saarland University                       --
# --                                               --
# --  Zeile 592 von Julia geadded, wegen           --
# --  Mutationsanalyse - inkorrekte snps werden    --
# --  nicht mehr aussortiert                       --
# ---------------------------------------------------


import sys
import os

class HairpinData:
    """Class for handling of Hairpin data from one amplicon"""

    filepath = ""
    sub_path = ""
    file_paths_additional = []
    filename_datafile = ""
    filename_datafile_snp = ""
    filename_summary_file = ""
    amplicon_type = ""
    metadata = []
    data     = []
    row_count        = 0
    result_row_count = 0
    col_count        = 0
    mapped_col_count = 0
    results          = []
    results_unsorted = []
    results_available = False

    additional_data_available = {"linker": False, "non_cpg": False, "snp": False}

    data_linker_match           = 0
    data_linker_under_threshold = 0
    data_linker_data_only       = 0
    data_linker_linker_only     = 0
    data_snp_correct   = 0
    data_snp_incorrect = 0

    # these are the column numbers (starting at 1):
    columns_to_capture = [1,   # ID
                          3,   # Sequence_Identity
                          4,   # CG_Methylation_Pattern
                          # 5,   # Mean_CG_Methylation
                          7,   # Conversion
                          8,   # Reference
                          9    # Sample
                         ]
    # the following are not the same numbers but the indexs in the above list (starting at 0)
    data_col = 2
    conversion_col = 3
    col_amplicon_type = 4
    col_sample = 5
    columns_to_capture.sort()

    # these are the column numbers (starting at 1) of the additional files:
    columns_to_capture_add = {"linker": [
        1,  # ID
        4,  # CG_Methylation_Pattern
        5,  # Mean_CG_Methylation
        7  # Conversion
    ], "non_cpg": [
        1,  # ID
        4  # CG_Methylation_Pattern
    ]}

    # information about this amplicon which we want
    # to print into a file later
    info_sample = ""
    info_amplicon = ""
    info_reads = ""
    info_cpg = 0
    info_04 = 0
    info_03 = 0
    info_01 = 0
    info_00 = 0
    info_na = 0
    info_05 = 0
    info_06 = 0
    info_08 = 0
    info_09 = 0
    # additional statistics (not considering x and 0)
    info_stats_04       = 0  # fully methylated (only 4)
    info_stats_043_041  = 0  # fully and half methylated (4 and 3 or 4 and 1)
    info_stats_0        = 0  # no methylation (only 0)
    info_stats_01       = 0  # half methylated left side (only 1)
    info_stats_03       = 0  # half methylated right side (only 3)
    info_stats_013      = 0  # half methylated left and right side (1 and 3)
    info_stats_0134     = 0  # full, half left and half right methylated (4, 3 and 1)
    info_stats_only_one = 0  # exactly one methylaten on one strand and zero or one on the other (040, 030, 010 and 013)
    info_stats_mosaic   = 0  # mosaic pattern on one side (404, 101, 303, 104, 401, 304 and 403, but not 301 or 103 (number of zeros can be more than one))
    info_stats_cont     = 0  # continuous pattern without mosaic (0440, 0110, 0330, 0140, 0410, 0340 and 0430 but not 0130 or 0310 (number of follow-up methylations has to be at least 2))

    def __init__(self, sub_path):
        """Reads specified file and evaluates the data"""

        # initialize variables
        self.sub_path = ""
        self.filepath = ""
        self.file_paths_additional = []
        self.metadata = []
        self.data      = []
        self.row_count = 0
        self.result_row_count = 0
        self.col_count = 0
        self.mapped_col_count = 0
        self.results = []
        self.results_unsorted= []
        self.results_available = False

        self.additional_data_available = {"linker": False, "non_cpg": False, "snp": False}

        self.data_linker_match           = 0
        self.data_linker_under_threshold = 0
        self.data_linker_data_only       = 0
        self.data_linker_linker_only     = 0
        self.data_snp_correct   = 0
        self.data_snp_incorrect = 0

        # initialize statistic variables
        self.info_sample = ""
        self.info_amplicon = ""
        self.info_reads = ""
        self.info_cpg = 0
        self.info_04 = 0
        self.info_03 = 0
        self.info_01 = 0
        self.info_00 = 0
        self.info_na = 0
        self.info_05 = 0
        self.info_06 = 0
        self.info_08 = 0
        self.info_09 = 0
        # additional statistics (not considering x and 0)
        self.info_stats_04       = 0
        self.info_stats_043_041  = 0
        self.info_stats_0        = 0
        self.info_stats_01       = 0
        self.info_stats_03       = 0
        self.info_stats_013      = 0
        self.info_stats_0134     = 0
        self.info_stats_only_one = 0
        self.info_stats_mosaic   = 0
        self.info_stats_cont     = 0


        # import configuration
        import imp
        configuration = imp.load_source('Configuration', 'configuration/configuration.py')
        self.filename_datafile     = configuration.filename_results
        self.filename_datafile_snp = configuration.filename_results_snp
        self.filename_summary_file = configuration.filename_results_summary

        # read the specified file
        self.sub_path = sub_path
        self.filepath = configuration.data_path_cpg.strip()
        self.file_paths_additional = configuration.data_paths_additional
        main_folder = configuration.data_path_main_folder
        self.filepath = main_folder + self.filepath
        for i in range(len(self.file_paths_additional)):
            self.file_paths_additional[i][1] = main_folder + self.file_paths_additional[i][1]

        if not self.read_file():
            return None

        # define the amplicon type of this file
        self.define_amplicon_type()

        # erase specified positions before mapping
        if not self.erase_special_positions():
            print("Error in HairpinData.__init__(): column count is odd.")
            return None

        # map lines
        self.map_lines()

        # erase specified positions after mapping
        if not self.delete_mapped_positions():
            return None

        # bind additional data to results: linker data,
        # non CpG data, SPN data
        self.bind_additional_data_to_results()

        # delete lines where SNP do not match 100%
        self.delete_incorrect_snps()

        # data processing is finished - now calculate the statistics
        self.calculate_statistics()

        # sort the results
        self.sort_results()

        self.results_available = True


    def read_file(self):
        """reads data from the specified file into HairpinData object"""
        # check whether filename has been specified
        if (self.filepath + self.sub_path) == "":
            print("Error in class HairpinData.readfile(): No filepath specified")
            return False

        filename = self.filepath + self.sub_path + self.filename_datafile
        # read the specified file
        try:
            f = open(filename, 'r')
        except IOError as e:
            print("Error in class HairpinData.read_file(): The file '" + filename + "' could not be found.")
            return False

        first_line = True
        self.row_count = 0
        # the columns of the input file we want to work with
        for line in f:
            if line.strip() != "":
                if not first_line:
                    self.row_count += 1
                    line_splitted = line.split("\t")
                    # line must have as least as many columns as the highest number we want to read
                    if len(line_splitted) >= self.columns_to_capture[len(self.columns_to_capture) - 1]:
                        tmp_data = []
                        for col in self.columns_to_capture:
                            tmp_data.append(line_splitted[col-1].strip())
                        self.data.append(tmp_data)
                    else:
                        print("Error: The following line cannot be read:")
                        print(line)
                        print()
                        return False
                else:
                    first_line = False
                    line_splitted = line.split()
                    for col in self.columns_to_capture:
                        self.metadata.append(line_splitted[col-1])
        self.col_count = self.mapped_col_count = len(self.data[0][self.data_col])

        return True


    def define_amplicon_type(self):
        """Read amplicon type from data"""
        self.amplicon_type = self.data[0][self.col_amplicon_type].lower()
        self.info_sample   = self.data[0][self.col_sample].lower()

        # cut the last two characters if they are equal to "hp"
        l = len(self.amplicon_type)
        if self.amplicon_type[l - 2:l] == "hp":
            self.amplicon_type = self.amplicon_type[0:l-2]
        return None


    def get_amplicon_type(self):
        """Returns the amplicon type"""

        return self.amplicon_type


    def erase_special_positions(self):
        """adjust input data by deleting columns according to cell type"""
        half_length = len(self.data[0][self.data_col]) / 2

        # import configuration
        import imp
        configuration = imp.load_source('Configuration', 'configuration/configuration.py')

        # change strings from config into lower case
        for i in range(len(configuration.remove_middle_position)):
             configuration.remove_middle_position[i] = configuration.remove_middle_position[i].lower()

        # cell types where defined positions should be removed
        delete_specified_positions = False
        positions = []
        for data in configuration.remove_special_positions:
            if data[0].lower() == self.amplicon_type.lower():
                delete_specified_positions = True
                positions = data[1]
                break
        if delete_specified_positions:
            remaining_positions = len(self.data[0][self.data_col]) - len(positions)
            if not(remaining_positions > 1 and remaining_positions % 2 == 0):
                # column count is even - should be odd before modification
                return False

            positions.sort(reverse=True)
            i = 0
            for full_line in self.data:
                line = full_line[self.data_col]
                for pos in positions:
                    line = line[0:(pos-1)] + line[pos:len(line)]
                self.data[i][self.data_col] = line
                i += 1
            self.mapped_col_count -= len(positions)

        # cell types where the middle position should be removed
        elif self.amplicon_type.lower() in configuration.remove_middle_position:
            if len(self.data[0][self.data_col]) % 2 == 0:
                # column count is even - should be odd before modification
                return False
            i = 0
            for full_line in self.data:
                line = full_line[self.data_col]
                new_line  = line[0:half_length]
                new_line += line[(half_length+1):len(line)]
                self.data[i][self.data_col] = new_line
                i += 1
            self.mapped_col_count -= 1

        # return whether column count is now even
        return len(self.data[0][self.data_col]) % 2 == 0


    def map_lines(self):
        """fold each line in the middle and calculate values for mapped positions"""

        for full_line in self.data:
            line = full_line[self.data_col]

            # map this line
            result = HairpinFunctionality.map_hairpin_positions(line)

            # if not all positions are mutated
            if (not (result.count("0") == 0 and
                     result.count("1") == 0 and
                     result.count("3") == 0 and
                     result.count("4") == 0)):
                # append results to results set
                self.results.append([result, full_line])

        self.mapped_col_count //= 2
        self.results_unsorted = self.results
        self.result_row_count = len(self.results)

        return None


    def delete_mapped_positions(self):
        """Delete specified columns of the mapped data"""

        # import configuration
        import imp
        configuration = imp.load_source('Configuration', 'configuration/configuration.py')

        for data in configuration.delete_mapped_columns:
            amplicon_type = data[0]

            if amplicon_type == self.amplicon_type:

                positions = data[1]
                if len(positions) > 0:
                    positions.sort(reverse=True)

                    # check whether position is in allowed range
                    if positions[0] > len(self.results[0][0]):
                        print("Error in HairpinData.delete_mapped_positions(): Column number to delete is higher than column count for amplicon type " + amplicon_type + ".")
                        return False
                    elif (positions[len(positions)-1] < 1):
                        print("Error in HairpinData.delete_mapped_positions(): Column number to delete is smaller than 1 for amplicon type " + amplicon_type + ".")
                        return False
                    else:
                        i = 0
                        for full_line in self.results:
                            line = full_line[0]
                            for pos in positions:
                                line = line[0:(pos-1)] + line[pos:len(line)]
                                self.results[i][0] = line
                                self.results_unsorted[i][0] = line
                            i += 1
                        self.mapped_col_count -= len(positions)
        return True


    def bind_additional_data_to_results(self):
        """Bind additional information to each line if available
            Linker:   Bind linker CpG pattern and delete all lines
                      which do not have a match either in data or
                      linker file.
            non CpGs: Bind non CpG methylation pattern to each line
            SNPs:     Bind SNP     methylation pattern to each line
        """

        # add id to the beginning of the data and sort after that id
        # (this id is removed at the end of this method)
        for i in range(len(self.results)):
            # add the id of this line to the front of the list. The
            # data now looks as follows: [id, mapped_CpGs, [data_line]]
            self.results[i].insert(0, self.results[i][1][0])
        self.results.sort()

        # import configuration
        import imp
        configuration = imp.load_source('Configuration', 'configuration/configuration.py')
        digits = configuration.digits_after_decimal_point_for_percentages

        # go through all additional data paths (linker, non CpGs, SNPs)
        main_data_sorted = False
        for filepath_additional in self.file_paths_additional:

            # check whether a file with the specified additional information exists
            filename = filepath_additional[1] + self.sub_path
            if filepath_additional[0] == "snp":
                filename += self.filename_datafile_snp
            else:
                filename += self.filename_datafile
            if filename == "" or (not os.path.exists(filename)):
                # if there is no file for this specific additional information, add a placeholder
                for i in range(len(self.results)):
                    self.results[i].append("")

                # remember that this additional data is not available
                self.additional_data_available[filepath_additional[0]] = False
                continue

            # read the file
            try:
                f = open(filename, 'r')
            except IOError as e:
                print("Error in class HairpinData.bind_additional_data_to_results(): The file '" + filename + "' could not be found.")
                return False
            additional_data = []
            first_line = True
            for line in f:
                if line.strip() != "":
                    if not first_line:
                        line_splitted = line.split()
                        if filepath_additional[0] == "snp":
                            # handle SNP-files differently
                            snp_id = line_splitted.pop(0)
                            snp_positions = line_splitted     # rest of the list without the id
                            if len(snp_positions) % 2 != 0:
                                print()
                                print("Error: Count of SNP positions in file " + filename + " is odd.")
                                raise InternalConflictException()
                            # fold the list in the middle and map positions on to each other
                            mapping_correct   = 0
                            mapping_incorrect = 0
                            for i in range(len(snp_positions) // 2):
                                l = snp_positions[i]
                                r = snp_positions[len(snp_positions)-i-1]
                                # ignore positions with mutations
                                if ((l == "a" or l == "c" or l == "g" or l == "t") and
                                    (r == "a" or r == "c" or r == "g" or r == "t")):
                                    # check whether the given mapping is a correct one
                                    if ((l == "a" and r == "t") or
                                        (l == "t" and (r == "a" or r == "g")) or
                                        (l == "c" and r == "g") or
                                        (l == "g" and (r == "c" or r == "t"))):
                                        mapping_correct += 1
                                    else:
                                        mapping_incorrect += 1
                            percentage_str = "_"
                            if mapping_correct + mapping_incorrect > 0:
                                percentage = float(mapping_correct) / (mapping_correct + mapping_incorrect)
                                percentage_str = str(round(percentage, digits))
                            snp_positions.append(percentage_str)
                            additional_data.append([snp_id, snp_positions])
                        else:
                            # line must have as least as many columns as the highest number we want to read
                            add_col_to_cap = self.columns_to_capture_add[filepath_additional[0]]
                            if len(line_splitted) >= add_col_to_cap[len(add_col_to_cap) - 1]:
                                tmp_data = []
                                for col in add_col_to_cap:
                                    # TODO: code might be optimized regarding performance here:
                                    tmp_id = line_splitted[col-1]
                                    # in case of non_cpg
                                    if filepath_additional[0] == "non_cpg":
                                        # convert all "NN" parts of the id (col=1) into "CG"
                                        if col == 1:
                                            tmp_id = tmp_id.replace("NN", "CG")
                                    tmp_data.append(tmp_id)
                                additional_data.append(tmp_data)
                            else:
                                print("Error: The following line (additional data file: " + filepath_additional[0] + ") cannot be read:")
                                print(line)
                                print()
                                raise InternalConflictException()
                    else:
                        first_line = False

            # To be able to speed up matching of IDs in CpG data and additional data,
            # we sort the data by IDs. Doing so we can run through each list only
            # once (with two pointers in each list), finding all matches (instead
            # of scanning the whole second list for each entry in the first list).

            # sort additional data (we don't have to change anything
            # in the list first, since the id is already the first element)
            additional_data.sort()

            # go through both lists and combine additional data with data
            # if ID is equal
            dp = 0   # data pointer
            ap = 0   # additional data pointer
            finished_scan = False

            # check whether ids are equal
            if filepath_additional[0] == "linker":
                self.data_linker_match           = 0
                self.data_linker_under_threshold = 0
                self.data_linker_data_only       = 0
                self.data_linker_linker_only     = 0
            new_results = []
            while not finished_scan and len(self.results) > 0:
                if self.results[dp][0] == additional_data[ap][0]:
                    # two lines have the same id
                    if filepath_additional[0] == "linker":
                        if (1 - float(additional_data[ap][2])) >= 0.8:
                            # mean CG methylation is under or equal to threshold
                            new_results.append([self.results[dp][0], self.results[dp][1], self.results[dp][2], [additional_data[ap][1], additional_data[ap][2]]])
                            self.data_linker_match += 1
                        else:
                            # mean CG methylation is above threshold
                            self.data_linker_under_threshold += 1
                    elif filepath_additional[0] == "non_cpg":
                        new_results.append([self.results[dp][0], self.results[dp][1], self.results[dp][2], self.results[dp][3], additional_data[ap][1]])
                    elif filepath_additional[0] == "snp":
                        new_results.append([self.results[dp][0], self.results[dp][1], self.results[dp][2], self.results[dp][3], self.results[dp][4], additional_data[ap][1]])
                    else:
                        # There should be no option left (if still another left
                        # stop the program).
                        raise InternalConflictException()
                    dp += 1
                    ap += 1
                elif self.results[dp][0] < additional_data[ap][0]:
                    # data line without matching additional data line
                    if filepath_additional[0] == "linker":
                        self.data_linker_data_only += 1
                    elif filepath_additional[0] == "non_cpg":
                        new_results.append([self.results[dp][0], self.results[dp][1], self.results[dp][2], self.results[dp][3], ""])
                    elif filepath_additional[0] == "snp":
                        new_results.append([self.results[dp][0], self.results[dp][1], self.results[dp][2], self.results[dp][3], self.results[dp][4], ""])
                    dp += 1
                elif self.results[dp][0] > additional_data[ap][0]:
                    # additional data line without matching data line
                    if filepath_additional[0] == "linker":
                        self.data_linker_linker_only += 1
                    ap += 1
                if dp == len(self.results) or ap == len(additional_data):
                    if filepath_additional[0] == "linker":
                        self.data_linker_data_only   += len(self.results) - dp
                        self.data_linker_linker_only += len(additional_data) - ap
                    # scan of data is finished
                    finished_scan = True

            self.results = new_results
            if filepath_additional[0] == "linker":
                # save data with linker information
                self.additional_data_available["linker"] = True
                self.result_row_count = self.data_linker_match
            elif filepath_additional[0] == "non_cpg":
                self.additional_data_available["non_cpg"] = True
            elif filepath_additional[0] == "snp":
                self.additional_data_available["snp"] = True

        # delete id from the beginning of the data 
        # (has been added at the beginning of this method)
        for i in range(len(self.results)):
            self.results[i].pop(0)


    def delete_incorrect_snps(self):
        """Delete subset from results where snp percentage is
           lower then 1.0 but not all positions are mutated"""

        if self.additional_data_available["snp"]:
            new_results = []
            self.data_snp_correct   = 0
            self.data_snp_incorrect = 0
            for i in range(len(self.results)):
                percentage = "_"
                if self.results[i][4] != "":
                    percentage = self.results[i][4][len(self.results[i][4])-1]
                if percentage == "_" or float(percentage) == 1.0:
                    new_results.append(self.results[i])
                    self.data_snp_correct   += 1
                else:
                    new_results.append(self.results[i])    #zeile löschen, wenn nicht passende snps aussortiert werden sollen
                    self.data_snp_incorrect += 1
            self.results = new_results
            self.row_count = len(self.results)
        return None


    def calculate_statistics(self):
        """ Calculates all statistics for the completely processed data"""

        # Calculate statistics centralized here

        self.info_00 = 0
        self.info_01 = 0
        self.info_03 = 0
        self.info_04 = 0
        self.info_05 = 0
        self.info_06 = 0
        self.info_08 = 0
        self.info_09 = 0

        self.info_stats_04      = 0
        self.info_stats_043_041 = 0
        self.info_stats_0       = 0
        self.info_stats_01      = 0
        self.info_stats_03      = 0
        self.info_stats_013     = 0
        self.info_stats_0134    = 0

        self.info_stats_only_one = 0
        self.info_stats_mosaic   = 0
        self.info_stats_cont     = 0

        for full_line in self.results:
            result = full_line[0]

            res_count_0 = result.count("0") # save these four variables for reuse
            res_count_1 = result.count("1")
            res_count_3 = result.count("3")
            res_count_4 = result.count("4")
            self.info_00 += res_count_0
            self.info_01 += res_count_1
            self.info_03 += res_count_3
            self.info_04 += res_count_4
            self.info_05 += result.count("5")
            self.info_06 += result.count("6")
            self.info_08 += result.count("8")
            self.info_09 += result.count("9")

            # Are all positions mutated?
            all_positions_mutated = False
            if (res_count_0 == 0 and
                res_count_1 == 0 and
                res_count_3 == 0 and
                res_count_4 == 0):
                all_positions_mutated = True

            # only add this line to the results, if not all positions are mutated
            if not all_positions_mutated:

                # calculate additional statistics
                if result.find("4") != -1:
                    # result line contains 4
                    if (result.find("3") != -1) and (result.find("1") != -1):
                        self.info_stats_0134 += 1
                    elif (result.find("3") != -1) or (result.find("1") != -1):
                        self.info_stats_043_041 += 1
                    else:
                        self.info_stats_04 += 1
                else:
                    if (result.find("3") != -1) and (result.find("1") != -1):
                        self.info_stats_013 += 1
                    elif result.find("3") != -1:
                        self.info_stats_03 += 1
                    elif result.find("1") != -1:
                        self.info_stats_01 += 1
                    else:
                        self.info_stats_0 += 1

                # define the structure of this hairpin
                structure = HairpinFunctionality.define_structure(result)
                if structure == "only_one":
                    self.info_stats_only_one += 1
                elif structure == "mosaic":
                    self.info_stats_mosaic += 1
                elif structure == "continuous":
                    self.info_stats_cont += 1
                elif structure == "no_meth":
                    pass

        return None

    def sort_results(self):
        """Sort the result rows by a specially defined order"""
        # sorting the output

        sorting_array = []
        for a in range(len(self.results)):
            sort_line = ""
            sum = 0
            numbers_count = self.mapped_col_count
            for i in range(self.mapped_col_count):

                char = self.results[a][0][i]
                if char == "1":
                    sort_line += "2"
                    sum += 6
                elif char == "3":
                    sort_line += "3"
                    sum += 6
                elif char == "4":
                    sort_line += "9"
                    sum += 42
                elif char == "0":
                    sort_line += "1"
                    sum += 1
                elif (char == "5" or
                      char == "6" or
                      char == "8" or
                      char == "9"):
                    sort_line += "0"
                    numbers_count -= 1
            mean = 0
            if numbers_count != 0:
                mean = sum / numbers_count

            if mean < 10:
                prefix_mean = "0"
            else:
                prefix_mean = ""

            if sum < 10:
                prefix_sum = "000"
            elif sum < 100:
                prefix_sum = "00"
            elif sum < 1000:
                prefix_sum = "0"
            else:
                print("Error in HairpinData.sort_results(): The sum exceeds a maximum total of 1000")
                assert(False)

            sort_line = prefix_mean + str(mean * 100) + "_" + prefix_sum + str(sum) + "_" + sort_line
            sorting_array.append((sort_line, self.results[a]))

        sorting_array.sort()
        for a in range(len(self.results)):
            self.results[a] = sorting_array[a][1]

        return None


    def print_results(self):
        """print results of mapping to standard output"""
        print()
        print('The following values have been mapped:')
        print()
        for x in self.results:
            print('  ', x[0])
        print()
        return None


    def write_results_to_file(self, output_filename):
        """Writes results into specified text file"""

        # check whether results have been calculated correctly
        if not self.results_available:
            print("Can't write results to file. There was an error in the result calculation.")
            return None

        output_filename = output_filename.strip()
        if output_filename == "":
            "Error in HairpinData.write_results_to_file(): No filename specified."
            return None

        try:
            f = open(output_filename, 'w')
        except IOError as e:
            print("Error in class HairpinData.write_results_to_file(): The file '" + output_filename + "' could not be written.")
            return None

        # write first line with metadata
        for i in (range(len(self.metadata))):
            if i == self.data_col:
                f.write("Mapped CpG Methylation Pattern")
            elif i == self.conversion_col:
                f.write("Conversion Read")
                if self.additional_data_available["linker"]:
                    f.write("\tConversion Linker")
            else:
                f.write(self.metadata[i])
            f.write("\t")

        if self.additional_data_available["linker"]:
            # for each digit of the linker CpG Methylation pattern
            # add one column (l1, l2, ...)
            for i in range(len(self.results[0][2][0])):
                f.write("l" + str(i+1) + "\t")

        # determine length of non-cpg data
        non_cpg_length = 0
        for i in range(len(self.results)):
            if self.results[i][3] != "":
                non_cpg_length = len(self.results[i][3])
                break
        if self.additional_data_available["non_cpg"] and non_cpg_length > 0:
            # for each digit of the non_CPG Methylation pattern
            # add one column (c1, c2, ...)
            for i in range(non_cpg_length):
                f.write("c" + str(i+1) + "\t")

        # determine length of SNP data
        snp_length = 0
        for i in range(len(self.results)):
            if self.results[i][4] != "":
                snp_length = len(self.results[i][4])
                break
        if self.additional_data_available["snp"] and snp_length > 0:
            # for each digit of the SNP Methylation pattern
            # add one column (snp1, snp2, ...)
            for i in range(snp_length):
                f.write("snp" + str(i+1))
                if i < snp_length - 1:
                    f.write("\t")
        f.write("\n")


        # write data
        for x in self.results:
            for i in range(len(x[1])):
                if i == self.data_col:
                    # mapped methylation pattern
                    f.write(x[0])
                elif i == self.conversion_col:
                    # conversion read
                    f.write(str(float(x[1][i])))
                    if self.additional_data_available["linker"]:
                        # conversion linker
                        f.write("\t" + str(1.0 - float(x[2][1])))
                else:
                    f.write(x[1][i])
                f.write("\t")
            # the linker information
            if self.additional_data_available["linker"]:
                for i in range(len(x[2][0])):
                    f.write(x[2][0][i] + "\t")
            # the non CpG information
            if self.additional_data_available["non_cpg"]:
                for i in range(non_cpg_length):
                    if x[3] != "":
                        f.write(x[3][i] + "\t")
                    else:
                        f.write("_\t")
            # the SNP information
            if self.additional_data_available["snp"]:
                for i in range(snp_length):
                    if x[4] != "":
                        f.write(x[4][i])
                    else:
                        f.write("_")
                    if i < snp_length - 1:
                        f.write("\t")
            f.write("\n")
        f.close()

        return None



    def create_image(self, output_filename):
        """Creates a png image file of the results"""

        # check whether results have been calculated correctly
        if not self.results_available:
            print("Can't create image file. There was an error in the result calculation.")
            return None

        # import configuration
        import imp
        configuration = imp.load_source('Configuration', 'configuration/configuration.py')

        # the dimensions of the final image file in pixels
        image_file_dim = configuration.heatmap_image_size

        from PIL import Image
        from PIL import ImageDraw

        def get_column_width(image_width, column_count, column_number):
            left_column_width = image_width // 5
            if column_number == 0:
                return left_column_width
            data_width = image_width - left_column_width
            data_width_plain = data_width - (column_count - 1)
            column_width = data_width_plain // column_count
            if data_width_plain % column_count >= column_number:
                column_width += 1
            return column_width

        # load color definitions from Configuration file
        color_frame = configuration.color_frame
        color_background = configuration.color_background
        color_frame_left_column = configuration.color_frame_left_column
        color_unmethylated      = configuration.color_unmethylated
        color_methylated_left   = configuration.color_methylated_left
        color_methylated_right  = configuration.color_methylated_right
        color_methylated_both   = configuration.color_methylated_both
        color_mutated           = configuration.color_mutated

        # ensure a minimum size for the image
        min_allowed_image_file_dim = (40, 70)
        if image_file_dim[0] < min_allowed_image_file_dim[0]:
            image_file_dim = (min_allowed_image_file_dim[0], image_file_dim[1])
        if image_file_dim[1] < min_allowed_image_file_dim[1]:
            image_file_dim = (image_file_dim[0], min_allowed_image_file_dim[1])

        self.result_row_count = len(self.results)
        data_dim  = (self.mapped_col_count, self.result_row_count)

        box_dim = (20, 10)
        column_size = (box_dim[0], data_dim[1] * box_dim[1])
        img_columns = []
        draw_columns = []
        for i in range(data_dim[0]):
            img_columns.append(Image.new("RGB", column_size, color_frame))
            draw_columns.append(ImageDraw.Draw(img_columns[i]))    # for drawing on column image

        # total count for (0:not methylated, 1:left methylated, 3:right methylated, 4:both methylated)
        total_count_1 = 0
        total_count_2 = 0
        total_count_3 = 0
        total_count_4 = 0
        line_count = 0
        for full_line in self.results:
            line = full_line[0]
            for i in range(self.mapped_col_count):
                # drawing the individual data points
                x1 = 0
                y1 = box_dim[1]*line_count
                x2 = box_dim[0]
                y2 = box_dim[1]*(line_count+1)
                box = (x1, y1, x2, y2)
                if line[i] == "0":
                    draw_columns[i].rectangle(box, fill=color_unmethylated)
                    total_count_1 += 1
                elif line[i] == "1":
                    draw_columns[i].rectangle(box, fill=color_methylated_left)
                    total_count_2 += 1
                elif line[i] == "3":
                    draw_columns[i].rectangle(box, fill=color_methylated_right)
                    total_count_3 += 1
                elif line[i] == "4":
                    draw_columns[i].rectangle(box, fill=color_methylated_both)
                    total_count_4 += 1
                else:
                    draw_columns[i].rectangle(box, fill=color_mutated)
            line_count += 1

        # create the final image file
        img = Image.new("RGB", image_file_dim, color_frame)
        draw = ImageDraw.Draw(img)  # for drawing on final image

        # for each column
        first_column_width = get_column_width(image_file_dim[0], data_dim[0], 0)
        first_column_distance = first_column_width // 4   # distance between first column and data
        col_position = first_column_width
        for i in range(data_dim[0]):
            # scale the column such that it fits into the final image
            col_width = get_column_width(image_file_dim[0], data_dim[0], i+1)
            img_columns[i] = img_columns[i].resize((col_width, image_file_dim[1]), Image.NEAREST)
            # copy the column into the final image
            box = (0, 0, col_width, image_file_dim[1])
            region = img_columns[i].crop(box)
            box = (col_position, 0, col_position + col_width, image_file_dim[1])
            img.paste(region, box)
            col_position += col_width + 1

        # adding the leftmost column
        x1 = 0
        y1 = 1
        x2 = first_column_width - first_column_distance - 1
        y2 = image_file_dim[1] - 1
        total = total_count_1 + total_count_2 + total_count_3 + total_count_4 * 1.0
        if total != 0:
            percentage = (total_count_1 / total, total_count_2 / total, total_count_3 / total, total_count_4 / total)
        else:
            print("Error in HairpinData.create_image(): Total is zero - division by zero is not possible.")
            return False

        # draw four boxes for each percentage one
        limit = 0.00000001
        if percentage[3] > limit:
            box = (x1, y1, x2, y2)
            draw.rectangle(box, fill=color_methylated_both) # fully methylated
        if percentage[2] > limit:
            box = (x1, y1, x2, y2 * (percentage[0] + percentage[1] + percentage[2]))
            draw.rectangle(box, fill=color_methylated_right) # methylated right
        if percentage[1] > limit:
            box = (x1, y1, x2, y2 * (percentage[0] + percentage[1]))
            draw.rectangle(box, fill=color_methylated_left) # methylated left
        if percentage[0] > limit:
            box = (x1, y1, x2, y2 * percentage[0])
            draw.rectangle(box, fill=color_unmethylated) # unmethylated

        # drawing lines separating the percentages
        first_line  = y2 * percentage[0]
        second_line = y2 * (percentage[0] + percentage[1])
        third_line  = y2 * (percentage[0] + percentage[1] + percentage[2])
        if first_line > 3:
            draw.line((0, first_line    , x2, first_line    ), fill = color_frame_left_column)
            draw.line((0, first_line + 1, x2, first_line + 1), fill = color_frame_left_column)
        if (second_line - first_line) > 2:
            draw.line((0, second_line    , x2, second_line    ), fill = color_frame_left_column)
            draw.line((0, second_line + 1, x2, second_line + 1), fill = color_frame_left_column)
        if (third_line - second_line) > 2:
            draw.line((0, third_line    , x2, third_line    ), fill = color_frame_left_column)
            draw.line((0, third_line + 1, x2, third_line + 1), fill = color_frame_left_column)

        # drawing frame around leftmost column
        draw.line((0, 0, 0, image_file_dim[1]), fill = color_frame_left_column) # left border line
        draw.line((1, 0, 1, image_file_dim[1]), fill = color_frame_left_column)
        draw.line((0, 0, x2, 0), fill = color_frame_left_column) # top border line
        draw.line((1, 1, x2, 1), fill = color_frame_left_column)
        draw.line((x2, 0, x2, image_file_dim[1]), fill = color_frame_left_column) # right border line
        draw.line((x2 - 1, 0, x2 - 1, image_file_dim[1]), fill = color_frame_left_column)
        draw.line((0, image_file_dim[1] - 1, x2, image_file_dim[1] - 1), fill = color_frame_left_column) # bottom border line
        draw.line((0, image_file_dim[1] - 2, x2, image_file_dim[1] - 2), fill = color_frame_left_column)

        # draw box between first column and data
        box = (x2 + 1, 0, first_column_width - 1, image_file_dim[1])
        draw.rectangle(box, fill=color_background)

        if configuration.heatmap_column_separator_transparent:
            # make box between first column and data transparent
            mask=Image.new('L', image_file_dim, color=255)
            draw=ImageDraw.Draw(mask)
            draw.rectangle(box, fill=0)
            img.putalpha(mask)

        # saving the final image
        img.save(output_filename.strip(), "PNG")

        return None


    def get_linker_match_information(self):
        """returns information about the matching process (data <--> linker data)"""
        return [
                self.data_linker_match,
                self.data_linker_under_threshold,
                self.data_linker_data_only,
                self.data_linker_linker_only
               ]


    def get_amplicon_meta_information(self):
        """Returns the meta information for the summary file. This list needs
           to have the same order as the list returned by get_amplicon_information()!"""

        # meta information for this amplicon
        meta_information = ["Reference", "Amplicon", "reads", "CpG", "4", "3", "1", "0", "na", "5", "6", "8", "9",
                            "m04", "m043_041", "m01", "m03", "m013", "m0134", "m0",
                            "only_one_methylation", "mosaic_pattern", "continuous_methylation",
                            "Data<->Linker: Conversion rate linker >= 80%)", "Conversion rate linker",
                            "Linker under threshold", "Linker not available", "SNPs perfect match"
                           ]

        # add a column for each CpG position of the linker
        if self.additional_data_available["linker"]:
            # for each digit of the linker CpG Methylation pattern
            # add one column (l1, l2, ...)
            for i in range(len(self.results[0][2][0])):
                meta_information.append("l" + str(i+1))

        # figure out length of non CpG methylation pattern
        non_cpg_length = 0
        for i in range(len(self.results)):
            if self.results[i][3] != "":
                non_cpg_length = len(self.results[i][3])
                break
        # add a column for each CpG position of the non CpG data
        if self.additional_data_available["non_cpg"] and non_cpg_length > 0:
            # for each digit of the non CpG Methylation pattern
            # add one column (l1, l2, ...)
            for i in range(non_cpg_length):
                meta_information.append("c" + str(i+1))

        # create line for result-file
        b = 1
        line = ""
        for info in meta_information:
            line += str(info)
            if b < len(meta_information):
                line += "\t"
            b += 1
        line += "\n"
        return line



    def get_amplicon_information(self):
        """Returns a list with information about this amplicon. This list needs
           to have the same order as the list returned by get_amplicon_meta_information()!"""

        self.info_amplicon = self.amplicon_type
        self.info_reads = self.result_row_count
        self.info_cpg = self.info_00 + self.info_01 + self.info_03 + self.info_04
        self.info_na  = self.info_05 + self.info_06 + self.info_08 + self.info_09

        # import digit information for rounding
        import imp
        configuration = imp.load_source('Configuration', 'configuration/configuration.py')
        digits = configuration.digits_after_decimal_point_for_percentages

        # conversion rate of the linker for all lines that are presented
        # both in the data files and in the linker files with same id
        data_linker_conversion_rate_linker = "no_linker_data"
        if self.additional_data_available["linker"]:
            # import configuration
            data_linker_conversion_rate_linker = round(float(self.data_linker_match) / (self.data_linker_match + self.data_linker_under_threshold), digits)

        # conversion rate of the linker from the summary.dat file from
        # the linker data folder
        conversion_rate_linker = "no_linker_data"
        mean_meth_level = self.get_linker_mean_methylation_level()
        if mean_meth_level:
            conversion_rate_linker = round(1.0 - mean_meth_level, digits)

        # percentage of linker under threshold from matching lines
        linker_under_threshold = "_"
        if self.additional_data_available["linker"]:
            linker_under_threshold = round(float(self.data_linker_under_threshold) / (self.data_linker_match + self.data_linker_under_threshold), digits)

        # percentage of lines where no linker data was available
        linker_not_available = "_"
        if self.additional_data_available["linker"]:
            linker_not_available = round(float(self.data_linker_data_only) / (self.data_linker_match + self.data_linker_under_threshold + self.data_linker_data_only), digits)

        # percentage of lines with SNP information that had a perfect match
        snp_perfect_match = "_"
        if self.additional_data_available["snp"]:
            snp_perfect_match = round(float(self.data_snp_correct) / (self.data_snp_correct + self.data_snp_incorrect), digits)

        amp_inf = [self.info_sample,
                   self.info_amplicon,
                   self.info_reads,
                   self.info_cpg,
                   self.info_04,
                   self.info_03,
                   self.info_01,
                   self.info_00,
                   self.info_na,
                   self.info_05,
                   self.info_06,
                   self.info_08,
                   self.info_09,

                   self.info_stats_04,
                   self.info_stats_043_041,
                   self.info_stats_01,
                   self.info_stats_03,
                   self.info_stats_013,
                   self.info_stats_0134,
                   self.info_stats_0,

                   self.info_stats_only_one,
                   self.info_stats_mosaic,
                   self.info_stats_cont,

                   data_linker_conversion_rate_linker,
                   conversion_rate_linker,
                   linker_under_threshold,
                   linker_not_available,
                   snp_perfect_match
                  ]

        # for each CpG position of the linker we add the percentage of
        # 1s out of the set {0,1}
        if self.additional_data_available["linker"]:
            # for each position of the linker CpG Methylation pattern
            # add up count of zeros and count of ones
            linker_length = len(self.results[0][2][0])
            counts = []
            for i in range(linker_length):
                counts.append([0,0])
            # for each result line
            for x in self.results:
                # for each CpG position 
                for i in range(linker_length):
                    # remember count for zeros and ones
                    if x[2][0][i] == "0":
                        counts[i][0] += 1
                    elif x[2][0][i] == "1":
                        counts[i][1] += 1
            for i in range(linker_length):
                if (counts[i][0] + counts[i][1]) != 0:
                    percentage = float(counts[i][1])/ (counts[i][0] + counts[i][1])
                    amp_inf.append(round(percentage, digits))
                else:
                    percentage = "na"

        # figure out length of non CpG methylation pattern
        non_cpg_length = 0
        for i in range(len(self.results)):
            if self.results[i][3] != "":
                non_cpg_length = len(self.results[i][3])
                break

        # for each non CpG position we add the percentage of
        # 1s out of the set {0,1}
        if self.additional_data_available["non_cpg"] and non_cpg_length > 0:
            # for each position of the non CpG Methylation pattern
            # add up count of zeros and count of ones
            counts = []
            for i in range(non_cpg_length):
                counts.append([0,0])
            # for each result line
            for x in self.results:
                if x[3] != "":
                    # for each non CpG position 
                    for i in range(non_cpg_length):
                        # remember count for zeros and ones
                        if x[3][i] == "0":
                            counts[i][0] += 1
                        elif x[3][i] == "1":
                            counts[i][1] += 1
            for i in range(non_cpg_length):
                if (counts[i][0] + counts[i][1]) != 0:
                    percentage = float(counts[i][1])/ (counts[i][0] + counts[i][1])
                    amp_inf.append(round(percentage, digits))
                else:
                    percentage = "na"



        # create line for result-file
        b = 1
        line = ""
        for info in amp_inf:
            line += str(info)
            if b < len(amp_inf):
                line += "\t"
            b += 1
        line += "\n"
        return line


    def get_linker_mean_methylation_level(self):
        """Reads the mean methylation level from the linkers summary.dat file"""

        # take first additional filepath which must be the linker filepath
        filename = self.file_paths_additional[0][1] + self.sub_path + self.filename_summary_file

        # read the specified file
        try:
            f = open(filename, 'r')
        except IOError as e:
            print("Error in function get_mean_methylation_level(): The file '" + filename + "' could not be found.")
            return False

        cur_line = 0
        # extract all available summary information
        information = {}
        for line in f:
            if line.strip() != "":
                ls = line.split("\t")
                if len(ls) == 2:
                    information[ls[0]] = ls[1]

        if "Mean methylation level" in information:
            return float(information["Mean methylation level"])
        else:
            return 2.0


    def get_additional_data_available(self, data_name):
        """ Returns whether there has been linker information found for this amplicon"""
        return self.additional_data_available[data_name]


class HairpinFunctionality:
    """The functionality needed for analysis of hairpin data"""

    @staticmethod
    def map_hairpin_positions(data):
        """Folds the provided string in the middle and returns a half-sized string with the mapped information"""

        # ensure that the provided string is at least of length 2
        if len(data) < 2: raise InvalidInputException()

        # ensure that the provided string is of even length
        if len(data) % 2 == 1: raise InvalidInputException()

        # ensure that only valid characters have been provided
        if len(data) != (data.count("0") + data.count("1") + data.count("x")): raise InvalidInputException()

        # calculate mapped string
        half_length = len(data) // 2
        left_side = data[0:half_length]
        right_side = data[half_length:len(data)]
        right_side = right_side[::-1]   # reverse order of string

        result = ""
        for i in range(half_length):
            result_char = 's'
            if left_side[i] == '0' and right_side[i] == '0':
                # both sides are unmethylated
                result_char = '0'
            elif left_side[i] == '1' and right_side[i] == '0':
                # left side is methylated, right side unmethylated
                result_char = '1'
            elif left_side[i] == '0' and right_side[i] == '1':
                # left side is unmethylated, right side is methylated
                result_char = '3'
            elif left_side[i] == '1' and right_side[i] == '1':
                # left side is unmethylated, right side is methylated
                result_char = '4'
            elif ((left_side[i] == '0' and right_side[i] == 'x') or
                  (left_side[i] == 'x' and right_side[i] == '0')):
                # one side is unmethylated, the other mutated
                result_char = '5'
            elif left_side[i] == '1' and right_side[i] == 'x':
                # left side is methylated, right side is mutated
                result_char = '6'
            elif left_side[i] == 'x' and right_side[i] == '1':
                # left side is methylated, right side is mutated
                result_char = '8'
            elif left_side[i] == 'x' and right_side[i] == 'x':
                # left side is methylated, right side is mutated
                result_char = '9'
            assert result_char != 's'   # if the result char is still 's', something went wrong
            result += result_char
        return result

    @staticmethod
    def define_structure(data):
        """Analyse whether strand has 'only_one', 'mosaic' or 'continuous' pattern"""
        # We want to relate the given data to one of the following groups:
        # no_meth    : No methylation at either sides.
        # only_one   : One side has exactly one methylation, the other side
        #              has zero or one methylations.
        #              Examples: 040, 100, 300, 130
        # mosaic     : At least one of the two sides shows a mosaic methylation
        #              pattern.
        #              Examples: 101, 303, 104, 134
        # continuous : The data is not mosaic and at least one side has a
        #              methylation directly followed by another one.
        #              Examples: 110, 133, 440, 410
        pass

        count_4 = data.count("4")
        count_3 = data.count("3")
        count_1 = data.count("1")

        # no methylation on either sides
        if count_4 == 0 and count_3 == 0 and count_1 == 0:
            return "no_meth"

        # only one methylation
        if ((count_4 == 1 and count_3 == 0 and count_1 == 0) or
            (count_4 == 0 and ((count_3 == 1 and count_1 == 1) or
                               (count_3 == 1 and count_1 == 0) or
                               (count_3 == 0 and count_1 == 1)))):
            return "only_one"

        # mosaic methylation pattern
        state_left  = "b"
        state_right = "b"
        # go through characters one after another and transit through
        # finit automaton. The following states are distinguished:
        # b, 0, 1, 10 and 101, b, 0, 3, 30 and 303
        for character in data:
            if character == "0" or character == "1" or character == "3" or character == "4":
                if character == "0":
                    # left side
                    if state_left == "b":
                        state_left = "0"
                    elif state_left == "1":
                        state_left = "10"
                    # right side
                    if state_right == "b":
                        state_right = "0"
                    elif state_right == "3":
                        state_right = "30"
                elif character == "1":
                    # left side
                    if state_left == "b" or state_left == "0":
                        state_left = "1"
                    elif state_left == "10":
                        return "mosaic"
                    # right side
                    if state_right == "b":
                        state_right == "0"
                    elif state_right == "3":
                        state_right = "30"
                elif character == "3":
                    # left side
                    if state_left == "b":
                        state_left = "0"
                    elif state_left == "1":
                        state_left = "10"
                    # right side
                    if state_right == "b" or state_right == "0":
                        state_right = "3"
                    elif state_right == "30":
                        return "mosaic"
                elif character == "4":
                    # left side
                    if state_left == "b" or state_left == "0":
                        state_left = "1"
                    elif state_left == "10":
                        return "mosaic"
                    # right side
                    if state_right == "b" or state_right == "0":
                        state_right = "3"
                    elif state_right == "30":
                        return "mosaic"

        # continuous methylation pattern
        state_left = "b"
        state_right = "b"
        # go through characters one after another and transit through
        # finit automaton. The following states are distinguished:
        # b, 0, 1 and 11, b 0, 3, 33
        for character in data:
            if character == "0" or character == "1" or character == "3" or character == "4":
                if character == "0":
                    if state_left == "b":
                        state_left = "0"
                    elif state_left == "1":
                        state_left = "0"
                    if state_right == "b":
                        state_right = "0"
                    elif state_right == "3":
                        state_right = "0"
                if character == "1":
                    if (state_left == "b") or (state_left == "0"):
                        state_left = "1"
                    elif state_left == "1":
                        return "continuous"
                    if (state_right == "b") or (state_right == "1"):
                        state_right = "0"
                if character == "3":
                    if (state_left == "b") or (state_left == "1"):
                        state_left = "0"
                    if (state_right == "b") or (state_right == "0"):
                        state_right = "3"
                    elif state_right == "3":
                        return "continuous"
                if character == "4":
                    if (state_left == "b") or (state_left == "0"):
                        state_left = "1"
                    elif state_left == "1":
                        return "continuous"
                    if (state_right == "b") or (state_right == "0"):
                        state_right = "3"
                    elif state_right == "3":
                        return "continuous"

        # now all cases are checked, code reaching therefore indicates
        # an error in the above code
        raise InternalConflictException()


class InvalidInputException(Exception):
    """Exception: The provided input is not in a valid form"""
    pass


class InternalConflictException(Exception):
    """Exception: There is an internal conflict, which cannot be resolved"""
    pass
