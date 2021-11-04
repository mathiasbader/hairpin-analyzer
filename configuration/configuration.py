# ---------------------------------------------------
# --                                               --
# --    Hairpin Analyzer v0.3                      --
# --                                               --
# --    --> Configuration file                     --
# --                                               --
# --                                               --
# --  2011-Aug-01 to 2011-Aug-23                   --
# --  by Mathias Bader (mail@mathiasbader.de)      --
# --  at Saarland University                       --
# ---------------------------------------------------
#
#
# --------------------------------------------------------------
# -- INFORMATION: Use this file to configure the software.    --
# -- Lines starting with a number sign (#) are ignored by the --
# -- software, their only purpose is to help the user to      --
# -- understand the meaning of the following lines. If you    --
# -- find any bugs or something remains unclear, please       --
# -- contact the author (contact information is given above). --
# --------------------------------------------------------------


# CpG data folder
# ---------------
# Specify were the CpG data you want to analyse is located.
# This is the main data to analyze. The other paths all
# specify data that will be related to this data:
data_path_cpg = "CpGs/"

# Additional data paths
# ---------------------
# Specify were additional data for analysis is located (linker, non_cpg and SNP):
data_paths_additional = [
    ["linker", "Conversion/"],
    ["non_cpg", "nonCpGs/"],
    ["snp", "SNPs/"]
]
# specify a top folder which will be used for all folders
# specified above
data_path_main_folder = "/media/data2TB/Ablage/userData/pascal/MiSeq_Data/BQ-Results/14.02.19_IAP-OxBis/"


# File names
# ----------
# name of the file that should be scanned containing the
# results (data, linker and non_CpG):
filename_results = "results.tsv"
# file for SNPs:
filename_results_snp = "SNPs.tsv"
# name of the file that should be scanned containing the
# summary for one result file:
filename_results_summary = "summary.dat"
# name of the file that should be created containing the
# summary for each amplicon type. The amplicon type will
# be added between filename and file extension
results_summary_filename       = "achieved_results_"
results_summary_file_extension = ".txt"

# Heatmap
# ------------------
# Define the image size of the heatmaps that are produced.
# The minimum size is 40x70 and will be enforced if you
# enter a smaller size.
heatmap_image_size = (220, 367)
# Define whether the are between the leftmost column in the
# heatmap and the other columns should be transparent or not
heatmap_column_separator_transparent = True
# All colors used for the heatmap are defined here
color_frame = "#dddddd"             # color between CpG site columns
color_background = "white"          # color between leftmost column and other columns (if not transparent)
color_frame_left_column = "#444444" # color of the frame around and the lines in the leftmost column
color_unmethylated = "#00ccff"  # color for unmethylated state on both sides
# if primer are designed 3'from restriction side
color_methylated_left  = "#76cc2f"  # color for methylation only on left side,
color_methylated_right = "#B9E6B8"  # color for methylation only on right side, 
# if primer are designed 5' from restriction side
#color_methylated_left  = "#B9E6B8"  # color for methylation only on left side
#color_methylated_right = "#76cc2f"  # color for methylation only on right side

color_methylated_both  = "#ff6600"  # color for methylation on both sides
color_mutated = "white"             # color for mutated state on either side or both
# for changing the colors, I recommend the following tool: http://www.colorpicker.com


# Deleting CpG positions before mapping
# -------------------------------------
# Define for each amplicon, which CpG positions should be deleted
# before any mapping happens. Amplicons having the same behavior
# are grouped together. Amplicons which are not specified here, are
# not changed.
#
# All amplicons with odd position count, where the
# center position should be removed:
remove_middle_position = [
    'l1',
    'mmetnhp',
    'b1',
    'igf2',
    'afp',
    'Oasl1HP',
    'tdg8',
    'tdg9',
    'tdg39',
    'tdg49',
    'tdg54'
]
# In the following list, you can define for an amplicon positions that should be removed
remove_special_positions = [
    ('ysat', [1]),
    ('msat', [1]),
    ('oct4', [9]),
    ('msatub', [9,10]),
    ('iapneu', [1,2]),
    ('tdg36',[9,16])
]


# Deleting columns from heatmap
# -----------------------------
# Define for each amplicon, which mapped columns should be
# removed before creating the heatmap
delete_mapped_columns = [
    ('iap', [1,2,3,4]),
    ('b1',  [1]),
    ('l1', [6,7]),
    ('iapez', [1,2,3,4]),
    ('snrpn', [1,2,3]) # wenn sequenzen zu kurz
]


# Rounding results which are percentages
# --------------------------------------
# When results are calculated in percentages, it might
# happen, that they have infinite number of digits after
# the decimal points. This variable defines how many
# digits should be output to the result file (the value
# will be rounded to the last digit, e.g. 0.435 and
# x=2 gives 0.44). Percentages are always given as numbers
# between 0 and 1.
digits_after_decimal_point_for_percentages = 4
