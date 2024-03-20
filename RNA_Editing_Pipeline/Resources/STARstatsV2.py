# DESCRIPTION: summerizes the STAR statistics for each sample
# ver 2 07.08.2018: match no _ at filename end (nameLog.final.out)

import os.path
import argparse
import sys

# must use this for importing external self commands
if __name__ == '__main__':
    sys.path.append(os.path.dirname(__file__))
from common.command_executer import execute

# # Create parser
# class SortingHelpFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
#     def add_arguments(self, actions):
#         actions = sorted(actions, key=attrgetter('option_strings'))
#         super(SortingHelpFormatter, self).add_arguments(actions)

# Create parser
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass

# region Constants
# define separator for csv
CSV_DELIMITER = ','
LOG_FILE_STR = 'STAR_StatsLog.txt'
STATS_FILE_STR = 'STAR_stats.csv'
STAR_LOG_FILE_STR = 'Log.final.out'
# endregion

# region Arguments Parsing
parser = argparse.ArgumentParser(formatter_class=MyFormatter,
                                 description='Generate a tab-delimited CSV file with alignment statistics (drawn from Log.final.out files). Rerunnning for plot generation does not overwrite CSV file by default.')
parser.add_argument('-r', '--root', dest='STAR_dir', action='store', metavar='STAR_dir', nargs='?', required=True, help='Path to STAR output directory')
parser.add_argument('-o', '--output', dest='output_dir', action='store', metavar='output_dir', nargs='?', required=True, help='Path to output directory')
parser.add_argument('-l', '--log', dest='log_path', action='store', metavar='log_path', nargs='?', default=None, help='Directory where log file should be created')
parser.add_argument('--force', dest='force_out', action='store_true', default=False, help='Force overwrite of CSV files even if they exist')
parser.add_argument('--follow_links', dest='follow_links', action='store_true', default=False, help='Flag: Should follow links in directory?')
parser.add_argument('--plot', dest='plot_stats', action='store_true', default=False, help='Flag: Should script plot alignment statistics?')
parser.add_argument('-po', dest='plot_output_dir', action='store', metavar='plot_output_dir', default='.',  nargs='?', help='Path to output directory for plots')
parser.add_argument('-g', dest='group_file', action='store', metavar='group_file', default=None, help='Optional when plotting statistics: divide graphs to facets by given groups. File should contain two columns: sample then group, tab delimited')
parser.add_argument('--prefix', dest='output_prefix', action='store', metavar='output_prefix', default=None, help='Add prefix to output files names')
parser.add_argument('-a', dest='all_plots', action='store_true', default=False, help='Optional when plotting statistics: create extra statistics graphs (not just read mapping distribution)')
arguments = parser.parse_args()
# endregion

# region Paths And Files
# get output file
if arguments.output_prefix:
    outfile = os.path.join(arguments.output_dir, arguments.output_prefix + "_" + STATS_FILE_STR)
else:
    outfile = os.path.join(arguments.output_dir, STATS_FILE_STR)

# open log file
if arguments.log_path:
    log_file_name = os.path.join(arguments.log_path, LOG_FILE_STR)
else:
    log_file_name = os.path.join(arguments.output_dir, LOG_FILE_STR)

# create output file directory if does not exist
if not os.path.exists(arguments.output_dir):
    try:
        os.makedirs(arguments.output_dir)
    except OSError as exc: # Guard against race condition
        if exc.errno != errno.EEXIST:
            raise
log_file = open(log_file_name, mode = 'w')
# endregion

# region Create Stats File
#if this file does not exist or is forced overwrite - calculte all again
if (not os.path.exists(outfile)) or arguments.force_out:
    logs = []
    # go over directories, get just the Log.final.out files, modify for with and without following links
    if arguments.follow_links:
        for root, dirs, files in os.walk(arguments.STAR_dir, followlinks = True):
            for file in files:
                if file.endswith(STAR_LOG_FILE_STR):
                    logs.append(os.path.join(root, file))  # add full file path to list of files
    else:
        for root, dirs, files in os.walk(arguments.STAR_dir):
            for file in files:
                if file.endswith(STAR_LOG_FILE_STR):
                    logs.append(os.path.join(root, file))  # add full file path to list of files

    log_file.write("STAR log files found: " + str(logs) + "\n")

    # create output list
    final = []
    first = True

    for stats in logs:
        # get lines from file
        stats = open(stats, 'r')
        lines = stats.readlines()
        stats.close()

        # filter just lines with "|"
        lines = [l.strip() for l in lines if "|" in l]
        # remove first 4 lines - not needed stats (runtime etc.)
        lines = lines[4:]

        # add header
        if first:
            header = ["Sample"]
            # split header (splittable lines only)
            header.extend([(l.split("|")[0]).strip() for l in lines])
            #replace "," by nothing
            header = [str.replace(r, ',', '') for r in header]
            final.append(CSV_DELIMITER.join(header))
            first = False

        # split all lines by | and take just second part
        sample = str.replace(os.path.basename(stats.name), STAR_LOG_FILE_STR, '')  # get file name only
        row = [sample]  # add sample
        row.extend([(l.split("|")[1]).strip() for l in lines])  # get statistics
        final.append(CSV_DELIMITER.join(row))

    # remove % from numbers
    final[1:] = [str.replace(r, '%', '') for r in final[1:]]

    log_file.write("Output file: " + outfile + "\n")

    # write output
    with open(outfile, "w") as f:
        f.write("\n".join(final) + "\n")

    log_file.write("DONE: statistics file generated\n")
else:
    log_file.write("Statistics file exists\n")
# endregion


# region Plot Stats
if arguments.plot_stats:
    command = ["Rscript", "/home/alu/twerski/Scripts/Nextflow/Levanon_lab_NEXTFLOW_PIPELINE/Resources/STAR_stats.R", "-i", outfile, "-o", arguments.plot_output_dir, "-l", log_file_name]
    if arguments.group_file:
        command.extend(["-g", arguments.group_file])
    if arguments.output_prefix:
        command.extend(["-p", arguments.output_prefix])
    if arguments.all_plots:
        command.append("-a")

    command = " ".join(command)
    log_file.write("Plotting alignment statistics: " + command + "\n")

    execute(command, "Cannot create plots")

    log_file.write("DONE: plots created\n")
# endregion

log_file.close()
