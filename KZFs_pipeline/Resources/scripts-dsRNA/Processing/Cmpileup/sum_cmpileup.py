#python

# region Imports
import argparse
import csv
import inspect
import logging
import os
import sys
from datetime import datetime
from pprint import pformat

# use this for importing external self commands
if __name__ == '__main__':
    sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from common.init_logging_dict import init_logging_dict
from common.files_general import *

# endregion

# region Logging
LOG_DIR = os.path.join("%(out_dir)s", "PoolCmpileup.%(time)s." +
                       "%(file_suffix)s.log")
RUN_PARAMS_DEBUG_MSG = "Running with Resource Files and Run Paths: %s"
FILES_FOUND_DEBUG_MSG = "Paths Found: %s"
FILES_NOT_FOUND_DEBUG_MSG = "No files ending with '%(suffix)s', starting with one of %(samples)s and containing '%(fragment)s' were found, quitting group"
READ_FILE_INFO_MSG = "Reading %s"
WRITING_OUTPUT_MSG = "Writing output into %s"
FILE_EXISTS_NO_FORCE_MSG = "File %s exists, force overwrite or delete"
DONE_MSG = "Done combining RNA Editing Index cmpileups, %s generated"
# endregion

DELIMITER = "\t"

# region Cmpileup definitions
CMPILEUP_REGION_HEADER_FIELDS = ["RegionChr", "RegionStart", "RegionEnd"]
CMPILEUP_POSITION_HEADER_FIELDS = ["PositionInRegion"]
CMPILEUP_REFERENCE_HEADER_FIELDS = ["Reference"]
CMPILEUP_BASE_COUNT_HEADER_FIELDS = ["TotalCoverage", "A", "C", "G", "T", "UnrecognizedBases", "LowQualityBases"]
# List keeps order, get dictionary of {column name : column index}
CMPILEUP_HEADER = {h: i for i, h in enumerate(
    CMPILEUP_REGION_HEADER_FIELDS + CMPILEUP_POSITION_HEADER_FIELDS + CMPILEUP_REFERENCE_HEADER_FIELDS + CMPILEUP_BASE_COUNT_HEADER_FIELDS)}

OUTPUT_HEADER = ["PositionChr", "PositionStart", "PositionEnd",
                 "Region"] + CMPILEUP_REFERENCE_HEADER_FIELDS + CMPILEUP_BASE_COUNT_HEADER_FIELDS


# endregion

def create_position_key(line, group_cols):
    '''
    Returns the key of each region, that is - the unchanging part of the cmpileup (all but the base columns)
    :param group_cols: list or set of group columns
    :param line: list of cmpileup fields
    :return: position str
    '''
    return DELIMITER.join(line[CMPILEUP_HEADER[h]] for h in group_cols)


def get_formatted_region_coordinates(region):
    ''' returns region coordinates in format chr:start-end'''
    return {
        "Region": region[CMPILEUP_HEADER["RegionChr"]] + ":" + region[CMPILEUP_HEADER["RegionStart"]] + "-" + region[
            CMPILEUP_HEADER["RegionEnd"]]}


def get_position_coordinates(region):
    ''' returns position location '''
    return {"PositionChr": region[CMPILEUP_HEADER["RegionChr"]],
            "PositionStart": int(region[CMPILEUP_HEADER["PositionInRegion"]]) - 1,
            "PositionEnd": region[CMPILEUP_HEADER["PositionInRegion"]]}

def get_formatted_reference(region):
    ''' returns formatted reference '''
    return {CMPILEUP_REFERENCE_HEADER_FIELDS[0] : region[CMPILEUP_HEADER[CMPILEUP_REFERENCE_HEADER_FIELDS[0]]].upper()}

def write_output(he_regions_dict, output_file):
    '''
    TODO fix output region - split back into columns
    :param group_col:
    :param he_regions_dict: regions dictionary to write
    :param output_file: path to output file
    :return: None
    '''
    logging.info(WRITING_OUTPUT_MSG % output_file)
    with open(output_file, 'w') as csvfile:
        # header treats group columns as one field, as that is the key : sums implementation
        # this also ensures our wanted order (to create a BED file)
        writer = csv.DictWriter(csvfile, lineterminator="\n", fieldnames=OUTPUT_HEADER)
        # write header
        writer.writeheader()
        for region, values in he_regions_dict.items():
            # get region information
            region = region.split(DELIMITER)
            values.update(get_formatted_reference(region))
            values.update(get_position_coordinates(region))
            values.update(get_formatted_region_coordinates(region))

            # write to file
            writer.writerow(values)


def add_position_to_sum(cols_to_sum, group_cols, he_regions_dict, line):
    '''
    Adds base count of line to file
    :param cols_to_sum:
    :param group_cols:
    :param he_regions_dict:
    :param line:
    :return:
    '''
    # get grouping variable (ES site)
    key = create_position_key(line, group_cols)
    # if this is the first time we see the key - create dictionary for key
    if not key in he_regions_dict:
        # insert initial sums
        he_regions_dict[key] = {s_col: int(line[CMPILEUP_HEADER[s_col]]) for s_col in cols_to_sum}
    else:
        # add sum to dictionary
        for s_col in cols_to_sum:
            he_regions_dict[key][s_col] += int(line[CMPILEUP_HEADER[s_col]])

def main(root_dir, file_suffix, output_file, log_path, force, sample_list=None, fragment=None, 
         group_cols = CMPILEUP_REGION_HEADER_FIELDS + CMPILEUP_POSITION_HEADER_FIELDS + CMPILEUP_REFERENCE_HEADER_FIELDS, 
         cols_to_sum = CMPILEUP_BASE_COUNT_HEADER_FIELDS):
    # init logging
    init_logging_dict(log_path if log_path else LOG_DIR % {'out_dir': root_dir,
                                                            'file_suffix': file_suffix,
                                                            'time': datetime.today().isoformat()},)
    args, _, _, values = inspect.getargvalues(inspect.currentframe())
    logging.debug(RUN_PARAMS_DEBUG_MSG % pformat(dict([(i, values[i]) for i in args])))

    if os.path.exists(output_file) and not force:
        logging.warning(FILE_EXISTS_NO_FORCE_MSG % output_file)
        return

    # find all files
    files = find_files(root_dir, file_suffix)
    # filter by string if needed
    if fragment:
        files = [f for f in files if fragment in f]
    # filter by sample list if needed
    if sample_list:
        files = [f for f in files if any([os.path.basename(f).startswith(s) for s in sample_list])]
    logging.info(FILES_FOUND_DEBUG_MSG % files)
    
    # if no files were found - issue warning and return
    if not files:
        logging.error(FILES_NOT_FOUND_DEBUG_MSG % {'suffix' : file_suffix, 'samples' : sample_list, 'fragment' : fragment})
        # return bad exit code
        raise IOError(FILES_NOT_FOUND_DEBUG_MSG % {'suffix' : file_suffix, 'samples' : sample_list, 'fragment' : fragment})

    # initialize
    he_regions_dict = {}

    # for each CSV file
    for fp in files:
        # info
        logging.info(READ_FILE_INFO_MSG % fp)
        # read file
        with open(fp, "r") as f:
            # read line
            line = f.readline().strip().split(DELIMITER)
            # while line exists
            while line[0]:
                # add line information
                add_position_to_sum(cols_to_sum, group_cols, he_regions_dict, line)
                # read next line (or EOF)
                line = f.readline().strip().split(DELIMITER)

    # write to file
    write_output(he_regions_dict, output_file)

    # notify of finish
    logging.info(DONE_MSG % output_file)
    
    return 0
    

if __name__ == '__main__':
    # Create parser
    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass


    parser = argparse.ArgumentParser(formatter_class=MyFormatter,
                                     description='Find repeats that are in opposite orientation at given regions, by families. Run on CentOS 7 for use of bedtools.2.27.1')
    parser.add_argument('-d', '--root_dir', dest='root_dir', action='store', metavar='root_dir',
                        nargs='?', required=True, help='The input directory')
    parser.add_argument('-s', '--file_suffix', dest='file_suffix', action='store', metavar='file_suffix',
                        nargs='?',default='_ucscHg38Alu.bed.gz_mpileup.cmpileup', help='Wanted file suffix')
    parser.add_argument('-o', '--output_file', dest='output_file', action='store', metavar='output_file', nargs='?',
                        required=True, help='Output file, default is to create or append in input dir')
    parser.add_argument('-l', '--log_path', dest='log_path', action='store', metavar='log_path', nargs='?',
                        default=None, help='Log file, default is to create in input dir')
    parser.add_argument('--force', dest='force_out', action='store_true', default=False,
                        help='Force overwrite of CSV files even if they exist')
    parser.add_argument('--fragment', dest='file_fragment', action='store', default=None,
                            help='Only include file with given string as a fragment of filename, fixed')
    parser.add_argument('--sample_list', dest='sample_list', action='store', nargs='+', default=None,
                            help='Only include samples from list of samples, fixed. Samples will be detected by concatenating sample name with file suffix')


    options = parser.parse_args()

    main(root_dir=options.root_dir,
         file_suffix=options.file_suffix,
         output_file=options.output_file,
         log_path=options.log_path,
         force=options.force_out,
         fragment=options.file_fragment,
         sample_list=options.sample_list
         )
