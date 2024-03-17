# python 3.9 

import argparse
from datetime import datetime
import gzip
import inspect
import logging
from multiprocessing import Pool
import os
from pprint import pformat
import re
import sys

if __name__ == '__main__':
    sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from common.files_general import *
from common.command_executer_py3 import execute
from common.gzip_operations import is_gzipped
from common.init_logging_dict import init_logging_dict
from common.BedRecord import BedFile

CHR = 0
START = 1
END = 2
STRAND = 5
STRAND_PLUS = "+"
STRAND_MINUS = "-"

def sort_bed_file(bed_file):
    return sorted(bed_file, key=lambda x: (x[0], x[1], x[2]))

def unique_bed_file(bed_file, uniq_on_indexes=None):
    if uniq_on_indexes:
        if isinstance(uniq_on_indexes, slice):
            bed_file = [record[uniq_on_indexes] for record in bed_file]
        else:
            bed_file = [[record[x] for x in uniq_on_indexes] for record in bed_file]
        
    return list(set(tuple(i for i in l) for l in bed_file))

def separate_bed_records(bed_records, uniq_on_indexes):
    logging.info("Sorting uniquely %s records on %s indexes" % (len(bed_records), uniq_on_indexes))
    # make unique and sort first (chr, start, end)
    sorted_unique_records = sort_bed_file(unique_bed_file(bed_records,
                                                          uniq_on_indexes=uniq_on_indexes))
    
    logging.info("Disjoining %s records" % len(sorted_unique_records))
    # initialize
    disjoint_sets = []
    # iterate over records
    for record in sorted_unique_records:
        # iterate over sets
        for disjoint_set in disjoint_sets:
            # get last records
            last_record = disjoint_set[-1]
            # if the last added records was before current (end of last is smaller than start of current) 
            # no book-ended features
            # OR on another chromosome
            if last_record[CHR] != record[CHR] or last_record[END] < record[START]:
                # append this records
                disjoint_set.append(record)
                break
        else:
            # did not add to any set
            # create a new one
            disjoint_sets.append([record])
            
    logging.info("Found %s disjoint sets" % len(disjoint_sets))
    logging.info("Disjoint sets size: %s" % [len(s) for s in disjoint_sets])
    
    return disjoint_sets

def write_beds(res, out_format):
    for i, disjoint_set in enumerate(res, start=1):
        write_dict_to_bed(dict_list=[{CHR : record[CHR], 
                                          START : record[START], 
                                          END : record[END]} for record in sort_bed_file(disjoint_set)], 
                              output_file=out_format % i, 
                              ordered_keys_list=[CHR, START, END], 
                              append=False)

def main(options, file_delim="\t"):
    # # init logging
    init_logging_dict(os.path.join(os.path.dirname(options.output_format), "DisjoinBedRecords%(suffix)s.%(time)s.log" % {'suffix': os.path.basename(options.output_format) % "Log", 'time': datetime.today().isoformat()}))
    args, _, _, values = inspect.getargvalues(inspect.currentframe())
    logging.debug("Running with Resource Files and Run Paths: %s" % pformat(dict([(i, values[i]) for i in args])))
    
    # logging.info("Loading BED files")
    # # iterate over files
    bed_records = []
    for bed_file in options.bed_input_files:
        open_fn = gzip.open if is_gzipped(bed_file) else open
        with open_fn(bed_file, "rt") as handle:
            # add new records
            bed_records.extend([l.strip().split("\t") for l in handle.readlines()])
    logging.debug("Read total of %s records" % len(bed_records))
    
    # convert start and end into int type
    for record in bed_records:
        for index in [START, END]:
            record[index] = int(record[index])
    
    if options.strand:      
        logging.info("Working on plus strand")
        # split into strands
        bed_records_plus = [r for r in bed_records if r[STRAND] == STRAND_PLUS]
        # disjoin 
        res = separate_bed_records(bed_records_plus, uniq_on_indexes=[CHR, START, END, STRAND])
        # write each result
        logging.debug("Writing BED files for plus sets")
        write_beds(res, out_format=options.output_format % "plus.%s")
        logging.debug("Finished writing BED files for plus sets")
        
        logging.info("Working on minus strand")
        # split into strands
        bed_records_minus = [r for r in bed_records if r[STRAND] == STRAND_MINUS]
        # disjoin
        res = separate_bed_records(bed_records_minus, uniq_on_indexes=[0,1,2,5])
        # write each result
        logging.debug("Writing BED files for plus sets")
        write_beds(res, out_format=options.output_format % "minus.%s")
        logging.debug("Finished writing BED files for plus sets")
    
    else:
        # disjoin
        res = separate_bed_records(bed_records, uniq_on_indexes=[CHR, START, END])
        # write each result
        logging.debug("Writing BED files of sets")
        write_beds(res, out_format=options.output_format)
        logging.debug("Finished writing BED files of sets")
    
    logging.info("Done creating disjoint sets")


if __name__ == "__main__":
    # Create parser
    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass


    parser = argparse.ArgumentParser(formatter_class=MyFormatter,
                                     description='Split BED into disjoint sets of non-overlapping, unique BED records')
    parser.add_argument('-i', '--input_files', dest='bed_input_files', action='store',
                        metavar='bed_input_files', nargs='+', required=True,
                        help='Path to BED files to process (can be gzipped)')
    parser.add_argument('-o', '--output_format', dest='output_format', action='store', metavar='output_format', 
                        nargs='?', required=True, help='A format for path of output directory of BED files. Must include a \%s for different sets.')
    parser.add_argument('--strand', dest='strand', action='store_true', default=False, 
                        help='Consider strand when splitting into sets, disjoined sets will be calculated for strands separately and will be output into different files.')

    # parser.add_argument('-p', '--processes', dest='processes', action='store', 
    #                     metavar='processes', nargs='?',
    #                     default=20, type=int, help='Maximal number of processes to run in parallel')
    # parser.add_argument('-m', '--min_len_bed', dest='min_len_bed', action='store', 
    #                     metavar='min_len_bed', nargs='?', default=0, type=int, 
    #                     help='Minimal segment (paired bases) length in filtered BED')

    options = parser.parse_args()
    
      
    # create dictionary of NM_id to gene
    main(options)