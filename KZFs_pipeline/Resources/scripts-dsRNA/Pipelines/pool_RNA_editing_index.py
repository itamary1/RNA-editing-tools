#! /bin/python

# region Imports
import argparse
import inspect
import logging
import os
import sys
from multiprocessing import Pool
from datetime import datetime
from pprint import pformat

# use this for importing external self commands
if __name__ == '__main__':
    sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from common.init_logging_dict import init_logging_dict
from common.files_general import *
from common.command_executer import execute
from Processing.Cmpileup import sum_cmpileup
from consts.Regions import *
from consts.Sites import *
# endregion

# region Logging
LOG_DIR = os.path.join("%(out_dir)s", "ProcessRNAEditingIndex%(suffix)s.%(time)s.log")
RUN_PARAMS_DEBUG_MSG = "Running with Resource Files and Run Paths: %s"
FILES_FOUND_DEBUG_MSG = "Paths Found: %s"
READ_FILE_INFO_MSG = "Reading %s"
WRITING_OUTPUT_MSG = "Writing output into %s"
FILE_EXISTS_NO_FORCE_MSG = "File %s exists, force overwrite or delete"
FORCE_RUN_MSG = "Force run: removing files: %s"
DONE_MSG = "Done combining RNA Editing Index cmpileups"
# endregion

# region Programs 
PYTHON_CMD="python3.6"
R_CMD="Rscript"
SHELL_CMD="sh"
# endregion


def init_logging(log_path):
    # init logging
    init_logging_dict(log_path)
    args, _, _, values = inspect.getargvalues(inspect.currentframe())
    logging.debug(RUN_PARAMS_DEBUG_MSG % pformat(dict([(i, values[i]) for i in args])))


def pool_regions(options):
    """Runs the main pipeline for pooling the StrandDerivingCountsPerRegion.csv files per region per sample.
    Args:
        input_dir (str): input directory
        output_dir (str): output directory
        log_path (str): log file path
        group_file (str): path of group classification of CSV file (Sample,Group )
    """
    #input_dir, output_dir, log_path, group_file

    # initialize log and write params
    init_logging(options.log_path if options.log_path else LOG_DIR % {'out_dir': options.input_dir,'suffix': "Regions",
                                                                       'time': datetime.today().isoformat()})

    # change file names
    logging.info("Changing top-level StrandDerivingCountsPerRegion.csv file names")
    execute(CHANGE_FILE_NAMES_CMD % {'input_dir': options.input_dir})

    # check if files exist
    files_with_suffix = find_files_by_suffix(input_dir=options.input_dir, suffix=SPLIT_OUT_SUFFIX)
    # if files exist, and force was activated, delete
    if options.force and files_with_suffix:
        logging.info(FORCE_RUN_MSG % (str(files_with_suffix)))
        remove_files(files_with_suffix)
        files_with_suffix = None

    # splitting files if needed (found + not force)
    if files_with_suffix:
        logging.info("Skipping splitting information by chromosome, found files with suffix %s" % SPLIT_OUT_SUFFIX)
        logging.debug("Files found: " + str(files_with_suffix))
    else: 
        logging.info("Splitting information by chromosomes")
        execute(SPLIT_BY_CHR_CMD % {'r': R_CMD, 'split_by_chr': SPLIT_BY_CHR_SCRIPT, 'input_dir' : options.input_dir})

    # pooling information
    logging.info("Pooling editing information by chromosomes")
    execute(POOL_BY_CHR_CMD % {'r': R_CMD, 'pool_by_chr': POOL_BY_CHR_SCRIPT, 
                                'input_dir' : options.input_dir, 'output_dir' : options.output_dir,
                                'group_file': options.group_file})

    # combine files
    logging.info("Combining editing info files to dataset-level")
    # pooling
    execute(COMBINE_FILES_CMD % {'output_dir' : options.output_dir})
    # counts (fix alt\fix chr - so a different treatment than previous files)
    execute(COMBINE_COUNTS_CMD % {'output_dir' : options.output_dir})

    logging.info("Done processing region data")


def get_groups(group_file):
    """Get set of groups in data to iterate over

    Args:
        group_file (str): path to CSV of format Sample,Group

    Returns:
        set: set of unique groups
    """
    with open(group_file, "r") as f:
        groups = {l.split(",")[1].strip() for l in f.readlines()[1:]}
    return groups


def get_samples(group, group_file):
    """Get list of samples classified as given group

    Args:
        group (str): group str as defined by the get_groups method
        group_file (str): path to CSV of format Sample,Group

    Returns:
        str: str of concatenated list of samples
    """
    with open(group_file, "r") as f:
        samples = {l.split(",")[0].strip() for l in f.readlines()[1:] if l.split(",")[1].strip() == group}
    return samples

def pool_sites_for_group(options, group, output_dir_cmpileups):
    # "flatten" cmpileup across groups by summing each group
    # call imported method    
    logging.info("Pooling cmpileup by group for %s" % group)
    try:
        sum_cmpileup.main(root_dir=options.input_dir,
                        file_suffix=options.cmpileup_suffix,
                        output_file=POOLED_OUT_FORMAT % {'output_dir' : output_dir_cmpileups, 'cmpileup' : options.cmpileup_suffix, 'group' : group}, 
                        log_path=options.log_path,
                        force=options.force,
                        sample_list=get_samples(group, options.group_file))
    except IOError:
        logging.error("Could not file files in %(input_dir)s ending with '%(suffix)s' and starting with one of %(samples)s" % {'suffix' : options.cmpileup_suffix, 'samples' : get_samples(group, options.group_file), 'input_dir' : options.input_dir})
        # exit function for this group
        return
    
    # convinience variable (cleaner code)
    out_format_dict = {'output_dir' : output_dir_cmpileups, 'cmpileup' : options.cmpileup_suffix, 'group' : group}
    logging.info("Finished pooling cmpileup by group for %s" % group)
    
    # split into SNP and non-SNP points for each file
    # if force run needed - delete files
    if options.force and (os.path.exists(POOLED_OUT_NO_SNP_FORMAT % out_format_dict) or 
                            os.path.exists(POOLED_OUT_W_SNP_FORMAT % out_format_dict)):
        logging.info(FORCE_RUN_MSG % (str([POOLED_OUT_NO_SNP_FORMAT % out_format_dict, POOLED_OUT_W_SNP_FORMAT % out_format_dict])))
        remove_files([POOLED_OUT_NO_SNP_FORMAT % out_format_dict, POOLED_OUT_W_SNP_FORMAT % out_format_dict])
    # if not force run but files exist - warn and continue
    if not options.force and (os.path.exists(POOLED_OUT_NO_SNP_FORMAT % out_format_dict) or os.path.exists(POOLED_OUT_W_SNP_FORMAT % out_format_dict)):
        logging.warning(FILE_EXISTS_NO_FORCE_MSG % str([POOLED_OUT_NO_SNP_FORMAT % out_format_dict, POOLED_OUT_W_SNP_FORMAT % out_format_dict]))
    else:
        # run
        logging.info("Splitting into SNP and non-SNP sites for %s" % group)    
        execute(SPLIT_CMPILEUP_BY_SNP_CMD % {'shell' : SHELL_CMD, 
                                             'split_summed_cmpileup_by_SNP' : SPLIT_CMPILEUP_BY_SNP_SCRIPT, 
                                             'output_dir' : output_dir_cmpileups, 
                                             'cmpileup' : options.cmpileup_suffix, 
                                             'group' : group,
                                             'snps' : options.snps_file})
        logging.info("Finished splitting into SNP and non-SNP sites for %s" % group) 
    
    # Pooling editing site information for site and region
    logging.info("Pooling editing site information for %s" % group)
    execute(POOL_ES_CMD % {'r' : R_CMD, 'pool_info_script' : POOL_ES_SCRIPT, 
                           'input_no_snps' : POOLED_OUT_NO_SNP_FORMAT % out_format_dict, 
                           'input_only_snps' : POOLED_OUT_W_SNP_FORMAT % out_format_dict, 
                           'output_dir' : options.output_dir, 
                           'group_name': group, 
                           'group_size' : len(get_samples(group, options.group_file))} + 
            # add force if needed
            (POOL_ES_FORCE_PARAM if options.force else "") + 
            # add new region id if needed
            (POOL_ES_ALTERNATIVE_REGION_ID % options.region_id_file if options.region_id_file else ""))
    logging.info("Finished pooling editing site information for %s" % group)
    
    logging.info("Done processing site data for %s" % group)
    


def pool_sites(options):
    """Runs the pipeline for pooling the cmpileup files of the RNA editing index.

    Args:
        input_dir (str): input directory
        output_dir (str): output directory
        log_path (str): log file path
        group_file (str): path of group classification of CSV file (Sample,Group )
        cmpileup_suffix (str): suffix of cmpileup files to sum
    """
    #input_dir, output_dir, log_path, group_file, cmpileup_suffix
    # initialize log and write params
    if not os.path.exists(options.output_dir):
        os.makedirs(options.output_dir)
    init_logging(options.log_path if options.log_path else LOG_DIR % {'out_dir': options.output_dir,'suffix': "Sites",
                                                                       'time': datetime.today().isoformat()})
    # get groups
    groups = get_groups(options.group_file)
    if groups:
        logging.info("Found groups: " + str(groups))
    else:
        # if no groups found - issue error and quit
        logging.error("No groups were found to work on, check group file (--group_file)")
        exit(1)

    # create directory before parallelization to avoid race conditions
    output_dir_cmpileups = os.path.join(options.output_dir, "CombinedCmpileups")
    if not os.path.exists(output_dir_cmpileups):
        logging.info("Creating output directory for CMpileup files at %s" % output_dir_cmpileups)
        os.makedirs(output_dir_cmpileups)

    # create commands
    pool_sites_run_list = []
    for group in groups:        
        pool_sites_run_list.append((options, group, output_dir_cmpileups))
    
    # run all in parallel
    with Pool(processes=options.processes) as pool:
        pool.starmap(pool_sites_for_group, pool_sites_run_list)
                                
    logging.info("Done processing site data")
    

if __name__ == '__main__':
    # Create parser
    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass


    parser = argparse.ArgumentParser(formatter_class=MyFormatter,
                                     description='Pool and recalculate the RNA editing index according to given parameters')
    # add a subparser
    subparsers = parser.add_subparsers(description='Pooling types')
    # add parser for index pooling
    region_parser = subparsers.add_parser("regions", formatter_class=MyFormatter, help="Pool region information")
    region_parser.add_argument('-i', '--input_dir', dest="input_dir", help='Path of RNA editing index directory (where per-sample results files reside)', required=True)
    region_parser.add_argument('-o', '--output_dir',dest="output_dir", help='Path of pooled results directory', required=True)
    region_parser.add_argument('-l', '--log_path', dest='log_path', action='store', metavar='log_path', nargs='?',
                                default=None, help='Log file, default is to create in output dir')
    region_parser.add_argument('--group_file',dest="group_file", nargs='?', help='CSV with sample to group. Will pool samples according to given groups, sample count per group will be calculated according to this file. Format: Sample,Group',required=True)
    region_parser.add_argument('--force',action="store_true", default=False, dest="force", help='Should force rerun splitting by chromosome?')
    region_parser.set_defaults(func=pool_regions)
    # add parser for site pooling
    site_parser = subparsers.add_parser("sites", formatter_class=MyFormatter, help="Pool site information")
    site_parser.add_argument('-i', '--input_dir', dest="input_dir", help='Path of RNA editing index directory (where cmpileup files reside)', required=True)
    site_parser.add_argument('-o', '--output_dir',dest="output_dir", help='Path of pooled results directory', required=True)
    site_parser.add_argument('-l', '--log_path', dest='log_path', action='store', metavar='log_path', nargs='?',
                                default=None, help='Log file, default is to create in input dir')
    site_parser.add_argument('--group_file',nargs='?', dest="group_file", help='CSV with sample to group. Will pool samples according to given groups, sample count per group will be calculated according to this file. Format: "Sample,Group", where "Sample" is the unique prefix of each sample and Group is a group name of your choice',required=True)
    site_parser.add_argument('-s', '--suff_cmpileup',nargs='?', dest="cmpileup_suffix", default="ucscHg38Alu.bed.gz_mpileup.cmpileup", help='cmpileup files suffix')
    site_parser.add_argument('--snps',nargs='?', dest="snps_file", help='Path of SNP BED file to remove from analysis', default="ReplaceSNPPATH/HomoSapiens/ucscHg38CommonGenomicSNPs150.bed.gz")
    site_parser.add_argument('-rid', '--region_id_file',nargs='?', dest="region_id_file", default=None, help='A CSV file to replace region coordinates as grouping variable during pooling (only affects last step). Must conform to the format "Region,NewRegionID", to ensure proper pooling. Can include a partial list of regions to pool (reduce computational load). This allows calculating pooled editing over sites from non-consecutive genomic regions, by renaming them to a single region ID. Note that an inner join is being made, so that *a genomic region can be joined in different combinations* and that *only transformed region IDs will be kept*.')
    site_parser.add_argument('-p', '--processes',action="store", type=int, default=3, dest="processes", help='Number of samples to run in parallel')
    site_parser.add_argument('--force',action="store_true", default=False, dest="force", help='Force rerun for pooling of cmpileups')
    site_parser.set_defaults(func=pool_sites)

    options = parser.parse_args()
    # call wanted program function
    options.func(options)