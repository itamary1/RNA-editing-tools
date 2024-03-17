import json
import os
import csv

def read_json(file):
    '''
    Reads JSON and returns dictionary
    :param file:
    :return: dict
    '''
    with open(file, "r") as json_file:
        return json.load(json_file)

def find_files(root_dir, suffix, follow_links = False):
    '''
    Find all files with given suffix in a directory recursively
    :param root_dir:
    :param suffix:
    :param follow_links:
    :return: lst
    '''
    found = []
    # go over directories, get just the Log.final.out files, modify for with and without following links
    for root, dirs, files in os.walk(root_dir, followlinks = follow_links):
        for file in files:
            if file.endswith(suffix):
                found.append(os.path.join(root, file))  # add full file path to list of files
    return found

def write_row_from_dict(row_dict, ordered_header_list, file_path, override = False):
    '''
    Write a single row from dictionary into (existing or non-existing) file.
    :param row_dict: dictionary with header : values
    :param ordered_header_list: list of dictionary keys according to wanted order
    :param file_path: CSV path
    :param override: should the file be overwritten and not appended to?
    :return: None
    '''
    # if override file or no header
    if override or not os.path.exists(file_path):
        with open(file_path, 'w') as csvfile:
            # create object
            writer = csv.DictWriter(csvfile, fieldnames=ordered_header_list)
            # write header and row
            writer.writeheader()
            writer.writerow(row_dict)
    else:
        # append to file
        with open(file_path, 'a') as csvfile:
            # create object
            writer = csv.DictWriter(csvfile, fieldnames=ordered_header_list)
            writer.writerow(row_dict)

def write_dict_to_csv(dict_list, ordered_keys_list, output_file, append = False):
    '''
    Writes list of dictionaries to CSV file
    :param dict_to_write: dictionary
    :param ordered_keys_list: wanted order of keys
    :param output_file: output path
    :param append: should the file be overwritten and not appended to?
    :return: None
    '''
    # if file exists and should be appended to
    if append and os.path.isfile(output_file):
        with open(output_file, 'a') as csvfile:
            # header treats group columns as one field, as that is the key : sums implementation
            # this also ensures our wanted order (to create a BED file)
            writer = csv.DictWriter(csvfile, lineterminator="\n", fieldnames=ordered_keys_list)
            writer.writerows(dict_list)
    else:
        # else write file
        with open(output_file, 'w') as csvfile:
            # header treats group columns as one field, as that is the key : sums implementation
            # this also ensures our wanted order (to create a BED file)
            writer = csv.DictWriter(csvfile, lineterminator="\n", fieldnames=ordered_keys_list)
            # write header
            writer.writeheader()
            writer.writerows(dict_list)

def remove_files(file_list):
    for f in file_list:
        try:
            os.remove(f)
        except OSError:
            pass

def listdir_fullpath(d):
    return [os.path.join(d, f) for f in os.listdir(d)]

def find_files_by_suffix(input_dir, suffix):
    return [f for f in listdir_fullpath(input_dir) if os.path.isfile(f) and f.endswith(suffix)]

def write_dict_to_tsv(dict_list, ordered_keys_list, output_file, append = False):
    '''
    Writes list of dictionaries to TSV file
    :param dict_to_write: dictionary
    :param ordered_keys_list: wanted order of keys
    :param output_file: output path
    :param append: should the file be overwritten and not appended to?
    :return: None
    '''
    # if file exists and should be appended to
    if append and os.path.isfile(output_file):
        with open(output_file, 'a') as csvfile:
            # header treats group columns as one field, as that is the key : sums implementation
            # this also ensures our wanted order (to create a BED file)
            writer = csv.DictWriter(csvfile, lineterminator="\n", fieldnames=ordered_keys_list,
                                    delimiter='\t', quoting=csv.QUOTE_NONE)
            writer.writerows(dict_list)
    else:
        # else write file
        with open(output_file, 'w') as csvfile:
            # header treats group columns as one field, as that is the key : sums implementation
            # this also ensures our wanted order (to create a BED file)
            writer = csv.DictWriter(csvfile, lineterminator="\n", fieldnames=ordered_keys_list,
                                    delimiter='\t', quoting=csv.QUOTE_NONE)
            # write header
            writer.writeheader()
            writer.writerows(dict_list)
            
def write_dict_to_bed(dict_list, ordered_keys_list, output_file, append = False):
    '''
    Writes list of dictionaries as BED file (tab-delimited, no header)
    :param dict_to_write: dictionary
    :param ordered_keys_list: wanted order of keys
    :param output_file: output path
    :param append: should the file be overwritten and not appended to?
    :return: None
    '''
    # if file exists and should be appended to
    if append and os.path.isfile(output_file):
        with open(output_file, 'a') as csvfile:
            # header treats group columns as one field, as that is the key : sums implementation
            # this also ensures our wanted order (to create a BED file)
            writer = csv.DictWriter(csvfile, lineterminator="\n", fieldnames=ordered_keys_list,
                                    delimiter='\t', quoting=csv.QUOTE_NONE)
            writer.writerows(dict_list)
    else:
        # else write file
        with open(output_file, 'w') as csvfile:
            # header treats group columns as one field, as that is the key : sums implementation
            # this also ensures our wanted order (to create a BED file)
            writer = csv.DictWriter(csvfile, lineterminator="\n", fieldnames=ordered_keys_list,
                                    delimiter='\t', quoting=csv.QUOTE_NONE)
            # do not write header
            writer.writerows(dict_list)