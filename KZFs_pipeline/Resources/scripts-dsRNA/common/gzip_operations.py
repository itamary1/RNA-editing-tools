# Function for wrapping subprocesses opening
import gzip
import shutil
import os

def is_gzipped(path):
    with open(path, "rb") as f:
        return f.read(2) == b'\x1f\x8b'

def gzip_and_remove_original(file_path):
    # zip file
    with open(file_path, 'rb') as f_in, gzip.open(file_path + ".gz", 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    #remove original
    os.remove(file_path)

def gunzip_file(file_path, suffix='gz'):
    # unzip zip file: remove suffix and '.' from end
    with gzip.open(file_path, 'rb') as f_in, open(file_path[:-(len(suffix)+1)], 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    return file_path[:-3]

def write_gzip_file(file_path, output):
    with gzip.open(file_path, "w") as f:
        f.write(output)
