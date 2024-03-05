import subprocess
import multiprocessing
import concurrent.futures
import random
import os
import re




"""when Near_site_folding_calculator run in parallel it will create temporery files with the pid as the uniq key
to prevent clashing between sequences so: !!!!! you cant run more then one parallellic run from the same procces !!!!"""

"""
A class for runing running Fold from RNAstructure package
"""
class Fold_runner:
    def __init__(self,dataPath,which_Fold):
        """constructor
        Args:
            which_Fold (str): the command to run Fold in your enviroment
            dataPath (str): a path to data_tables directory of download with the RNAstructure package
        """
        self.which_Fold = which_Fold
        self.dataPath = dataPath

    def run_Fold(self, seq, Fname="-", convertToRNA=True, temprature=303.15):
        """will run fold on give seq, will return its output if no out-path given
        otherwise will write out to the out path

        Args:
            seq (str): rna sequnce
            Fname (str, optional): path to out file. Defaults to "-" == stdout(captured).
            convertToRNA (bool, optional): will replace t with u. Defaults to True.
            temprature (float, optional): temprature of folding. Defaults to 303.15.
        Returns:
            str: stdout - will be empty if u choose writing to file
        """
        if convertToRNA:
            seq = seq.replace("T", "U")
        p = subprocess.run(
            [f"export DATAPATH={self.dataPath}; echo {seq} | {self.which_Fold} '-' {Fname} -mfe -q --temperature {temprature} -k"], capture_output=True, text=True, shell=True
        )
        return p.stdout

    """good only for bracet format"""

    def parse_Fold_file(self, file):
        """will parse fold out file

        Args:
            file (str): path to fold out file

        Returns:
            tuple: (energy, structure) - (0,0) if the seqeunce were unfolded
        """
        with open(file, "r") as ff:
            data = ff.readlines()
        data = [line.strip("\n") for line in data]
        return self.parse_fold_out("\n".join(data))

    def parse_fold_out(self, data):
        """will get fold output as a string, and will return tuple of (energy, structure) - (0,0) if the seqeunce were unfolded"""
        data = data.split("\n")
        energy_line = data[0].strip("\n")
        if "ENERGY" in energy_line:
            energy = -float(energy_line.split("=")[1].split("-")[1])
            structure = data[2].strip("\n")
            return energy, structure
        else:
            return (0, 0)

    """run Fold and return output"""



class Near_site_folding_calculator:
    """this class porpose is to calculate the delta-G energy of the dsRNA (only) where the editing site reside
    - for  a given seq with editing point in the middle
    """
    def __init__(self,dataPath,which_Fold,which_bpRNA,temp_dir="/tmp/"):
        """constructor

        Args:
            which_Fold (str): the command to run Fold in your enviroment
            dataPath (str): a path to data_tables directory of download with the RNAstructure package
            which_bpRNA (str): a command to run bpRNA - you may need have perl with Graph.pm.
            temp_dir (str, optional): _description_. Defaults to "/tmp/".
        """
        self.temp_dir = os.path.join(temp_dir, "Fold_tmp_fixed/")
        # create dir for temp files
        if not os.path.exists(self.temp_dir):
            os.umask(0)
            os.mkdir(self.temp_dir, mode=0o777)
        self.bp = which_bpRNA
        self.fold_runner = Fold_runner(dataPath,which_Fold)

    def run_bpRNA(self, file):
        os.chdir(self.temp_dir)
        p = subprocess.run([self.bp, file], capture_output=True, text=True)
        assert not p.stdout, "bpRNA cant run file: " + file
        # create out file name
        out_f = file.split("/")[-1]
        out_f = out_f.removesuffix(os.path.splitext(out_f)[1]) + ".st"
        # return absolut path to out f
        return self.temp_dir + out_f

    def parseStfile(self, file, wind):
        with open(file, "r") as bpf:
            data = bpf.readlines()
        cre = re.compile(r"(\s\d+\.\.\d+\s)")
        for line in data:
            if "segment" in line:
                l = cre.split(line)
                part1 = l[1].strip()
                part1 = {
                    "start": int(part1.split(".")[0]),
                    "end": int(part1.split(".")[-1]),
                }
                part2 = l[3].strip()
                part2 = {
                    "start": int(part2.split(".")[0]),
                    "end": int(part2.split(".")[-1]),
                }
                # if the editing site is in that segmant
                if (part1["start"] <= wind + 1 <= part1["end"]) or (
                    part2["start"] <= wind + 1 <= part2["end"]
                ):
                    coords = (
                        part1["start"],
                        part1["end"],
                        part2["start"],
                        part2["end"],
                    )
                    seqs = (l[2], l[-1].strip("\n"))
                    return (coords, seqs)
        # if we our editng site is in any segment
        return (0, 0)

    def remove_files(self, files):
        [os.remove(file) for file in files]

    def calculate_substrcture_one_seq(self, seq, wind, out_fName=None, temprature=303.15):
        if out_fName:
            out_fName = self.temp_dir + out_fName + ".dbn"
        else:
            out_fName = self.temp_dir + str(random.randint(0, 100000)) + ".dbn"
        # check out file with that name dosent already exist
        assert not os.path.exists(out_fName), "fold out file already exist file = " + out_fName
        self.fold_runner.run_Fold(seq, out_fName, temprature=temprature)
        # check ford run succesfuily
        assert os.path.exists(out_fName), "fold didnt create outfile, \n fname = {} \n seq = {}\n".format(
            out_fName, seq
        )
        # get structure and enenrgy from fold out
        energy1, struct1 = self.fold_runner.parse_Fold_file(out_fName)
        # if its unfoldable return zeroes
        if not energy1:
            self.remove_files(out_fName)
            return {
                "big_f_e": 0,
                "big_s": 0,
                "small_seq": 0,
                "small_seq_coords": 0,
                "small_f_e": 0,
                "small_s": 0,
            }

        # run bpRNA
        bpRNA_file = self.run_bpRNA(out_fName)
        # parse bprna to get dsRNAseq+coords
        coords, near_seqs = self.parseStfile(bpRNA_file, wind)
        # remove bp file in background
        self.remove_files([bpRNA_file, out_fName])
        # if our site is not in any segment
        if not near_seqs:
            return {
                "big_f_e": energy1,
                "big_s": struct1,
                "small_seq": 0,
                "small_seq_coords": 0,
                "small_f_e": 0,
                "small_s": 0,
            }
        # run fold to variable(dsRNAseq+7*"n")
        new_seq = near_seqs[0] + "N" * 7 + near_seqs[1]
        result = self.fold_runner.run_Fold(new_seq, temprature=temprature)
        # parse energy2+struct2
        energy2, struct2 = self.fold_runner.parse_fold_out(result)
        # return energy1, structure1, dsRNAseq, energy2, struct2 as dict
        return {
            "big_f_e": energy1,
            "big_s": struct1,
            "small_seq": new_seq,
            "small_seq_coords": coords,
            "small_f_e": energy2,
            "small_s": struct2,
        }

    """
    will calculate substrcture for multiple sequences in a list, runing for each seq in parallel

    Args:
        file (str): path to fold out file

    Returns:
        tuple: (energy, structure) - (0,0) if the seqeunce were unfolded
    """
    def calculate_substrcture_seq_list(self, seqs, temprature=303.15):
        assert len(seqs[0]) % 2 == 1, "seq len should be odd not even"
        wind = int((len(seqs[0]) - 1) / 2)
        mpid = str(os.getpid())
        inp_l = [
            {"seq": seq, "wind": wind, "out_fName": mpid + str(i), "temprature": temprature}
            for i, seq in enumerate(seqs)
        ]
        with concurrent.futures.ProcessPoolExecutor() as executor:
            results = executor.map(self.calculate_substrcture_one_seq_oneIn, inp_l)
        return list(results)

    '''
    just an helper function for calculate_substrcture_seq_list for being able to run parallelic to give executor.map one argumnt.
    there is a support of multy arguments now so I dont need it anymore and should update my code but its works so I keep it that way 
    '''
    def calculate_substrcture_one_seq_oneIn(self, inp):
        return (
            self.calculate_substrcture_one_seq(inp["seq"], inp["wind"], inp["out_fName"], inp["temprature"])
            if "out_fName" in inp
            else self.calculate_substrcture_one_seq(inp["seq"], inp["wind"], temprature=inp["temprature"])
        )


