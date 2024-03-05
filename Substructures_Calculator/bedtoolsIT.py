__author__ = "Itamar T"

"""this module created to enable fast non-temp files impelmention for running bedtools"""

from Bio import SeqIO
import subprocess


class genome_reader:
    """this calss porpose is to enable one time loading of the genome, and use it multiple times to get fasta seqs"""

    def __init__(self, file_p):
        """consturctor - will read fasta file

        Args:
            file_p (str): path to fasta file
        """
        self.data_dict = SeqIO.to_dict(SeqIO.parse(file_p, "fasta"))

    def get_fasta(self, chr, start, end=None):
        """will return fasta seq from start to end if end given,
            otherwise,  will return the the postiaonal nucleotide if only one argument is given for position

        Args:
            chr (str): chr
            start (int): seq wanted start postion if range is given, else - position of wanted nuc
            end (int, optional): end postion if range is given. Defaults to None.

        Returns:
            str: seq
        """
        # if end and round and len(self.data_dict[chr].seq) < end:
        #     pass
        # else:
        return (
            str(self.data_dict[chr].seq)[start : end + 1]
            if end
            else str(self.data_dict[chr].seq)[start : start + 1]
        )

    def get_chromosooms_names(self):
        """will return names of chromosoms in our genome file"""
        return self.data_dict.keys()


def intersect_aStdin(bed_var, bed_file, bedtoolsOptions):
    bedtools_command = "bedtools intersect " + bedtoolsOptions + ' -a "stdin"' + " -b " + bed_file
    command = "echo -e '" + bed_var + "' | " + bedtools_command
    btOut = subprocess.run(command, capture_output=True, text=True, shell=True)
    return btOut.stdout if not btOut.stderr else ""


def intersect_bStdin(bed_var, bed_file, bedtoolsOptions):
    bedtools_command = "bedtools intersect " + bedtoolsOptions + " -a " + bed_file + ' -b "stdin"'
    command = "echo -e '" + bed_var + "' | " + bedtools_command
    btOut = subprocess.run(command, capture_output=True, text=True, shell=True)
    return btOut.stdout if not btOut.stderr else ("ERROR " + btOut.stderr)


if __name__ == "__main__":
    print("\n\n>>>>>>>>>>>>>>>>>>>>>running example<<<<<<<<<<<<<<<<<<<")
    print("\n\ncode :")
    print('gr=genome_reader("/private4/gabayo2/Fasta/S-cerevisiae/sacCer3.fa")')
    print("gr.get_fasta('chrI',100,120)")
    print("\noutput : ")
    gr = genome_reader("/private4/gabayo2/Fasta/S-cerevisiae/sacCer3.fa")
    print(gr.get_fasta("chrI", 100, 120))
    print(
        "NOTE! every instance of genome_reader will have the genome one the RAM, so you better make as few instace as possible"
    )

    print("\n\ncode :")
    print(r"bed_var=Scaffold100019\t157348\t157549\nScaffold100019\t158094\t158095")
    print('intersect_aStdin(bed_var, "/private8/Projects/zohar/RNAstructure/squid/editing_site.bed","-wa")')
    print("\noutput : ")
    bed_var = "Scaffold100019\t157348\t157549\nScaffold800019\t158094\t158095"
    print(intersect_aStdin(bed_var, "/private8/Projects/zohar/RNAstructure/squid/editing_site.bed", "-wa"))
