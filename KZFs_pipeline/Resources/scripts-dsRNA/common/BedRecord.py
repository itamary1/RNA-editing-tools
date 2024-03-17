#!/usr/bin/env python3.8
"""
This module attempts to wrap BED files and work with this data structure
"""
__author__ = 'Roni'


from enum import Enum
class Sign(str, Enum):
    PLUS = "+"
    MINUS = "-"

    @classmethod
    def sign(cls, sign_str):
        if sign_str == Sign.MINUS:
            return Sign.MINUS
        elif sign_str == Sign.PLUS:
            return Sign.PLUS
        raise ValueError("Unrecognized strand sign: %s" % sign_str)

    # def __str__(self):
    #     return self.value
class UnstrandedBedRecord:
    def __init__(self, record_lst):
        self._chr = record_lst[0]
        self._start = int(record_lst[1])
        self._end = int(record_lst[2])
        self._interval = self._chr + ":" + str(self._start) + "-" + str(self._end)

    def __repr__(self):
        return self.interval

    def __len__(self):
        return self.end - self.start

    def __eq__(self, other):
        return self.chr == other.chr and self.start == other.start and self.end == other.start

    def __ne__(self, other):
        return self.chr != other.chr or self.start != other.start or self.end != other.start
    
    def __lt__(self, other):
        # less than
        # bed is smaler
        return self.chr == other.chr and self.end < other.start
    
    def __le__(self, other):
        # less than or equal to
        return self.chr == other.chr and self.end <= other.start
    
    def __gt__(self, other):
        # greater than
        return self.chr == other.chr and other.end < self.start
    
    def __ge__(self, other):
        # greater than or equal to
        return self.chr == other.chr and other.end <= self.start
    
    def __hash__(self):
        return hash((self.chr, self.start, self.end))

    # python-ic constructor overloading using class methods
    @classmethod
    def bed_record_from_str(cls, bed_str, bed_sep="\t"):
        return UnstrandedBedRecord(bed_str.split(bed_sep))

    # region Properties
    @property
    def chr(self):
        return self._chr
    @property
    def start(self):
        return self._start
    @property
    def end(self):
        return self._end
    # endregion

    def bed(self):
        return "\t".join([self.chr, str(self.start), str(self.end)])

    @property
    def interval(self):
        return self._interval

class StrandedBedRecord(UnstrandedBedRecord):
    def __init__(self, record_lst):
        self._name = record_lst[3]
        self._score = record_lst[4]
        self._strand = Sign.sign(record_lst[5])
        
        # invoking the __init__ of the parent class
        super().__init__(record_lst)

    def __repr__(self):
        return self.interval + "(" + self.strand + ")"

    def __eq__(self, other):
        return self.chr == other.chr and self.start == other.start and self.end == other.start and \
               self.name == other.name and self.score == other.score and self.strand == other.strand
    
    def __ne__(self, other):
        return self.chr != other.chr or self.start != other.start or self.end != other.start or \
               self.name != other.name and self.score != other.score and self.strand != other.strand
    
    def __hash__(self):
        return hash((self.chr, self.start, self.end, self.name, self.score, self.strand))

    def __len__(self):
        return self.end-self.start
        
    # python-ic constructor overloading using class methods
    @classmethod
    def bed_record_from_str(cls, bed_str, bed_sep="\t"):
        return StrandedBedRecord(bed_str.split(bed_sep))

    # region Properties
    @property
    def name(self):
        return self._name
    @property
    def score(self):
        return self._score
    @property
    def strand(self):
        return self._strand
    # endregion

    def bed(self):
        return "\t".join([self.chr, str(self.start), str(self.end), self.name, self.score, self.strand])

    def sign(self):
        return self.strand
    
    def nullify_name_and_score(self):
        # make into nothing to make life easier
        return StrandedBedRecord([self.chr, self.start, self.end, "0", "0", self.strand])
        
    def unstrand(self):
        return UnstrandedBedRecord([self.chr, self.start, self.end])
    

class BedFile:
    """Represents a BED file with relevant methods
    """
    def __init__(self, bed_record_list=None, stranded=True):
        self._bed_list = bed_record_list if bed_record_list else []
        self._record_class = StrandedBedRecord if stranded else UnstrandedBedRecord
        self._stranded = stranded
        

    # python-ic constructor overloading using class methods
    @classmethod
    def bed_file_from_str_list(cls, bed_str_list, stranded=True):
        return BedFile([StrandedBedRecord.bed_record_from_str(line.strip()) for line in bed_str_list], stranded=True) \
            if stranded else \
            BedFile([UnstrandedBedRecord.bed_record_from_str(line.strip()) for line in bed_str_list], stranded=False)

    @property
    def bed_list(self):
        return self._bed_list

    def add_record(self, value):
        self._bed_list.append(value)

    def add_record_from_list(self, bed_list):
        self.add_record(self._record_class(bed_list))

    def add_record_from_str(self, bed_str):
        self.add_record(self._record_class.bed_record_from_str(bed_str))

    def add_record_list(self, rec_list):
        self._bed_list.extend(rec_list)
        
    def __repr__(self):
        return str([record.__repr__ for record in self._bed_list])

    def __add__(self, other):
        return BedFile(self.bed_list + other.bed_list)
    
    def __len__(self):
        return len(self.bed_list)
    
    def __getitem__(self, offset):
        return self.bed_list[offset]
    
    def __getslice__(self, low, high):
        return BedFile(self.bed_list[low:high], stranded=self._stranded)

    def keep_chr(self, chromosome):
        return BedFile([b for b in self.bed_list if b.name == chromosome])

    def keep_name(self, name):
        return BedFile([b for b in self.bed_list if b.name == name]) if self._stranded else self

    def keep_score(self, score):
        return BedFile([b for b in self.bed_list if b.score == score]) if self._stranded else self

    def keep_strand(self, strand):
        return BedFile([b for b in self.bed_list if b.strand == strand]) if self._stranded else self

    def bed_str_list(self):
        return [b.bed() for b in self.bed_list]

    def bed(self):
        return "\n".join(self.bed_str_list())

    def sort(self):
        ic(self._bed_list[:3])
        # sorts on chromosome, start, end
        self._bed_list.sort(key=lambda bed_record: (bed_record.chr, bed_record.start,bed_record.end))
         
    def keep_unique(self):
        # keep only unique records
        return BedFile.bed_file_from_str_list(bed_str_list=list(set(self.bed_str_list())))
    
    def nullify_name_and_score(self):
        # make self and score 0 to make life easy
        return BedFile([record.nullify_name_and_score() for record in self._bed_list], stranded=True) if self._stranded else self
        
    def unstrand_bed(self):
        return BedFile([record.unstrand() for record in self._bed_list], stranded=False) if self._stranded else self