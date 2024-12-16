from TRgenerator import TR_multiMotif
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

##################################
# generate test data
##################################

test_data = TR_multiMotif(['AATGG','CCATT','CTT'], 13000, 0.1, seed)
test_seq = Seq(test_data.sequence)
test_record = SeqRecord(test_seq, id = "testSeq", description = '')

with open("testData/testData.fasta", "w") as handle:
    SeqIO.write(test_record, handle, "fasta")
df = test_data.annotation
df.to_csv('testData/testData_annotation.tsv', sep = '\t', index = False)
df = test_data.annotation_woMut
df.to_csv('testData/testData_annotation_woMut.tsv', sep = '\t', index = False)
