import sys
sys.path.insert(0, r'C:/Users/1/PycharmProjects/3DPredictor-tests/3Dpredictor/source')
from RNASeqReader import RNAseqReader
from ChiPSeqReader import ChiPSeqReader
from shared import Interval

# rna_reader = RNAseqReader(r"C:\Users\admin\3DPredictor-tests\3Dpredictor\input\K562\RNA-seq\rna-seqPolyA.tsvpre.txt", name="RNA")
rna_reader = RNAseqReader(r"C:\Users\1\PycharmProjects\3DPredictor-tests\3Dpredictor\input\K562\RNA-seq\rna-seqPolyA.tsvpre.txt", name="RNA")
rna_reader.read_file(rename={"Gene stable ID": "gene",
                                              "Gene start (bp)": "start",
                                              "Gene end (bp)": "end",
                                              "Chromosome/scaffold name": "chr",
                                              "FPKM": "sigVal"}, encoding='utf-8', sep='\t')
print("Data before inversion:")
print(rna_reader.chr_data['chr11'].iloc[10:20, [rna_reader.chr_data["chr11"].columns.get_loc("chr"),
                                                rna_reader.chr_data["chr11"].columns.get_loc("start"),
                                                rna_reader.chr_data["chr11"].columns.get_loc("end"),
                                                rna_reader.chr_data["chr11"].columns.get_loc("sigVal"),
                                                rna_reader.chr_data["chr11"].columns.get_loc("gene")]])
inversion = Interval("chr11", 207600, 280000)     # 1) 207501, 207505  2) 207500, 207600  3) 207600, 209700  4) 207500, 215200
# rna_reader.inverse_region_RNA(inversion, r"C:\Users\admin\3DPredictor-tests\tests\deletion_RNASeq_2\mart_export_txt")   # 5) 207500, 236600  6) 207500, 253000
# 7) 204623, 235000 8) 207600, 280000
rna_reader.inverse_region_RNA(inversion, r"C:\Users\1\PycharmProjects\3DPredictor-tests\3Dpredictor\tests\tests\deletion_RNASeq_2\mart_export_txt")
print("Inversion interval:", inversion.start - inversion.end)
print('start',inversion.start, '\n', 'end',inversion.end )
print("Data after inversion:")
print(rna_reader.chr_data['chr11'].iloc[10:20, [rna_reader.chr_data["chr11"].columns.get_loc("chr"),
                                                rna_reader.chr_data["chr11"].columns.get_loc("start"),
                                                rna_reader.chr_data["chr11"].columns.get_loc("end"),
                                                rna_reader.chr_data["chr11"].columns.get_loc("sigVal"),
                                                rna_reader.chr_data["chr11"].columns.get_loc("gene")]])
