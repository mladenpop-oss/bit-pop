import sys

input_file = sys.argv[1] if len(sys.argv) > 1 else "D:/CAMI_low/cami_reads_100k_random.fastq"
r1_file = "D:/CAMI_low/cami_reads_R1.fastq"
r2_file = "D:/CAMI_low/cami_reads_R2.fastq"

r1_count = 0
r2_count = 0

with open(input_file, 'r') as fin, open(r1_file, 'w') as fout1, open(r2_file, 'w') as fout2:
    while True:
        header1 = fin.readline()
        if not header1:
            break
        seq1 = fin.readline()
        plus1 = fin.readline()
        qual1 = fin.readline()
        
        header2 = fin.readline()
        seq2 = fin.readline()
        plus2 = fin.readline()
        qual2 = fin.readline()
        
        if header1.strip().endswith('/1'):
            fout1.write(header1 + seq1 + plus1 + qual1)
            r1_count += 1
            fout2.write(header2 + seq2 + plus2 + qual2)
            r2_count += 1
        else:
            fout2.write(header1 + seq1 + plus1 + qual1)
            r2_count += 1
            fout1.write(header2 + seq2 + plus2 + qual2)
            r1_count += 1

print(f"R1: {r1_count} reads")
print(f"R2: {r2_count} reads")
