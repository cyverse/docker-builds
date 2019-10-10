#!/usr/bin/python


def fasta_reader(infile):
        read_file = open(infile, "rU")
        seq_id, seq = None, []
        for line in read_file:
            line = line.rstrip()
            if line.startswith(">"):
                if seq_id:
                    yield (seq_id, ''.join(seq))
                seq_id, seq = line, []
            else:
                seq.append(line)
        if seq_id:
            yield (seq_id, ''.join(seq))


def qual_filter(MinQual, QThr, ReadQual, QualFact):
    ReadQual = list(ReadQual)
    ReadQualNum = map(ord, ReadQual)
    ReadQualNum = [ele-QualFact for ele in ReadQualNum]
    LowQualCount = sum(i <= QThr for i in ReadQualNum)
    ReadLen = len(ReadQual)
    result = (float(LowQualCount)*100)/float(ReadLen)
    if result >= MinQual:
        return True
    else:
        return False


def success_flag(success_path):
        create_file = open(success_path, 'w')
        create_file.write("success")
        create_file.close()