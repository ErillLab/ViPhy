# This is a sample Python script.


def nucleotideToAminoacid(DNAseq):
    """
    Translate a nucleotide sequence into six amino acid sequences.
    DNAseq: nucleotide sequence, should be string
    Returns a single amino acid sequence that come from combining the six amino acid sequences, should be string
    """
    translationTable = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }

    proteinSeq1 = ""
    proteinSeq2 = ""
    proteinSeq3 = ""
    proteinSeq4 = ""
    proteinSeq5 = ""
    proteinSeq6 = ""

    # First direction
    for i in range(0, len(DNAseq), 3):
        codon1 = DNAseq[i:i + 3]
        codon2 = DNAseq[i + 1:i + 4]
        codon3 = DNAseq[i + 2:i + 5]
        if (len(codon1)==3):
            proteinSeq1 += translationTable[codon1]
        if (len(codon2)==3):
            proteinSeq2 += translationTable[codon2]
        if (len(codon3)==3):
            proteinSeq3 += translationTable[codon3]

    # Second direction
    DNAseqR = revertSeq(DNAseq)

    for j in range(0, len(DNAseqR), 3):
        codon1 = DNAseqR[j:j + 3]
        codon2 = DNAseqR[j + 1:j + 4]
        codon3 = DNAseqR[j + 2:j + 5]
        if (len(codon1) == 3):
            proteinSeq4 += translationTable[codon1]
        if (len(codon2) == 3):
            proteinSeq5 += translationTable[codon2]
        if (len(codon3) == 3):
            proteinSeq6 += translationTable[codon3]

    proteinSeq = proteinSeq1 + proteinSeq2 + proteinSeq3 + proteinSeq4 + proteinSeq5 + proteinSeq6

    return proteinSeq

def revertSeq(DNAseq):
    """
    Translate a nucleotide sequence into six amino acid sequences.
    DNAseq: nucleotide sequence, should be string
    Returns the opposite sequence to the input sequence, should be a string
    """
    reversSeq= DNAseq[::-1]
    rDNAseq = ""
    for i in reversSeq:
        if i == 'A':
            rDNAseq += 'T'
        elif i == 'C':
            rDNAseq += 'G'
        elif i == 'G':
            rDNAseq += 'C'
        elif i == 'T':
            rDNAseq += 'A'
        else:
            rDNAseq += i

    return rDNAseq

def fastaToString(DNAfile):
    """
    Takes a FASTA file with a nucleotide sequence and transforms it to a string.
    DNAfile: file containing nucleotide sequence, should be FASTA format
    Returns the nucleotide sequence obtained from the fasta file and the first line of the fasta file.
    """
    file = open(DNAfile)
    fastaContent = file.read()
    fastaList = fastaContent.split("\n")
    fastaName = fastaList[0]
    fastaList = fastaList[1:]
    DNAseq = ""
    for line in fastaList:
        DNAseq += line

    file.close()
    return DNAseq, fastaName

def exportFasta(proteinSeq, fastaName):
    """
    Saves the protein sequence into a fasta file
    proteinSeq: combination of the six possible protein sequences. Obtained from "nucleotideToAminoacid" function
    fastaName: first line from the original fasta file
    """
    fileContent = ""
    file = open("proteinSeq.fasta", "w")
    for i in proteinSeq:
        fileContent += i
        if (len(fileContent) % 63) == 0:
            fileContent += '\n'
    fileContent = fastaName + " " + fileContent
    file.write(fileContent)
    file.close()


# Main function
if __name__ == '__main__':
    textFile = 'sequence.fasta'
    DNAseq, fastaName = fastaToString(textFile)
    proteinSeq = nucleotideToAminoacid(DNAseq)
    exportFasta(proteinSeq, fastaName)



