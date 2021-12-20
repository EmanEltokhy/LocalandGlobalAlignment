# firstString = "ACCGTT"
# secondString = "AGTTCA"
import sys
from Bio.SubsMat import MatrixInfo
blosum = MatrixInfo.blosum62
Nucleotides = ["A", "C", "G", "T"]
def optionality():
    print("Welcome to our Local Alignment System")
    print("Do you want to make local alignment for 1) DNA sequences or 2) Protein sequences or 3) Exit the system?")
    choose = int(input("Enter your choice, please: "))
    if choose == 1:
        check_input()
    elif choose == 2:
        localProtein()
    else:
        sys.exit()

def diagonalmatch(pair, matrix):
    if pair not in matrix:
        return matrix[(tuple(reversed(pair)))]
    else:
        return matrix[pair]

def validateDNA(DNAseq):
    sequence = DNAseq.upper()
    for n in sequence:
        if n not in Nucleotides:
            return False
    return True

def check_input():
    firstString = input("Enter the 1st DNA sequence: ").upper()
    secondString = input("Enter the 2nd DNA sequence: ").upper()
    if validateDNA(firstString) and validateDNA(secondString):
        localDNA(firstString, secondString)
    else:
        print("You've entered an invalid sequence, please try again :) ")
        check_input()

def localDNA(firstString, secondString):
    gap = -1
    mm = -2
    m = 1
    matrix = [[0 for j in range(len(firstString) + 1)] for i in range(len(secondString) + 1)]
    notations = [['' for j in range(len(firstString) + 1)] for i in range(len(secondString) + 1)]
    score = 0
    diagonal = 0
    for i in range(len(firstString) + 1):
        matrix[0][i] = 0
    for i in range(len(secondString) + 1):
        matrix[i][0] = 0
    for i in range(1, len(secondString) + 1):
        for j in range(1, len(firstString) + 1):
            if firstString[j - 1] == secondString[i - 1]:
                diagonal = m
            else:
                diagonal = mm
            maximum = max(matrix[i - 1][j] + gap, matrix[i][j - 1] + gap, matrix[i - 1][j - 1] + diagonal)
            if maximum > 0:
                matrix[i][j] = maximum
                if maximum > score:
                    score = maximum
                    dim = [i, j]
            else:
                matrix[i][j] = 0
            if maximum == matrix[i][j - 1] + gap:
                notations[i][j] = '-'
            elif maximum == matrix[i - 1][j] + gap:
                notations[i][j] = '|'

            elif maximum == matrix[i - 1][j - 1] + diagonal:
                notations[i][j] = '\\'
    fs = ""
    ss = ""
    i = dim[0]
    j = dim[1]
    while i > 0 and j > 0:
        if (notations[i][j] == "\\"):
            i = i - 1
            j = j - 1
            fs = fs + firstString[j]
            ss = ss + secondString[i]
        elif notations[i][j] == '|':
            i = i - 1
            fs = fs + '-'
            ss = ss + secondString[i]
        elif notations[i][j] == '-':
            j = j - 1
            fs = fs + firstString[j]
            ss = ss + '-'
    while i > 0:
        i = i - 1
        fs = fs + '-'
        ss = ss + secondString[i]
    while j > 0:
        j = j - 1
        fs = fs + firstString[j]
        ss = ss + '-'
    print(matrix)
    for i in range(len(notations)):
        for j in range(len(notations[0])):
            print(notations[i][j], end='')
        print()
    print(ss[::-1])
    print(fs[::-1])

def localProtein():
    first_protein_seq = input("Enter the 1st protein sequence: ").upper()
    second_protein_seq = input("Enter the 2nd protein sequence: ").upper()
    gap = int(input("Enter the gap penalty: "))
    matrix = [[0 for j in range(len(first_protein_seq) + 1)] for i in range(len(second_protein_seq) + 1)]
    notations = [['' for j in range(len(first_protein_seq) + 1)] for i in range(len(second_protein_seq) + 1)]
    score = 0
    diagonal = 0
    for i in range(len(first_protein_seq) + 1):
        matrix[0][i] = 0
    for i in range(len(second_protein_seq) + 1):
        matrix[i][0] = 0
    for i in range(1, len(second_protein_seq) + 1):
        for j in range(1, len(first_protein_seq) + 1):
            diagonal = diagonalmatch((first_protein_seq[j - 1], second_protein_seq[i - 1]), blosum)
            # VDSCY   VESLCY   -8
            maximum = max(matrix[i - 1][j] + gap, matrix[i][j - 1] + gap, matrix[i - 1][j - 1] + diagonal)
            if maximum > 0:
                matrix[i][j] = maximum
                if maximum > score:
                    score = maximum
                    dim = [i, j]
            else:
                matrix[i][j] = 0
            if maximum == matrix[i][j - 1] + gap:
                notations[i][j] = '-'
            elif maximum == matrix[i - 1][j] + gap:
                notations[i][j] = '|'

            elif maximum == matrix[i - 1][j - 1] + diagonal:
                notations[i][j] = '\\'
    fs = ""
    ss = ""
    i = dim[0]
    j = dim[1]
    while i > 0 and j > 0:
        if (notations[i][j] == "\\"):
            i = i - 1
            j = j - 1
            fs = fs + first_protein_seq[j]
            ss = ss + second_protein_seq[i]
        elif notations[i][j] == '|':
            i = i - 1
            fs = fs + '-'
            ss = ss + second_protein_seq[i]
        elif notations[i][j] == '-':
            j = j - 1
            fs = fs + first_protein_seq[j]
            ss = ss + '-'
    while i > 0:
        i = i - 1
        fs = fs + '-'
        ss = ss + second_protein_seq[i]
    while j > 0:
        j = j - 1
        fs = fs + first_protein_seq[j]
        ss = ss + '-'
    print(matrix)
    for i in range(len(notations)):
        for j in range(len(notations[0])):
            print(notations[i][j], end='')
        print()
    print(ss[::-1])
    print(fs[::-1])

optionality()