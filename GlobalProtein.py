from Bio.SubsMat import MatrixInfo
blosum = MatrixInfo.blosum62

first_protein_seq = input("Enter the 1st protein sequence: ")
second_protein_seq = input("Enter the 2nd protein sequence: ")
gap = int(input("Enter the gap penalty: "))

def diagonalmatch(pair, matrix):
    if pair not in matrix:
        return matrix[(tuple(reversed(pair)))]
    else:
        return matrix[pair]


matrix = [[0 for j in range(len(first_protein_seq) + 1)] for i in range(len(second_protein_seq) + 1)]
notations = [['' for j in range(len(first_protein_seq) + 1)] for i in range(len(second_protein_seq) + 1)]
diagonal = 0
for i in range(len(first_protein_seq) + 1):
    matrix[0][i] = i * gap
for i in range(len(second_protein_seq) + 1):
    matrix[i][0] = i * gap
for i in range(1, len(second_protein_seq) + 1):
    for j in range(1, len(first_protein_seq) + 1):
        diagonal = diagonalmatch((first_protein_seq[j - 1], second_protein_seq[i - 1]), blosum)
        # VDSCY   VESLCY   -8
        maximum = max(matrix[i - 1][j] + gap, matrix[i][j - 1] + gap, matrix[i - 1][j - 1] + diagonal)
        matrix[i][j] = maximum
        if maximum == matrix[i - 1][j] + gap:
            notations[i][j] = '|'
        elif maximum == matrix[i][j - 1] + gap:
            notations[i][j] = '-'
        elif maximum == matrix[i - 1][j - 1] + diagonal:
            notations[i][j] = '\\'

fs = ""
ss = ""
i = len(second_protein_seq)
j = len(first_protein_seq)

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
print(notations)
print(ss[::-1])
print(fs[::-1])
