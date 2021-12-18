firstString = "ACCGTT"
secondString = "AGTTCA"
gap = -1
mm = -1
m = 1
matrix = [[0 for j in range(len(firstString) + 1)] for i in range(len(secondString) + 1)]
notations = [['' for j in range(len(firstString) + 1)] for i in range(len(secondString) + 1)]
score = 0
diagonal = 0
for i in range(len(firstString) + 1):
    matrix[0][i] = i*-1
for i in range(len(secondString) + 1):
    matrix[i][0] = i*-1
for i in range(1, len(secondString) + 1):
    for j in range(1, len(firstString) + 1):
        if firstString[j - 1] == secondString[i - 1]:
            score = 1
            diagonal = m
        else:
            score = 0
            diagonal = mm
        maximum = max(matrix[i - 1][j] + gap, matrix[i][j - 1] + gap, matrix[i - 1][j - 1] + diagonal)
        matrix[i][j] = maximum
        if maximum == matrix[i - 1][j] + gap:
            notations[i][j] = '|'
        elif maximum == matrix[i][j - 1] + gap:
            notations[i][j] = '-'
        elif maximum == matrix[i - 1][j - 1] + diagonal:
            notations[i][j] = '\\'
fs=""
ss=""
i = len(secondString)
j = len(firstString)
while i > 0 and j > 0:
    if(notations[i][j] == "\\"):
        i = i-1
        j = j-1
        fs = fs + firstString[j]
        ss = ss + secondString[i]
    elif notations[i][j] == '|':
        i = i-1
        fs = fs + '-'
        ss = ss + secondString[i]
    elif notations[i][j] == '-':
        j = j - 1
        fs = fs + firstString[j]
        ss = ss + '-'
while i > 0:
    i = i-1
    fs = fs + '-'
    ss = ss + secondString[i]
while j > 0:
    j = j - 1
    fs = fs + firstString[j]
    ss = ss + '-'
print(matrix)
print(notations)
print(ss[::-1])
print(fs[::-1])