import warnings
from Bio import BiopythonDeprecationWarning

warnings.simplefilter('ignore', BiopythonDeprecationWarning)
from Bio.SubsMat.MatrixInfo import blosum62
import itertools

# take query sequence
print("enter the query sequence")
query = input()

# word length
print("enter the word length")
word_len = int(input())


# word threshold
# print("enter the word threshold")
# word_thres = int(input())

# hsp threshold
# print("enter the hsp threshold")
# hsp_thres = int(input())

# step 1 remove low complexity and repeated sequence
def remove_repeats():
    count = 0
    new_seq = ""
    seq1 = query[0]
    seq2 = query[1]
    var1 = ''
    var2 = ''
    # QCEgcgcgcGHI
    for i in range(2, len(query)):
        if query[i] == seq1 and query[i + 1] == seq2:
            count += 1
            var1 = seq1
            var2 = seq2
        seq1 = seq2
        seq2 = query[i]
    if count >= 3:
        for i in range(len(query)):
            if query[i] != var1 and query[i] != var2:
                new_seq += query[i]
        return new_seq
    else:
        return query


print("sequence without repeats " + remove_repeats())

# open file
databasefile = open('database_sequence.txt')

# step 2 splitting sequence into words
new_query = remove_repeats()
word = []
for i in range(0, len(new_query)):
    if i + word_len - 1 < len(new_query):
        word.append(new_query[i:i + word_len])
print("words ", end="")
print(word)


# calculate score
def cal_score(seq1, y, matrix):
    score = 0
    for i in range(len(seq1)):
        seq = (seq1[i], y[i])
        if seq not in blosum62:
            reverse_pair = tuple(reversed(seq))
            score += blosum62[reverse_pair]
        else:
            score += blosum62[seq]
    return score


# step 3 neighbourhood words
def neighbourwords():
    AminoAcid = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V',
                 'B', 'Z', 'X']
    str2 = ""
    score = 0
    # blosum62_dic = blosum62
    compination_dic = {}
    for i in range(len(word)):
        for j in range(0, 2):
            for l in range(len(AminoAcid)):
                if j == 0:
                    str2 = word[i]
                    str2 = AminoAcid[l] + str2[1] + str2[2]
                    score = cal_score(word[i], str2, blosum62)
                    compination_dic[str2] = score
                elif j == 1:
                    str2 = word[i]
                    str2 = str2[0] + AminoAcid[l] + str2[2]
                    score = cal_score(word[i], str2, blosum62)
                    compination_dic[str2] = score
                elif j == 2:
                    str2 = word[i]
                    str2 = str2[0] + str2[1] + AminoAcid[l]
                    score = cal_score(word[i], str2, blosum62)
                    compination_dic[str2] = score
    return compination_dic


print(neighbourwords())


# filter seeds for specfic threshold
def threshold(t=17):
    seeds = {}
    neighwords = neighbourwords()
    key = list(neighwords.keys())
    value = list(neighwords.values())
    for i in range(len(key)):
        if value[i] >= t:
            seeds[key[i]] = value[i]
    return seeds


print(threshold())


def exactmatch():
    databaseseq = databasefile
    word_hit = {}
    seeds = threshold()
    key = list(seeds.keys())
    for seed in key:
        list_Seq = []
        for seq in databaseseq:
            if seed in seq:
                list_Seq.append(seq[:-2])
        word_hit[seed] = list_Seq
        databaseseq.seek(0)
        if word_hit[seed] == []:
            del word_hit[seed]

    return word_hit


print(exactmatch())


# extend seeds
def extention(hspthres=15):
    database = databasefile.readlines()
    word_exact_hit = exactmatch()
    key = word_exact_hit.keys()
    value = word_exact_hit.values()
    seed_score = threshold()
    for hit in len(query):
        for j in range(word_exact_hit, len(value)):

            for wordhit in key:
                max_hits = []
                current_score = 0
            max_score = 0
            # extend to right
            for j in range(len(query)):
                current_score += cal_score(query[j], database[wordhit], blosum62)
                wordhit = wordhit + 1
            if current_score >= max_score:
                max_score = max_score + current_score
            elif max_score - current_score > hspthres:
                break

                # extend to left
                current_score = max_score
                for k in range(wordhit - 1, -1, -1):
                    # blossom
                    current_score += cal_score(query[k], database[wordhit], blosum62)
                    wordhit = wordhit - 1
                if current_score >= max_score:
                    max_score = max_score + cal_score(query[k], database[wordhit[2] + wordhit[1]])
                elif max_score - current_score > hspthres:
                    break
                max_hits.append(max_score)

    return max(max_hits)

print(extention(15))
# close file
databasefile.close()
