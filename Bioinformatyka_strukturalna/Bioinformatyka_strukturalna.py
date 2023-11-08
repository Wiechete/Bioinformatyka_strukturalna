from enum import Enum
import numpy as np
import sys

if len(sys.argv) != 3:
    print("Zla komenda sprobuj: python program.py sekwencja1.fasta sekwencja2.fasta")
    sys.exit(1)

file_name1 = sys.argv[1]
file_name2 = sys.argv[2]

try:
# czytanie pliku
    with open(file_name1, 'r') as file1:
        next(file1)
        aminoacid1 = file1.read()

    with open(file_name2, 'r') as file2:
        next(file2)
        aminoacid2 = file2.read()

    print("Sekwencja 1:")
    print(aminoacid1)

    print("Sekwencja 2:")
    print(aminoacid2)
    
# error jesli nie ma pliku
except FileNotFoundError:
    print("Jeden z plikow nie istnieje.")

# wpisywanie dlugosci/inicjalizowanie
length1 = len(aminoacid1)
length2 = len(aminoacid2)
Matrix = np.zeros((length1 + 1, length2 + 1))

# Inicjalizacja macierzy traceback do sledzenia sciezki
Traceback = np.zeros((length1 + 1, length2 + 1), dtype=int)

# Inicjalizacja wartosci kary za luki (gap penalty)
gap_penalty = -2

for j in range(length1 + 1):
    Matrix[j][0] = gap_penalty * j

for j in range(length2 + 1):
    Matrix[0][j] = gap_penalty * j

# Wypelnienie macierzy wynikowej (algorytm Needlemana-Wunscha)
for i in range(1, length1 + 1):
    for j in range(1, length2 + 1):
        match = Matrix[i - 1][j - 1] + (1 if aminoacid1[i - 1] == aminoacid2[j - 1] else -1)
        delete = Matrix[i - 1][j] + gap_penalty
        insert = Matrix[i][j - 1] + gap_penalty

        # Wybor maksymalnej wartosci
        Matrix[i][j] = max(match, delete, insert)

        # sledzenie sciezki (backtracking)
        if Matrix[i][j] == match:
            Traceback[i][j] = 1  # Przesuniecie na ukos
        elif Matrix[i][j] == delete:
            Traceback[i][j] = 2  # Przesuniecie w gore
        else:
            Traceback[i][j] = 3  # Przesuniecie w lewo

# Wartosc optymalnego dopasowania znajduje sie w Matrix[length1][length2]

print("Matrix:")
for i in range(length1 + 1):
    for j in range(length2 + 1):
        print(Matrix[i][j], end="\t")
    print() 

print("Traceback:")
for i in range(length1 + 1):
    for j in range(length2 + 1):
        print(Traceback[i][j], end="\t")
    print() 
    
# Wyszukiwanie optymalnej sciezki (backtracking) w macierzy Traceback
i, j = length1, length2
aligned_seq1 = ""
aligned_seq2 = ""

while i > 0 and j > 0:
    if Traceback[i][j] == 1:  # Przesuniecie na ukos
        aligned_seq1 = aminoacid1[i - 1] + aligned_seq1
        aligned_seq2 = aminoacid2[j - 1] + aligned_seq2
        i -= 1
        j -= 1
    elif Traceback[i][j] == 2:  # Przesuniecie w gore
        aligned_seq1 = aminoacid1[i - 1] + aligned_seq1
        aligned_seq2 = "-" + aligned_seq2
        i -= 1
    else:  # Traceback[i][j] == 3, Przesuniecie w lewo
        aligned_seq1 = "-" + aligned_seq1
        aligned_seq2 = aminoacid2[j - 1] + aligned_seq2
        j -= 1

# W razie potrzeby doprowadz do konca wyrazow, jeeli i lub j nie wynosza 0
while i > 0:
    aligned_seq1 = aminoacid1[i - 1] + aligned_seq1
    aligned_seq2 = "-" + aligned_seq2
    i -= 1

while j > 0:
    aligned_seq1 = "-" + aligned_seq1
    aligned_seq2 = aminoacid2[j - 1] + aligned_seq2
    j -= 1

# Wyswietlenie wynikow
with open("wynik.txt", "w") as file:
    # Zapisanie wynikow do pliku
    file.write("Optimal Alignment:\n")
    file.write(aligned_seq1 + "\n")
    file.write(aligned_seq2 + "\n")
    file.write("Alignment Score: " + str(Matrix[length1][length2]) + "\n")
