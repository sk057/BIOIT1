#skaiciuoja atstuma tarp dvieju seku skirtumu moduliu budu
# ...+|tikimybe1*100proc-tikimybe2*100proc|+... kiekvienam kodonui

def readFile (location,l):
    freqDict = {}
    file = open(location, mode='r', encoding='utf-16')
    for line in file:
        line = line.strip()
        freqDict[line[0:l]] = round(float(line[l+2:len(line)]),5)
    return freqDict

freq1codon = readFile('frequencies/bacterial3codons.txt',3)
freq1dicodon = readFile('frequencies/bacterial3dicodons.txt',6)
freq2codon = readFile('frequencies/bacterial4codons.txt',3)
freq2dicodon = readFile('frequencies/bacterial4dicodons.txt',6)

sum = 0

for codon1 in freq1codon:
    if codon1 in freq2codon:
        sum += abs(freq1codon[codon1]*100-freq2codon[codon1]*100)
    else: sum+=freq1codon[codon1]*100

for codon2 in freq2codon:
    if codon2 not in freq1codon:
        sum -= freq2codon[codon1]*100

sumDicodon = 0

for dicodon1 in freq1dicodon:
    if dicodon1 in freq2dicodon:
        sumDicodon += abs(freq1dicodon[dicodon1]*100-freq2dicodon[dicodon1]*100)
    else: sumDicodon+=freq1dicodon[dicodon1]*100

for dicodon2 in freq2dicodon:
    if dicodon2 not in freq1dicodon:
        sumDicodon -= freq2dicodon[dicodon2]*100

print(round(sum,5))
print(round(sumDicodon,5))

# for codon in freq1dicodon:
#     print(codon, freq1dicodon[codon])