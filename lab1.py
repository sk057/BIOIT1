from Bio import SeqIO
record = SeqIO.read("data/bacterial4.fasta", "fasta")

#randa start ir stop kodonu indexus ir kodonu daznius
def readSequence(frame, seq, dictCodons, totalCodons):
  starts= []
  stops = []
  for i in range(frame, len(seq), 3):
    totalCodons += 1
    codon = record.seq[i:i+3]
    if codon == "ATG":   #start kodonu indexai
      starts.append(i)
    if codon == "TAG" or codon == "TAA" or codon == "TGA":  #stop kodonu indexai
      stops.append(i)
    if (len(codon)==3):
      if codon in dictCodons:
          dictCodons[codon] += 1
      else:
          dictCodons[codon] = 1
  return starts, stops, dictCodons, totalCodons

def readSequenceDicodons(frame, seq, dictDicodons, totalDicodons):
  for i in range(frame, len(seq), 6):
    totalDicodons += 1
    dicodon = record.seq[i:i+6]
    if (len(dicodon)==6):
      if dicodon in dictDicodons:
          dictDicodons[dicodon] += 1
      else:
          dictDicodons[dicodon] = 1
  return dictDicodons, totalDicodons

#1. Pateiktoje sekoje surastu visas start ir stop kodonų poras, tarp kurių nebutu stop kodono
def findPairs(starts, stops):
  startStopPairs = []
  pairFound = False
  lastStopIndex = -1
  for startIndex in starts:
    pairFound = False
    for stopIndex in stops:
      if (startIndex<stopIndex and not pairFound
      and startIndex>lastStopIndex and stopIndex>lastStopIndex #2. Kiekvienam stop kodonui parinkti toliausiai nuo jo esanti start kodoną
      and stopIndex-startIndex>100): #3. Atfiltruokite visus fragmentus ("tai butu baltymų koduojancios sekos"), kurie trumpesni nei 100 simboliu.
        startStopPairs.append((startIndex,stopIndex))
        lastStopIndex = stopIndex
        pairFound = True
  return startStopPairs

#5. Parasykite funkcijas, kurios ivertintu kodonu ir dikodonu daznius
def findFrequency(dictCodons, totalCodons, dictDicodons, totalDicodons):
  for codon in dictCodons:
    dictCodons[codon] = dictCodons[codon] / totalCodons
  for dicodon in dictDicodons:
    dictDicodons[dicodon] = dictDicodons[dicodon] / totalDicodons
  return dictCodons, dictDicodons

totalCodons = 0
totalDicodons = 0
dictCodons = {}
dictDicodons = {}
# 6 skaitymo remeliai
starts1, stops1, dictCodons, totalCodons = readSequence(0,record.seq, dictCodons, totalCodons)
startStopPairs1 = findPairs(starts1, stops1)
dictDicodons, totalDicodons = readSequenceDicodons(0, record.seq, dictDicodons, totalDicodons)

starts2, stops2, dictCodons, totalCodons = readSequence(1,record.seq, dictCodons, totalCodons)
startStopPairs2 = findPairs(starts2, stops2)
dictDicodons, totalDicodons = readSequenceDicodons(1, record.seq, dictDicodons, totalDicodons)

starts3, stops3, dictCodons, totalCodons = readSequence(2,record.seq, dictCodons, totalCodons)
startStopPairs3 = findPairs(starts3, stops3)
dictDicodons, totalDicodons = readSequenceDicodons(2, record.seq, dictDicodons, totalDicodons)

starts4, stops4, dictCodons, totalCodons = readSequence(0,record.seq.reverse_complement(), dictCodons, totalCodons)
startStopPairs4 = findPairs(starts4, stops4)
dictDicodons, totalDicodons = readSequenceDicodons(0, record.seq.reverse_complement(), dictDicodons, totalDicodons)

starts5, stops5, dictCodons, totalCodons = readSequence(1,record.seq.reverse_complement(), dictCodons, totalCodons)
startStopPairs5 = findPairs(starts5, stops5)
dictDicodons, totalDicodons = readSequenceDicodons(1, record.seq.reverse_complement(), dictDicodons, totalDicodons)

starts6, stops6, dictCodons, totalCodons = readSequence(2,record.seq.reverse_complement(), dictCodons, totalCodons)
startStopPairs6 = findPairs(starts6, stops6)
dictDicodons, totalDicodons = readSequenceDicodons(2, record.seq.reverse_complement(), dictDicodons, totalDicodons)

dictCodons, dictDicodons = findFrequency(dictCodons, totalCodons, dictDicodons, totalDicodons)

# print frequencies to file (through terminal)
# for codon in dictCodons:
#   print(codon, dictCodons[codon])
for dicodon in dictDicodons:
  print(dicodon, dictDicodons[dicodon])

# (258, 780)
# (300, 780) skip
# (609, 780) skip
# (612, 780) skip
# (783, 2403)
