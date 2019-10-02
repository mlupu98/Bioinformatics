
import sys
import numpy as np
from collections import defaultdict

def referenceStr(reference):
    f = open(reference, 'r')
    first_line = True
    ref = ''
    for line in f:
        if first_line:
            first_line = False
            continue  # We skip the first line, since it
            # only contains the name of the chromosome the reads
            # came from.
        line = line.strip()
        ref += line  # We append each line to the output reference string.
    #ref +="AAAAAAAAAAAXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx"

    print(len(ref))
    return ref

def referenceDict(ref):
    refDict = defaultdict(lambda: [])

    newElem = ref[0:10]
    newRef = ref[10:]

    refDict[newElem].append(0)

    pos = 1
    for element in newRef:
        newElem = newElem[1:]
        newElem += element
        refDict[newElem].append(pos)
        pos += 1

    return refDict

# opens the reads file and place them in a list
def splitReads(reads):
    f = open(reads).read().splitlines()

    tupleLst = []

    firstLine = True

    for lines in f:
        if firstLine:
            firstLine = False
            continue

        tupleLst.append(lines.split(',')[0])
        tupleLst.append(lines.split(',')[1])

    return tupleLst

'''def checkMatch(refStr, startStr, read):

    retVal = False
    iStart = startStr

    quality = 0

    for i in range(0, 10): # this used to be 10 but doesnt work with 10, why?
        #print(i)
        #print(refStr[startStr], read[i])
        if refStr[iStart] == read[i]:
            quality += 1


    if quality >= 9: # can change this value to be more strict on the number of incorrect values
        retVal = True

    return retVal'''

def checkMatch(refStr, startStr, read):

    retVal = 0
    iStart = startStr

    quality = 0

    for i in range(0, 10): # this used to be 10 but doesnt work with 10, why?
        if refStr[iStart+i] == read[i]:
            quality += 1


    if quality >= 9: # can change this value to be more strict on the number of incorrect values
        retVal = 1

    return retVal

def reverseReads(reads):

    for elem in reads:
        reads.append(elem[::-1])


def checkDic(refStr, refDic, reads):

   qualityCheck = 0

   nucleotideArr = []

   for i in range(0, 9999):
       nucleotideArr.append([])

   counter = 0
   #newReads = reverseReads(reads)

   found = False


   counter = 0
   for readElem in reads:
       for i in range(0, 5):
           start = i * 10
           startRev = 40 - i * 10
           newRead = readElem[::-1]

           '''
           if readElem[start:start + 10] in refDic:
               counter += 1
               for elem in refDic[readElem[start:start+10]]:
                   readStart = elem - i * 10
                   for j in range(0, 5):
                       jump = j * 10
                       if checkMatch(refStr, readStart + jump, readElem[jump:jump + 10]):
                           for k in range(0, 10):
                               nucleotideArr[readStart + jump + k].append(readElem[jump + k])
                           refDic[readElem[start:start+10]].append(refDic[readElem[start:start+10]][0])
                           del refDic[readElem[start:start+10]][0]

           elif newRead[startRev:startRev+10] in refDic:
               counter += 1
               for elem in refDic[readElem[start:start+10]]:
                   readStart = refDic[newRead[startRev:startRev + 10]][0] - (40 - i * 10)
                   for j in range(0, 5):
                       jump = j * 10
                       if checkMatch(refStr, readStart + jump, newRead[jump:jump + 10]):
                           for k in range(0, 10):
                               nucleotideArr[readStart + jump + k].append(newRead[jump + k])
           break
           '''

           if readElem[start:start + 10] in refDic:
               counter += 1
               for elem in refDic[readElem[start:start + 10]]:
                   readStart = elem - i * 10
                   for j in range(0, 5):
                       jump = j * 10
                       if checkMatch(refStr, readStart + jump, readElem[jump:jump + 10]):
                           for k in range(0, 10):
                               nucleotideArr[readStart + jump + k].append(readElem[jump + k])
                           refDic[readElem[start:start + 10]].append(refDic[readElem[start:start + 10]][0])
                           del refDic[readElem[start:start + 10]][0]

           elif newRead[startRev:startRev + 10] in refDic:
               counter += 1
               for elem in refDic[readElem[start:start + 10]]:
                   readStart = refDic[newRead[startRev:startRev + 10]][0] - (40 - i * 10)
                   for j in range(0, 5):
                       jump = j * 10
                       if checkMatch(refStr, readStart + jump, newRead[jump:jump + 10]):
                           for k in range(0, 10):
                               nucleotideArr[readStart + jump + k].append(newRead[jump + k])
           break

   return nucleotideArr



#def placePairs(refStr, refDic, readsTuples):

 #   return True


def largest(A, C, G, T):

    retVal = 'X'

    if A == 0 and C == 0 and G == 0 and T == 0:
        retVal = 'X'
    elif A >= C and A >= G and A >= T:
        retVal = 'A'
    elif C > A and C > G and C > T:
        retVal = 'C'
    elif G > A and G > C and G > T:
        retVal = 'G'
    elif T > A and T > C and T > G:
        retVal = 'T'

    return retVal


def findCommon(nucleotideArr):

    commonArr = []
    A = 0
    C = 0
    G = 0
    T = 0

    for elem in nucleotideArr:
        for nucleo in elem:
            if nucleo == 'A':
                A += 1
            elif nucleo == 'C':
                C += 1
            elif nucleo == 'G':

                G += 1
            elif nucleo == 'T':
                T += 1
        commonArr.append(largest(A, C, G, T))
        A = 0
        C = 0
        G = 0
        T = 0

    return commonArr


def findSNP(refStr, nucleotideArr):

    commonArr = findCommon(nucleotideArr)

    SNP = []
    count = 0

    print(">SNP")
    for i in range(0, len(nucleotideArr)): #len(nucleotideArr)): should be common arr
        if commonArr[i] == 'X':
            continue
        if refStr[i] != commonArr[i]:
            count += 1
            newItem = [refStr[i], commonArr[i], i]
            SNP.append(newItem)

    realSNP = []

    for j in range(1, len(SNP)-1):
        if (SNP[j][2] - SNP[j-1][2]) < 3 or (SNP[j+1][2] - SNP[j][2]) < 3:
            print("HEY")
            continue
        else:
            realSNP.append(SNP[j])

    for i in range(0, len(SNP)):
        print(SNP[i][0] + ',' + SNP[i][1] + ',' + str(SNP[i][2]))

    #for i in range(0, len(realSNP)):
     #   print(realSNP[i][0] + ',' + realSNP[i][1] + ',' + str(realSNP[i][2]))


def main():
    reference = sys.argv[1]
    reads = sys.argv[2]

    refStr = referenceStr(reference)

    refDic = referenceDict(refStr)


    readLst = splitReads(reads) #this is a list

    nucleoArr = checkDic(refStr, refDic, readLst)

    #print(refDic)

    print(">hw1_W_2_chr_1")

    findSNP(refStr, nucleoArr)

    print(">STR")
    print(">CNV")
    print(">ALU")
    print(">INV")
    print(">INS")
    print(">DEL")

    return




if __name__ == "__main__":
    main()

