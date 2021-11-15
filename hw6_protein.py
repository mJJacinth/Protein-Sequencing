"""
Protein Sequencing Project
Name:
Roll Number:
"""

from os import read, remove
import hw6_protein_tests as test
from textwrap import wrap
from itertools import zip_longest

project = "Protein" # don't edit this

### WEEK 1 ###

'''
readFile(filename)
#1 [Check6-1]
Parameters: str
Returns: str
'''
def readFile(filename):
    str=""
    f=open(filename,"r")
    text=f.read().splitlines()
    for i in text:
        str+=i
    return str


'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):
    lst=[]
    rna=dna.replace("T","U")
    for i in range(startIndex,len(dna),3):
        temp=rna[i:i+3]
        if rna[i:i+3]=="UAA" or rna[i:i+3]=="UAG" or rna[i:i+3]=="UGA":
            lst.append(temp)
            break
        else:
            lst.append(temp)
    return lst     

'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
def makeCodonDictionary(filename):
    dictionary={}
    lst=[]
    import json
    f=open(filename,"r")
    file=json.load(f)
    for i in file:
        for j in file[i]:
            new=j.replace("T","U")
            dictionary[new]=i
    return dictionary
        

'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''
def generateProtein(codons, codonD):
    lst=[]
    for i in codons:
        if i  in codonD:
            lst.append(codonD[i])
            if codonD[i]=="Stop":
                break
    if lst[0] !="Start":
        lst[0]="Start"        
    return lst


'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):
    dna_file=readFile(dnaFilename)
    codon_file=makeCodonDictionary(codonFilename)
    lst=[]
    count=0 
    i=0
    while i!= len(dna_file):
        code=dna_file[i:i+3]
        if code=="ATG":
            rna=dnaToRna(dna_file,i)
            protein=generateProtein(rna,codon_file)
            lst.append(protein)
            i=i+(3*len(rna))
        else:
            i+=1
            count+=1
    return lst
    
def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")


### WEEK 2 ###

'''
commonProteins(proteinList1, proteinList2)
#1 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs
Returns: 2D list of strs
'''
def commonProteins(proteinList1, proteinList2):
    common_pro=[]
    for i in proteinList1:
        if i in proteinList2:
            common_pro.append(i)
    return common_pro


'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''
def combineProteins(proteinList):
    combine_pro=[]
    for i in proteinList:
        for j in i:
            combine_pro.append(j)
    return combine_pro


'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):
    amino_acids_dict={}
    for i in aaList:
        if i not in amino_acids_dict:
            amino_acids_dict[i]=1
        else:
            amino_acids_dict[i]+=1
    return amino_acids_dict


'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    lst=[]
    list1=combineProteins(proteinList1)
    list2=combineProteins(proteinList2)
    dic1=aminoAcidDictionary(list1)
    dic2=aminoAcidDictionary(list2)
    diff_lis1=list(set(list2)-set(list1))
    diff_lis2=list(set(list1)-set(list2))
    for i in diff_lis1:
        dic1[i]=0
    for j in diff_lis2:
        dic2[j]=0
    len1=len(list1)
    freq1={}
    for acid1 in dic1:
        freq1[acid1]=dic1[acid1]/len1
    len2=len(list2)
    freq2={}
    for acid2 in dic2:
        freq2[acid2]=dic2[acid2]/len2
    for k,l in dic1.items():
        lst2=[]
        if k!="Start" and k!="Stop":
            if k in dic2.keys():
                freq=abs(freq1[k]-freq2[k])
                if freq>cutoff:
                    lst2.append(k)
                    lst2.append(freq1[k])
                    lst2.append(freq2[k])
                    lst.append(lst2)         
    return lst
'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    print("The commonalities:")
    for i in sorted(commonalities):
        commonProteins=""
        let =i[1:len(i)-1]
        count=0
        for j in let:
            commonProteins+=j
            count+=1
            if count!=len(let):
                commonProteins+="-"
    print(commonProteins)
    for i in commonalities:
        if i=="Start" or i=="Stop":
            pass
        else:
            print(i)
    print("DNA sequences")
    for item in differences:
        print(item[0],":",round(item[1]*100,2),"% in seq1,",round(item[2]*100,2),"% in seq2")
    return       


def runWeek2():
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)


### WEEK 3 ###

'''
makeAminoAcidLabels(proteinList1, proteinList2)
#2 [Hw6]
Parameters: 2D list of strs ; 2D list of strs
Returns: list of strs
'''
def makeAminoAcidLabels(proteinList1, proteinList2):
    pro_lis1,pro_lis2=combineProteins(proteinList1),combineProteins(proteinList2)
    lst=[]
    for i,j in zip_longest(pro_lis1,pro_lis2):
        if i not in lst and i!=None:
            lst.append(i)
        if j not in lst and j!= None:
            lst.append(j)
    return sorted(lst)


'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
    comb_list = combineProteins(proteinList) 
    pro_dict= aminoAcidDictionary(comb_list) 
    freq_list=[] 
    for i in labels:
        if i in pro_dict:
            freq_list.append(pro_dict[i]/ len(comb_list))
        else: 
            freq_list.append(0)  
    return freq_list
    


'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
    import matplotlib.pyplot as plt
    return


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):
    return


'''
runFullProgram()
#6 [Hw6]
Parameters: no parameters
Returns: None
'''
def runFullProgram():
    return


### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":
    # print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    # test.week1Tests()
    # print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    # runWeek1()

    ## Uncomment these for Week 2 ##
   
    # print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    # test.week2Tests()
    # print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    # runWeek2()
   

    ## Uncomment these for Week 3 ##
    print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    test.week3Tests()
    print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()
    
