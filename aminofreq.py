def codons_extract(DNA, frame): 
    codons =[]
    for i in range(int((len(DNA) - frame)/3)):
        codon = ""
        for j in range(3):
            codon += (DNA[3*i + j + frame])
        codons.append(codon)
        # codons.append(DNA[3*i + frame] + DNA[3*i + 1 + frame] + DNA[3*i + 2 + frame])
    return codons

print(codons_extract("ATAATGGATCGCCG", 4))

def finder(codon, start, stop):
    first = -1
    last = -1
    for i in range(len(codon)):
        if codon[i] in start:
            first = i
            break
    for i in range(first+1,len(codon)):
        if codon[i] in stop:
            last = i
            break
    return first,last

def protein_extract(codon,start,stop):
    protein = []
    first,last = finder(codon,start,stop)
    if first==-1 or last ==-1:
        return protein
    else:
        protein = codon[first : last+1]
    return protein

print(protein_extract(['CTA','ATG','GTG','GTC','ATG','AAT','TAG','GGT'],['GTG','ATG'],['TAG','TGA']))

# def protein_extract(codon, start, stop): Williams lösn på protein_extract
#     protein = []
#     a = 0
#     for i in range(len(codon)):
#         if codon[i] in start:
#             protein.append(codon[i])
#             break
#         else:
#              a +=1
    
#     for i in range(a+1,len(codon)):
#         if codon[i] in stop:
#             protein.append(codon[i])
#             break
#         else: 
#             protein.append(codon[i])
#     b = 0
#     for i in range(len(codon)):
#         if codon[i] in stop:
#             b = i
#     if b < a:
#         protein = []     
        
#     return protein

def aa_count(codon, dic):
    # codon = list(dict.fromkeys(codon)), if wanting to count duplicates only once: Uncomment this line.
    amino = {}
    for i in range(len(codon)): # i = 0, in range 1
        if codon[i] in dic:
            if dic[codon[i]] in amino:
                amino[dic[codon[i]]] += 1 
            else:
                amino[dic[codon[i]]] = 1
    print(amino)
    pass

aa_count(['TTT','TTT', 'BRR'], {'ATG':'M','TTT':'F','TTC':'F','GTC':'V','GTG':'V','AAC':'N','TAT':'Y'})     

def read_dna(ident, path):
    f = open(path)
    str = ""
    lista = f.read().splitlines()

    while '' in lista:
        lista.remove('')

    for i in range(len(lista)):
        lista[i] = lista[i].upper()
    
    ident = ident.upper()
    start = len(lista)
    stop = len(lista)
    for i in range(len(lista)):
        if ident in lista[i]:
            start = i+1
            break
    
    for i in range(start, len(lista)):
        if ((lista[i])[0]) != "A" and ((lista[i])[0]) != "T" and ((lista[i])[0]) != "C" and ((lista[i])[0]) != "G" :
            stop = i
            break
    for i in range(start, stop):
        str = str + lista[i]
    return str 

print(read_dna("sequence_3","Lab1/lab1/examples/example4.fna"))

def composition(ident, path, dic, start, stop):
    if read_dna(ident, path) == False:
        return
    for i in range(3):
        lista = protein_extract(codons_extract(read_dna(ident ,path),i), start, stop)
        if protein_extract(codons_extract(read_dna(ident ,path),i), start, stop) == []:
            print(f"Frame {i}: no valid protein")

        else:
            print(f"Frame {i}:")
            print(f"protein: {protein_extract(codons_extract(read_dna(ident ,path),i), start, stop)}")
            print(f"amino acid count: {aa_count(lista, dic)}")
genetic_code = {'GCT':'A','GCC':'A','GCA':'A','GCG':'A', 
               'CGT':'R','CGC':'R','CGA':'R','CGG':'R','AGA':'R','AGG':'R', 
               'AAT':'N','AAC':'N', 
               'GAT':'D','GAC':'D', 
               'TGT':'C','TGC':'C', 
               'CAA':'Q','CAG':'Q', 
               'GAA':'E','GAG':'E', 
               'GGT':'G','GGC':'G','GGA':'G','GGG':'G', 
               'CAT':'H','CAC':'H', 
               'ATT':'I','ATC':'I','ATA':'I', 
               'CTT':'L','CTC':'L','CTA':'L','CTG':'L','TTA':'L','TTG':'L', 
               'AAA':'K','AAG':'K', 
               'ATG':'M', 
               'TTT':'F','TTC':'F', 
               'CCT':'P','CCC':'P','CCA':'P','CCG':'P', 
               'TCT':'S','TCC':'S','TCA':'S','TCG':'S','AGT':'S','AGC':'S', 
               'ACT':'T','ACC':'T','ACA':'T','ACG':'T', 
               'TGG':'W', 
               'TAT':'Y','TAC':'Y', 
               'GTT':'V','GTC':'V','GTA':'V','GTG':'V'}