def read_dna(ident,path):
    file = open(path)
    file = file.readlines()
    return file

print(read_dna('Sequence1','examples/example1.fna'))