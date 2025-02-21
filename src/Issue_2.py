import RNA
# me la he inventado
seq = "AUGUAGUAGUAUCUGUACUAGUGUAUGAUUCUAUGAUG"

(ss, mfe) = RNA.fold(seq)

print("{}\n{} [ {:6.2f} ]".format(seq, ss, mfe))