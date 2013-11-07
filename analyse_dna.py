__author__ = 'markus'

import timeit


def read_fasta(filename):
    """ Reads a sequence in Fasta format """
    with open(filename, 'r') as fp:
        header = ""
        seq = ""
        for line in fp:
            if line == "":
                break
            if line.startswith('>'):
                header = line[1:].strip()
            else:
                seq += line.strip().upper()
    return header, seq


def calc_unique_kmer_len(seq_len):
    residual = seq_len / 4
    kmer_len = 1
    while residual >= 1.0:
        residual /= 16
        kmer_len += 2
    if residual * 4 <= 1.0:
        kmer_len -= 1
    return kmer_len


def find_kmers_pos(seq, kmer_len):
    seq_len = len(seq)
    kmers_pos = {}

    for i in range(seq_len - len(kmer_len)):
        curr_kmer = seq[i:i+kmer_len]
        if 'N' in curr_kmer:
            continue

        kmers_pos.get(curr_kmer, []).append(i)

        #if curr_kmer in unique_kmers:
        #    unique_kmers[curr_kmer][0][0] += 1
        #    unique_kmers[curr_kmer][1].append(i)
        #    non_unique += 1
        #else:
        #    unique_kmers[curr_kmer] = ([1], [i])
    #print(non_unique)
    kmers_pos = {k: v for k, v in kmers_pos.items() if len(v) > 2}
    return kmers_pos


def find_longest_non_unique(seq, kmers_pos):
    max_length = 0
    extend_alternate = 0
    found_kmer_pos = []

    #loops through all the kmers in the dict
    for kmer in kmers_pos.iterkeys():
        kmer_len = len(kmer)

        #Step kmer through kmer and extend it as long as the pair is equal
        #Do that for all the current kmers and always update the length if
        # it is longer then the current max
        for i, pos in enumerate(kmers_pos[kmer][:-1]):
            for pos2 in kmers_pos[kmer][i+1:]:
                found_length = extend_kmers(seq, kmer_len, pos, pos2)

                if found_length > max_length:
                    found_kmer_pos = [pos, pos2]
                    max_length = found_length
                elif found_length == max_length:
                    found_kmer_pos.append(pos, pos2)

    return max_length, found_kmer_pos


def extend_kmers(seq, start_len, pos, pos2):
    not_unique_left = True
    not_unique_right = True

    prefix = 0
    suffix = 0

    while not_unique_left or not_unique_right:
        while not_unique_left:
            prefix += 1
            if seq[pos-prefix:pos+start_len+suffix] == seq[pos2-prefix:pos2+start_len+suffix]:
                not_unique_right = True
            else:
                not_unique_left = False
                prefix -= 1

        while not_unique_right:
            suffix += 1
            if seq[pos-prefix:pos+start_len+suffix] == seq[pos2-prefix:pos2+start_len+suffix]:
                not_unique_left = True
            else:
                not_unique_right = False
                suffix -= 1

    return prefix + start_len + suffix

if __name__ == '__main__':
    header, seq = read_fasta("MusChr01.fa.txt")
    kmers_pos = find_kmers_pos(seq, 14)
    max_length, found_kmers_pos = find_longest_non_unique(seq, kmers_pos)

    #t = timeit.Timer('header, seq = ad.read_fasta("MusChr01.fa.txt")',
    #                 setup='import analyse_dna as ad; ')
    #seconds = t.timeit(1)
    #print(seconds)
    #print(str(seconds/60))
    #unique = ad.find_unique_seq(seq)
    #unique = find_unique_seq(seq)
    #print(len(unique))
    #print(unique)
    #print(header)
    #print(len(seq))
    #print(find_unique_kmer_len(249 * 10 ** 6))
