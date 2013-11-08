__author__ = 'markus'

import sys
import os
import math
from itertools import zip_longest
from multiprocessing import Process, Queue

#seq = ''


def read_fasta(filename):
    """ Reads a sequence in Fasta format """
    seq = ''
    with open(filename, 'r') as fp:
        header = ""
        for line in fp:
            if line == "":
                break
            if line.startswith('>'):
                header = line[1:].strip()
            else:
                seq += line.strip().upper()
    return seq


def calc_unique_kmer_len(seq_len):
    residual = seq_len / 4
    kmer_len = 1
    while residual >= 1.0:
        residual /= 16
        kmer_len += 2
    if residual * 4 <= 1.0:
        kmer_len -= 1
    return kmer_len


#decodes the numbers back into the original kmer
def to_kmer(num):
    translation = {5: 'A', 7: 'C', 4: 'T', 1: 'G'}
    kmer = ''
    while num:
        digit = num % 10
        num //= 10
        kmer += translation.get(digit)

    return kmer[::-1]


#encodes a kmer to numbers, to save storage
def to_digits(kmer):
    return int(''.join(map(str, [ord(char) % 10 for char in kmer.upper()])))


def find_kmers_pos(kmer_len, unique):
    """Finds all the positions of the given kmer length
    and saves them into a dictionary

    Keyword arguments:
    seq -- The sequence that contains all the kmers
    kmer_len -- length of the kmers we want to find
    unique -- If True only returns unique kmers
              If False returns all kmers occuring > 2.
    """
    global seq
    seq_len = len(seq)
    kmers_pos = {}

    for i in range(seq_len - kmer_len):
        curr_kmer = seq[i:i+kmer_len]
        if 'N' in curr_kmer:
            continue

        curr_kmer = to_digits(curr_kmer)
        kmers_pos[curr_kmer] = kmers_pos.get(curr_kmer, [])
        kmers_pos[curr_kmer].append(i)

    print("Filtering Dictionary...")
    if unique:
        kmers_pos = {k: v for k, v in kmers_pos.items() if len(v) == 1}
    else:
        kmers_pos = {k: v for k, v in kmers_pos.items() if len(v) > 2}
    return kmers_pos


def find_longest_non_unique(kmers_pos):
    """Finds the longest non unique kmer
    Loops through all the kmers and withing one kmer position-list,
    makes all n choose 2 (n=length of the list) kmer 'comparisons'
    or 'extensions'. Where-as a comparison Extends both kmers on each
    side as long as they match. If the found length is longer than the
    before found length, max_length gets updated.

    Keyword arguments:
    seq -- The sequence that contains all the kmers
    kmers_pos -- A dictionary with kmers plus its starting positions (k: [pos1, pos2, etc])
    """
    max_length = 0
    found_kmer_pos = []
    c = 0
    nr_different_kmers = len(kmers_pos)
    global seq
    print('Process {0} started. Shared sequence memory reference: {2}'.format(os.getpid(), len(seq), id(seq)))
    #loops through all the kmers in the dict
    for kmer in list(kmers_pos.keys()):
        c += 1

        kmer_len = len(kmer)

        #Step kmer through kmer and extend it as long as the pair is equal
        #Do that for all the current kmers and always update the length if
        # it is longer then the current max
        for i, pos in enumerate(kmers_pos[kmer][:-1]):
            for pos2 in kmers_pos[kmer][i+1:]:
                found_length = extend_kmers(kmer_len, pos, pos2)

                if found_length > max_length:

                    print("Length updated to: " + str(found_length))
                    print("k-mer " + str(c) + " out of " + str(nr_different_kmers) + " in the dictionary.\n")
                    found_kmer_pos = [pos, pos2]
                    max_length = found_length
                elif found_length == max_length:
                    print("Length " + str(found_length) + " found again.")
                    print("k-mer " + str(c) + " out of " + str(nr_different_kmers) + " in the dictionary.")
                    found_kmer_pos.extend((pos, pos2))
        del kmers_pos[kmer]
    return max_length, found_kmer_pos


def extend_kmers(start_len, pos1, pos2):
    """Extend given kmers, as long as they are unique
    It extends both of the kmers on the left sight, then on the right side and whenever
    there was an extension on either side on the other side again. Only if there was no
    extension on both sides it will return the length of the current, extended kmer

    Keyword arguments:
    seq -- The sequence that contains all the kmers
    start_len -- the length of both kmers before any extending
    pos, pos2 -- the starting positions of the kmers to extend and compare
    """
    global seq
    not_unique_left = True
    not_unique_right = True

    prefix = 0
    suffix = 0
    seq_len = len(seq)

    while not_unique_left or not_unique_right:
        while not_unique_left:
            prefix += 1
            pos1_len = pos1 - prefix
            pos2_len = pos2 - prefix
            if (pos1_len >= 0 and pos2_len >= 0 and
                        seq[pos1_len] != 'N' and
                        seq[pos1_len] == seq[pos2_len]):
                not_unique_right = True
            else:
                not_unique_left = False
                prefix -= 1

        while not_unique_right:
            suffix += 1
            pos1_len = pos1 + start_len + suffix
            pos2_len = pos2 + start_len + suffix
            if (pos1_len < seq_len and pos2_len < seq_len and
                        seq[pos1_len] != 'N' and
                        seq[pos1_len] == seq[pos2_len]):
                not_unique_left = True
            else:
                not_unique_right = False
                suffix -= 1

    return prefix + start_len + suffix


if __name__ == '__main__':
    seq = ''
    print("Reading file...")
    seq = read_fasta("MusChr01.fa.txt")
    print('Done.')
    #print(seq)

    if 1 == len(sys.argv) or sys.argv[1] == "unique":
        if len(sys.argv) == 3:
            kmers_len = int(sys.argv[2])
        else:
            kmers_len = 9

        print("Finding positions of all " + str(kmers_len) + "-mers...")
        kmers_pos = find_kmers_pos(kmers_len, True)
        print('Positions of the unique 9-mers:')
        print(kmers_pos)
        print('Kmers:')
        for pos in kmers_pos.values():
            print(str(seq[pos[0]:pos[0]+kmers_len]))

    else:
        if len(sys.argv) == 3:
            kmers_len = int(sys.argv[2])
        else:
            kmers_len = 14

        nr_proc = 4
        if len(sys.argv) == 4:
            nr_proc = int(sys.argv[3])

        print("Finding positions of all " + str(kmers_len) + "-mers...")
        kmers_pos = find_kmers_pos(kmers_len, False)
        print("Finding longest non-unique k-mer with " + str(nr_proc) + " Processes...")
        kmers_pos = zip_longest(*[iter(kmers_pos.items())]*math.ceil(len(kmers_pos)/nr_proc))
        result_queue = Queue()
        processes = []
        for chunk in zip(kmers_pos):
            chunk = dict(item for item in chunk[0] if item is not None)
            p = Process(target=find_longest_non_unique, args=(chunk, result_queue))
            processes.append(p)
            p.start()

        for p in processes:
            p.join()

        print('All jobs finished. Selecting best results...')

        max_length = 0
        found_kmers_pos = []
        for result in [result_queue.get() for p in processes]:
            print(result[0])
            if result[0] > max_length:
                max_length = result[0]
                found_kmers_pos = result[1]
            elif result == max_length:
                found_kmers_pos.append(result[1])

        print("Length of longest non-unique k-mer:")
        print(max_length)
        print("Starting positions of said k-mers:")
        print(list(set(found_kmers_pos)))
        print("K-mer: " + seq[found_kmers_pos[0]:found_kmers_pos[0]+max_length])
