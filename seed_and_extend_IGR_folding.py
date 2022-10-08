__author__ = 'Alexander Frotscher'
__email__ = 'alexander.frotscher@student.uni-tuebingen.de'

import argparse
import math
import subprocess
from sys import stdin, stdout
import pandas as pd
import numpy as np
from subprocess import Popen, PIPE

import win32con


def read_fasta(file_name: str) -> list:
    """
    Reads the given FastA file and returns a list of pairs.

    Parameters
    ----------
    file_name : str
        Path or filename of the FastA file to be read.

    Returns
    -------
    List of pairs of headers and respective sequences.
    """

    sequences = []

    with open(file_name, 'r') as reader:
        current_seq = ['', '']

        for line in reader:
            if line.startswith('>'):
                if len(current_seq[1]) >= 1:
                    sequences.append(tuple(current_seq))
                current_seq[0] = line.rstrip()
                current_seq[1] = ''
            else:  # reading a line containing the sequence
                current_seq[1] += line.rstrip()
        sequences.append(tuple(current_seq))  # add the last tuple
    return sequences


def extend_window_right(depth: list, end_pos: int, expand_threshold: int, maximum_growth: int) -> int:
    """
    This function extends the core to the right according to the depth 

    Parameters
    ----------
    depth: list
        The read depth at each position of the IGR
    end_pos: int
        The ending position of the core
    expand_threshold: int
        The threshold that determines whether the sequence is increased
    maximum_growth: int
        The maximum number of nucleotides the core is allowed to grow
    Returns
    -------
    The final ending position of the possible ncRNA
    """
    if end_pos == len(depth):  # the last base of the sequence can not extend the IGR
        return end_pos
    else:
        final_pos = end_pos
        for i in range(end_pos + 1, len(depth)):  # increase the final ending position if possible
            if final_pos - end_pos > math.floor(
                    0.5 * maximum_growth):  # ncRNA is not larger than defined, two segments -> left, right -> 0.5
                break
            elif depth[i] >= expand_threshold:
                final_pos += 1
            else:
                break
        return final_pos


def extend_window_left(depth: list, start_pos: int, expand_threshold: int, maximum_growth: int) -> int:
    """
    This function extends the core to the left according to the depth 

    Parameters
    ----------
    depth: list
        The read depth at each position of the IGR
    start_pos: int
        The starting position of the core
    expand_threshold: int
        The threshold that determines whether the sequence is increased
    maximum_growth: int
        The maximum number of nucleotides the core is allowed to grow
    Returns
    -------
    The final starting position of the possible ncRNA
    """
    if start_pos == 0:  # the first base of the sequence can not extend the IGR
        return start_pos
    else:
        final_pos = start_pos
        for i in range(start_pos - 1, 0, -1):  # decrease the final starting position if possible
            if start_pos - final_pos > math.floor(
                    0.5 * maximum_growth):  # ncRNA is not larger than defined, two segments -> left, right -> 0.5
                break
            if depth[i] >= expand_threshold:
                final_pos -= 1
            else:
                # maybe calc energy
                break
        return final_pos


def find_core(igr: str, depth: list, window_size: int, depth_cutoff: int, expand_threshold: int,
              maximum_growth: int, mode: str) -> list:
    """
    This function finds the core regions in the IGR 

    Parameters
    ----------
    igr: str
        The IGR to be searched
    depth: list
        The read depth at each position of the IGR
    window_size: int
        The size of the window for scanning the IGR
    depth_cutoff: int
        The threshold to determine a new ncRNA
    expand_threshold: int
        The threshold that determines whether the sequence is increased
    maximum_growth: int
        The maximum number of nucleotides the core is allowed to grow
    mode: str
        The mode for protein or RNA
    Returns
    -------
    The positions of all cores that could correspond to an ncRNA
    """
    my_cores = []
    i = 0  # start position
    if sum(depth) > depth_cutoff:
        while i < len(igr):
            j = i + window_size
            score = sum(depth[i:j])
            if score >= depth_cutoff:
                start_pos = extend_window_left(depth, i, expand_threshold, maximum_growth)
                end_pos = extend_window_right(depth, j, expand_threshold, maximum_growth)
                if mode == 'r':
                    # p = subprocess.Popen('C:\Program Files (x86)\ViennaRNA Package\RNAfold.exe', stdin=PIPE,
                    #                     stdout=PIPE)
                    # ans = p.communicate(igr[start_pos:end_pos].encode())
                    with open('fold.txt', 'w+') as writer:
                        writer.write(f'>my RNA\n{igr[start_pos:end_pos]}')
                    ans = subprocess.check_output(
                        f'C:\Program Files (x86)\ViennaRNA Package\RNAfold.exe fold.txt')
                    # ans = ans.decode().split(' ')[-1].strip()[:-1]
                    ans = ans.decode().split(' ')[-1].strip().replace('(', '').replace(')', '')
                    if float(ans) < -20.0:
                        start_pos = extend_window_left(depth, i, expand_threshold/2, maximum_growth)
                        end_pos = extend_window_right(depth, j, expand_threshold/2, maximum_growth)
                my_cores.append([start_pos, end_pos])
                i = end_pos
            else:
                i += 1
    return my_cores


def main():
    # check for mode
    if args.category == 'r':
        if args.window is None:
            window_size = 21
        if args.growth is None:
            maximum_growth = 550
    else:
        if args.window is None:
            window_size = 100
        if args.growth is None:
            maximum_growth = 1500

    # setup the depth as a list
    my_cores = []
    df1 = pd.read_csv(args.depth, engine='c')
    genes = df1['gene_name'].unique()
    sequences = read_fasta(args.sequence)
    for gen, seq in zip(genes, sequences):
        if gen not in seq[0]:
            print(
                f'The gen {gen} is not in the fasta file or they do not have the same order. It is excluded from the analysis')
        else:
            if len(seq[1]) >= window_size:
                bases = df1[df1['gene_name'] == gen]
                bases = bases.iloc[:, 2:]
                for sample_and_time in bases:
                    depth = bases[sample_and_time].to_numpy()
                    if len(seq[1]) != len(depth):
                        print(
                            f'Sequence length is not the same as coverage length! {gen} is not taken into account for {sample_and_time}. ')
                        if len(seq[1]) - 1 == len(depth):
                            print(f'The problem is samtools. Let me help you with that!'
                                  f' Appends depth count that is the mean over the depths for first base!')
                            np.insert(depth, 0, np.mean(depth))
                            cores = find_core(seq[1], depth, window_size, args.threshold, args.expand, maximum_growth,
                                              args.category)
                    else:
                        cores = find_core(seq[1], depth, window_size, args.threshold, args.expand, maximum_growth,
                                          args.category)
                    if len(cores) > 0:
                        for i in range(len(cores) - 1, 0, -1):
                            core = cores[i]
                            if core[0] <= cores[i - 1][1]:
                                cores[i - 1][1] = core[1]
                                cores.pop(i)
                        for core in cores:
                            core.append(sum(depth[core[0]:core[1]]))
                        my_cores.append((sample_and_time, gen, seq[1], cores))
    with open(args.output, 'w') as w:
        for sample, gen, sequence, cores in my_cores:
            for core in cores:
                w.write(f'>{gen} in {sample} from {core[0]} to {core[1]} has the score {core[2]}\n')
                for i in range(core[0], core[1], +80):
                    if i + 80 < len(sequence):
                        w.write(sequence[i:i + 80] + '\n')
                    else:
                        w.write(sequence[i:len(sequence)] + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Determine sequences in IGRs over multiple timestamps or samples with seed-and-extend.')
    parser.add_argument('-s', '--sequence', required=True, metavar='*.fasta',
                        help='The IGRs to be read')
    parser.add_argument('-d', '--depth', required=True, metavar='*.csv',
                        help='The read depth of the corresponding IGRs for every base. Read in as a matrix where the '
                             'columns correspond to the samples with IGR annotation and the rows are the depths.')
    parser.add_argument('-c', '--category', required=True, choices=['p', 'r'], metavar='',
                        help='The sequences the program is meant to find. p for protein, r for RNA')
    parser.add_argument('-t', '--threshold', required=True, metavar='', type=int,
                        help='The depth threshold over the complete window in order to determine a gene/ncRNA in the '
                             'IGR. Should be the median of the sum off the depth counts per base for the '
                             'comparable genes, e.g., tRNA depth counts for detecting ncRNA and CDS counts for proteins')
    parser.add_argument('-w', '--window', metavar='', type=int,
                        help='The window size for scanning the IGR. Default value for ncRNA is 21 and for protein 100')
    parser.add_argument('-e', '--expand', default=30, metavar='', type=int,
                        help='The depth that has to be achieved in order to extend the window and gain larger '
                             'sequences. Can be the median of the depth counts for the comparable genes, e.g., '
                             'tRNA depth counts for detecting ncRNA and CDS counts for proteins. Default value is 30 ')
    parser.add_argument('-g', '--growth', metavar='', type=int,
                        help='The number of nucleotides the window is allowed to grow. The larger the value the '
                             'larger the found sequences are. Default value for ncRNA is 550 and for proteins 1500')
    parser.add_argument('-o', '--output', required=True, metavar='',
                        help='The output file containing the found sequences')
    args = parser.parse_args()
    main()
