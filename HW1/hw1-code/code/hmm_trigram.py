'''''''''''''''''''''''''''''''''''''''''''''''''''''''''
Probabilistic Graphical Models

Homework 1 - 4.2 HMM with trigram features

Haekyu Park
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''

'''
Import packages
'''
import os
from collections import defaultdict


def count_word_grams(count_path):
    '''
    Count word, tag, word-tag, 2-GRAM, 3-GRAM
    '''

    word_cnt = defaultdict(lambda: 0)
    tag_cnt = defaultdict(lambda: 0)
    word_tag_cnt = defaultdict(lambda: 0)
    bigram_cnt = defaultdict(lambda: 0)
    trigram_cnt = defaultdict(lambda: 0)

    with open(count_path) as f:
        while True:

            # Read line
            line = f.readline()
            if line == '':
                break

            # Count word
            if 'WORDTAG' in line:
                # ex) 92 WORDTAG O reading
                tokens = line.split()
                cnt = int(tokens[0])
                tag = tokens[-2]
                word = tokens[-1]
                word_cnt[word] += cnt
                word_tag_cnt[(word, tag)] += cnt
            #
            # Count tag
            elif '1-GRAM' in line:
                # ex) 41072 1-GRAM I-GENE
                tokens = line.split()
                cnt = int(tokens[0])
                tag = tokens[-1]
                tag_cnt[tag] += cnt
            #
            # Count 2-GRAM
            elif '2-GRAM' in line:
                # ex) 24435 2-GRAM I-GENE I-GENE
                tokens = line.split()
                cnt = int(tokens[0])
                tag1 = tokens[-2] # y_{i-1}
                tag2 = tokens[-1] # y_{i}
                bigram_cnt[(tag1, tag2)] += cnt
            #
            # Count 3-GRAM
            elif '3-GRAM' in line:
                # ex) 11 3-GRAM I-GENE I-GENE STOP
                tokens = line.split()
                cnt = int(tokens[0])
                tag1 = tokens[-3] # y_{i-2}
                tag2 = tokens[-2] # y_{i-1}
                tag3 = tokens[-1] # y_{i}
                trigram_cnt[(tag1, tag2, tag3)] += cnt

    return word_cnt, tag_cnt, word_tag_cnt, bigram_cnt, trigram_cnt


def q_global(y_i, y_i_2, y_i_1, bigram_cnt, trigram_cnt):
    '''
    Compute q parameter
    '''
    return trigram_cnt[(y_i_2, y_i_1, y_i)] / bigram_cnt[(y_i_2, y_i_1)]


def e_global(x, y, word_cnt, word_tag_cnt, tag_cnt):
    '''
    Compute emission where x is a word and y is a tag
    '''
    if x in word_cnt:
        return word_tag_cnt[(x, y)] / tag_cnt[y]
    else:
        return word_tag_cnt[('_RARE_', y)] / tag_cnt[y]


def viterbi(sentence, word_cnt, tag_cnt, bigram_cnt, trigram_cnt):
    '''
    Run Viterbi algorithm to solve the MLE problem
    - sentence = x_1, x_2, ..., x_n
    '''

    # Define functions to compute q and e
    def q(y_i, y_i_2, y_i_1):
        return q_global(y_i, y_i_2, y_i_1, bigram_cnt, trigram_cnt)

    def e(x, y):
        return e_global(x, y, word_cnt, word_tag_cnt, tag_cnt)


    # Get the number of words
    n = len(sentence)

    # Initialize tags
    y = {}

    # Initialize K
    K = {}
    K[-1] = ['*']
    K[0] = ['*']
    possible_tags = list(tag_cnt.keys())
    for k in range(1, n + 1):
        K[k] = possible_tags[:]

    # Initialize pi (the product of q values)
    pi = {(0, '*', '*'): 1}

    # Initialize bp
    bp = {}

    # Run Viterbi algorithm

    # 1. Compute pi
    for k in range(1, n + 1):
        x_k = sentence[k - 1]

        for u in K[k - 1]:
            for v in K[k]:
                # Find out the maximum pi and
                # w which makes the maximum pi
                max_pi = -1000
                max_pi_w = None

                for w in K[k - 2]:

                    curr_pi = pi[(k - 1, w, u)] * q(v, w, u) * e(x_k, v)
                    if curr_pi > max_pi:
                        max_pi = curr_pi
                        max_pi_w = w

                pi[(k, u, v)] = max_pi
                bp[(k, u, v)] = max_pi_w

    # 2. Get the tail tags, i.e., y_{n-1} and y_n
    argmax_u = None
    argmax_v = None
    max_pi = -1000

    for u in K[n - 1]:
        for v in K[n]:
            curr_pi = pi[(n, u, v)] * q('STOP', u, v)
            if curr_pi > max_pi:
                max_pi = curr_pi
                argmax_u = u
                argmax_v = v

    y[n - 1] = argmax_u
    y[n] = argmax_v

    # 3. Get the remaining tags
    for k in range(n - 2, 0, -1):
        y[k] = bp[(k + 2, y[k + 1], y[k + 2])]

    # Get the tags in order
    ys = []
    for k in range(n):
        ys.append(y[k + 1])

    return ys


def test_tag(test_path, out_path, word_cnt, tag_cnt, bigram_cnt, trigram_cnt):
    '''
    Read the test set and save the tag
    '''

    with open(test_path) as f_read:
        with open(out_path, 'w') as f_write:

            sentence = []

            while True:
                line = f_read.readline()
                if line == '':
                    break

                if line == '\n':
                    ys = viterbi(sentence, word_cnt, tag_cnt, bigram_cnt, trigram_cnt)
                    for x, y in zip(sentence, ys):
                        f_write.write('{} {}\n'.format(x, y))
                    sentence = []
                    f_write.write('\n')
                else:
                    word = line.replace('\n', '')
                    sentence.append(word)


def evaluate_viterbi():
    '''
    Run eval_gene_tagger.py, the given code
    to evaluate the viterbi tagger.
    '''
    os.system('sh run_eval_viterbi.sh')


if __name__ == '__main__':
    # 4.2.1
    count_path = '../hmm/gene.new.counts'
    word_cnt, tag_cnt, word_tag_cnt, bigram_cnt, trigram_cnt = count_word_grams(count_path)

    # 4.2.2
    test_path = '../hmm/gene.test'
    out_path = '../hmm/gene_test.p2.out'
    test_tag(test_path, out_path, word_cnt, tag_cnt, bigram_cnt, trigram_cnt)
    evaluate_viterbi()
