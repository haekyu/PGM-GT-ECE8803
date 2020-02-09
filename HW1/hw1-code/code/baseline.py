'''''''''''''''''''''''''''''''''''''''''''''''''''''''''
Probabilistic Graphical Models

Homework 1 - 4.1 Baseline

Haekyu Park
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''

'''
Import packages
'''
import os
from collections import defaultdict



'''
Functions for 4.1.1
'''


def parse_word_counts(training_path, new_training_path):
    '''
    Parse the word counts by replacing infrequent words
    in the original training data file with _RARE_.
    '''

    # Read the training dataset and count the words
    word_cnt = defaultdict(lambda: 0)
    with open(training_path) as f:
        while True:
            line = f.readline()
            if line == '':
                break
            if line == '\n':
                continue
            word = line.split()[0]
            word_cnt[word] += 1

    # Generate a new training dataset with _RARE_ symbol.
    with open(training_path) as f_read:
        with open(new_training_path, 'w') as f_write:
            while True:

                # Read line
                line = f_read.readline()
                if line == '':
                    break

                if line == '\n':
                    f_write.write(line)
                    continue

                # If infrequent word, replace the word into _RARE_
                word = line.split()[0]
                if word_cnt[word] < 5:
                    tag = line.split()[1]
                    line = '{} {}\n'.format('_RARE_', tag)

                f_write.write(line)


def run_count_freqs():
    '''
    Run count_freqs.py, the given code
    to count out the frequency of words
    '''
    os.system('sh run_gene_count.sh')



'''
Functions for 4.1.2
'''


def count_word_tag(count_path):
    '''
    Count the tags and word-tag pairs
    '''

    word_tag_cnt = defaultdict(lambda: 0)
    tag_cnt = defaultdict(lambda: 0)
    word_tag_dict = defaultdict(lambda: {})

    with open(count_path) as f:
        while True:
            line = f.readline()
            if line == '':
                break

            tokens = line.split()
            if 'WORDTAG' in tokens:
                word = tokens[-1]
                tag = tokens[-2]
                cnt = int(tokens[0])
                word_tag_cnt[(word, tag)] += cnt
                word_tag_dict[word][tag] = True
            elif '1-GRAM' in tokens:
                tag = tokens[-1]
                cnt = int(tokens[0])
                tag_cnt[tag] += cnt

    for word in word_tag_dict:
        word_tag_dict[word] = list(word_tag_dict[word].keys())

    return word_tag_cnt, tag_cnt, word_tag_dict


def calculate_emission(word_tag_cnt, tag_cnt):
    '''
    Compute emission parameters
    '''

    emission = defaultdict(lambda: 0)
    for word, tag in word_tag_cnt:
        emission[(word, tag)] = word_tag_cnt[(word, tag)] / tag_cnt[tag]

    return emission



'''
Functions for 4.1.3
'''


def simple_gene_tagger(count_path, test_path, test_tag_path, word_tag_dict, emission):
    '''
    Tag for each word with the new counts
    including _RARE_ words.
    '''

    with open(test_path) as f_read:
        with open(test_tag_path, 'w') as f_write:
            while True:
                # Read the word in the test set
                word = f_read.readline()
                if word == '':
                    break
                if word == '\n':
                    f_write.write(word)
                    continue
                word = word.replace('\n', '')

                # _RARE_
                write_word = word
                if word not in word_tag_dict:
                    word = '_RARE_'

                # Find out the best tag
                max_ems = -1000
                max_tag = ''
                for cand_tag in word_tag_dict[word]:
                    ems = emission[(word, cand_tag)]
                    if max_ems < ems:
                        max_ems = ems
                        max_tag = cand_tag
                f_write.write('{} {}\n'.format(write_word, max_tag))


def save_test_tag(test_tag_path, test_word_tag, test_word_sequence):
    '''
    Save the test tag
    '''
    with open(test_tag_path, 'w') as f:
        for word in test_word_sequence:
            if word == '\n':
                f.write(word)
            else:
                tag = test_word_tag[word]
                f.write('{} {}\n'.format(word, tag))


def evaluate_baseline():
    '''
    Run eval_gene_tagger.py, the given code
    to evaluate the baseline tagger.
    '''
    os.system('sh run_eval_baseline.sh')



if __name__ == '__main__':

    # 4.1.1
    training_path = '../hmm/gene.train'
    new_training_path = '../hmm/gene.new.train'

    parse_word_counts(training_path, new_training_path)
    run_count_freqs()

    # 4.1.2
    count_path = '../hmm/gene.new.counts'
    word_tag_cnt, tag_cnt, word_tag_dict = count_word_tag(count_path)
    emission = calculate_emission(word_tag_cnt, tag_cnt)

    # 4.1.3
    test_path = '../hmm/gene.test'
    test_tag_path = '../hmm/gene_test.p1.out'
    simple_gene_tagger(count_path, test_path, test_tag_path, word_tag_dict, emission)
    evaluate_baseline()
