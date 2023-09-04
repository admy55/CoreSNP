#!/usr/bin/env python3
# Chunchao Wang <admy55@gmail.com> & Tingyu Dou <tydou622@gmail.com>
# modified on 2023-07-21
import argparse
import gzip
import itertools
import logging
import os
import random
import re
import subprocess
from collections import namedtuple, Counter, deque
import numpy as np


def check_plink():
    try:
        subprocess.run(('plink', '--version'), stdout=subprocess.PIPE)
        plink_path = 'plink'
    except OSError:
        logger.error('Cannot find plink in your PATH.')
        if not os.path.exists('plink'):
            logger.error('Alternatively, put plink in the same directory with this script.')
            raise SystemExit

        os.chmod('plink', 0o755)
        plink_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'plink')
        try:
            subprocess.run((plink_path, '--version'), stdout=subprocess.PIPE)
            logger.info(f'Using "{plink_path}" instead.')
        except OSError:
            logger.error('Please make sure plink has executable permission or is complete.')
            raise SystemExit

    return plink_path


def get_args():
    parser = argparse.ArgumentParser(description='Generate core SNPs from a vcf file.')
    parser.add_argument(
        '-v',
        '--vcf',
        required=True,
        help='[filename] the input vcf file',
    )
    parser.add_argument(
        '-i',
        '--include',
        help='[filename] a file contained SNPs that must be included in the core set',
    )
    parser.add_argument(
        '-e',
        '--exclude',
        help='[filename] a file contained SNPs that would never be included in the core set',
    )
    parser.add_argument(
        '-x',
        '--flexing',
        type=int,
        default=1,
        choices=range(1, 6),
        help='the minimal number of candidate SNPs at each round (default: %(default)s)',
    )
    parser.add_argument(
        '-m',
        '--minimal',
        type=int,
        default=1,
        choices=range(1, 3),
        help='the minimal number of differential SNPs for each pairwise samples (default: %(default)s)',
    )
    parser.add_argument(
        '-c',
        '--count',
        type=int,
        default=1,
        choices=range(1, 11),
        help='the number of core sets this program generates (default: %(default)s)',
    )
    parser.add_argument(
        '-g',
        '--missing',
        type=float,
        default=0.2,
        help='the threshold of missing call frequency to filter variants (default: %(default)s)',
    )
    parser.add_argument(
        '-o',
        '--out',
        default='result',
        help='directory of output (default: %(default)s)',
    )
    parser.add_argument(
        '-l',
        '--log',
        default='coreSNP.log',
        help='filename of log file (default: %(default)s)',
    )
    parser.add_argument(
        '-M',
        '--more-info',
        action='store_true',
        help='print more info to log when running (default: %(default)s)',
    )

    return parser.parse_args()


def get_logger(file):
    logger = logging.getLogger('coreSNP')
    if args.more_info is True:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    formatter = logging.Formatter(
        '[%(asctime)s] %(name)s::%(funcName)s | %(levelname)s | %(message)s',
        datefmt='%Y/%m/%d %H:%M:%S',
    )

    fh = logging.FileHandler(file, 'w', encoding='utf-8')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)

    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    # ch.setFormatter(formatter)

    logger.addHandler(ch)
    logger.addHandler(fh)

    return logger


args = get_args()
if not os.path.exists(args.out):
    os.makedirs(args.out)

logger = get_logger(f'{args.out}/{args.log}')
plink = check_plink()


class FilePrefix:
    """
    Definition for prefix of filename.
    """

    def __init__(self, prefix):
        self.prefix = prefix

    @property
    def s1(self):
        return f'{self.prefix}.s1'

    @property
    def s2(self):
        return f'{self.prefix}.s2'

    @property
    def s3(self):
        return f'{self.prefix}.s3'

    @property
    def s4(self):
        return f'{self.prefix}.s4'

    @property
    def z(self):
        return self.prefix


RePattern = namedtuple(
    'pattern',
    (
        'ref',
        'alt',
        'missing_hetero',
        'missing_vcf',
        'missing_txt',
        'vcf_comment',
        'vcf_title',
        'space',
        'cap_double_w',
    ),
)
re_pattern = RePattern(
    ref=re.compile('0/0'),
    alt=re.compile('1/1'),
    missing_hetero=re.compile('0/1|\./\.'),
    missing_vcf=re.compile('\./\.'),
    missing_txt=re.compile('-'),
    vcf_comment=re.compile('^##'),
    vcf_title=re.compile('#CHROM'),
    space=re.compile('\s+'),
    cap_double_w=re.compile('(\w\w)'),
)

BI2SINGLE = {
    'AA': 'A',
    'GG': 'G',
    'CC': 'C',
    'TT': 'T',
    'AG': 'R',
    'GA': 'R',
    'CT': 'Y',
    'TC': 'Y',
    'AC': 'M',
    'CA': 'M',
    'GT': 'K',
    'TG': 'K',
    'GC': 'S',
    'CG': 'S',
    'AT': 'W',
    'TA': 'W',
    '00': '-',
}


def vcf2bed(file):
    """
    Preparing data for IBS, locus missing, individual missing, and MAF.

    Parameters:
        file: a .vcf or .vcf.gz file
    """

    script = (
        plink,
        '--vcf',
        file,
        '--geno',
        str(args.missing),
        '--double-id',
        '--allow-extra-chr',
        '--freq',
        '--missing',
        '--ibs-matrix',
        '--make-bed',
        '--recode',
        'vcf-iid',
        'bgz',
        '--out',
        f'{args.out}/data',
    )
    subprocess.run(script, stdout=subprocess.PIPE)


def get_major_allele(file):
    """
    Get major allele for SNP with a .bim file.

    Parameters:
        file: a .bim file

    Returns:
        major_allele: a dict, {SNP ID : major allele}
    """

    major_allele = {}
    with open(file, mode='r', encoding='utf-8') as f:
        for line in f:
            _, snp_id, *_, major = line.strip().split('\t')
            major_allele.setdefault(snp_id, major)

    return major_allele


def record_positions(pattern, text):
    """
    Return positions where re pattern located in a string.
    
    Parameters:
        pattern: a re.compile pattern
        text: a str string
    
    Returns:
        positions: a set, [int], zero-based
    """

    positions = set()
    matches = pattern.finditer(text)
    try:
        first_match = next(matches)
    except StopIteration:
        # logger.debug(f'Non missing exists.')
        pass
    else:
        positions = set(x.start() for x in itertools.chain((first_match, ), matches))

    return positions


def binary_matrix(text):
    """
    Replace REF to '0', ALT to '1', Missing & Hetero to '-'.

    Parameters:
        text: a str string

    Returns:
        text: a str string
    """

    text = re_pattern.ref.sub('0', text)
    text = re_pattern.alt.sub('1', text)
    text = re_pattern.missing_hetero.sub('-', text)

    return text


def fill_by_major(file):
    """
    Fill a genotype matrix with major allele.

    Parameters:
        file: a .vcf or .vcf.gz file
    
    Returns:
        sample_list: a list, [str]
        snp2genotype: a dict, {SNP ID : [genotype]}
    """

    major_snp = get_major_allele(f'{args.out}/data.bim')

    if file.endswith('gz'):
        fh = gzip.open(file, mode='rb')
    else:
        fh = open(file, mode='r', encoding='utf-8')

    snp2genotype = {}
    for line in fh:
        if file.endswith('gz'):
            line = line.decode()

        if re_pattern.vcf_comment.match(line):
            continue
        elif re_pattern.vcf_title.match(line):
            _, samples = line.split('FORMAT')
            sample_list = samples.strip().split('\t')
            continue

        _, _, snp_id, ref, *_, genotype = line.rstrip().split('\t', maxsplit=9)

        # Generate a matrix contains [0, 1, -].
        genotype = binary_matrix(genotype)

        # Replace missing with major allele.
        if ref == major_snp[snp_id]:
            allele = '0'
        else:
            allele = '1'
        genotype = re_pattern.missing_txt.sub(allele, genotype)

        snp2genotype[snp_id] = genotype.split('\t')

    return (sample_list, snp2genotype)


def read_snp_id(file):
    """
    Get SNP ID from a text file.

    Parameters:
        file: a text file containing SNP ID, one ID per line

    Returns:
        snp_id: a set, {SNP ID}
    """

    if not os.path.exists(file):
        raise FileNotFoundError(f'[{file}] does not exist, please check the filename or path.')

    snp_id = set()
    with open(file, mode='r', encoding='utf-8') as file_obj:
        for line in file_obj:
            if line.strip():
                snp_id.add(line.strip())

    return snp_id


def core_selection(prefix, samples, snp2geno):
    """
    Selecting a core SNP matrix.

    Parameters:
        prefix: prefix of filename
        samples: a list, [str]
        snp2geno: a dict, {SNP ID : genotype}
    """

    snp_file = f'{args.out}/{prefix.s1}.snp'
    haplo_file = f'{args.out}/{prefix.s1}.haplo'

    population_size = len(samples)
    snp_list = list(snp2geno.keys())
    removed_snp = set()

    # Check if user has provided a SNP list that must be included or excluded.
    included_snp = set()
    if args.include is not None:
        included_snp = read_snp_id(args.include)

    excluded_snp = set()
    if args.exclude is not None:
        excluded_snp = read_snp_id(args.exclude)
        removed_snp = excluded_snp.copy()

    both_in_include_and_exclude = included_snp.intersection(excluded_snp)
    if both_in_include_and_exclude:
        logger.warning(
            f'{both_in_include_and_exclude} is shared by [{args.include}] and [{args.exclude}], and it will be excluded.'
        )
        # included_snp = set(snp for snp in included_snp if snp not in both_in_include_and_exclude)
        for snp in both_in_include_and_exclude:
            included_snp.remove(snp)

    # Generate haplotype with included SNP.
    num = 0
    haplo_counter = 1
    primary_core = []
    core_haplo = []
    if included_snp:
        primary_core = list(included_snp)
        for snp in included_snp:
            num += 1
            logger.info(f'Including SNP #{num}: {snp}')
            removed_snp.add(snp)

            if not core_haplo:
                core_haplo = snp2geno[snp]
            else:
                core_haplo = [x + y for x, y in zip(core_haplo, snp2geno[snp])]

        haplo_counter = len(set(core_haplo))
        logger.debug(f'Haplotype count: {haplo_counter}')
        logger.debug(f'Haplotype frequency: {",".join(map(str, Counter(core_haplo).values()))}')

    # Select core SNP.
    while snp_list:
        shannon = {}
        combined_haplo = {}
        haplo_group = {}
        for snp in snp_list:
            if snp in removed_snp:
                continue

            if not core_haplo:
                haplotype = snp2geno[snp]
            else:
                haplotype = [x + y for x, y in zip(core_haplo, snp2geno[snp])]

            # Frequency of each haplotype.
            hap_freq = Counter(haplotype).values()
            if len(hap_freq) == haplo_counter:
                removed_snp.add(snp)
                continue

            combined_haplo[snp] = haplotype
            haplo_group[snp] = len(hap_freq)
            frequency = np.array(tuple(hap_freq)) / population_size
            shannon[snp] = -1 * (frequency * np.log(frequency)).sum()

        # Quit if no more SNP is added.
        if not shannon:
            break

        # Select a SNP.
        mx = max(Counter(shannon).values())
        candidate_snp = [k for k, v in shannon.items() if v == mx]

        # Update snp list, SNPs which could contribute a new haplotype were kept.
        snp_list = [*shannon]
        logger.debug(f'SNP number: {len(snp_list)}')
        if len(candidate_snp) < args.flexing:
            snp_list = sorted(shannon, key=shannon.get, reverse=True)
            candidate_snp.extend(snp_list[len(candidate_snp):args.flexing])
        logger.debug(f'Candidate SNPs: {candidate_snp}')

        # shannon_of_candidate_snp = [shannon[x] for x in candidate_snp]
        # logger.debug(f'Shannon of candidate SNPs: {shannon_of_candidate_snp}')

        num += 1
        selected_snp = random.choice(candidate_snp)
        logger.info(f'Selecting SNP #{num}: {selected_snp}')
        removed_snp.add(selected_snp)

        primary_core.append(selected_snp)
        core_haplo = combined_haplo[selected_snp]
        haplo_counter = haplo_group[selected_snp]
        logger.debug(f'Haplotype count: {haplo_counter}')
        logger.debug(f'Haplotype frequency: {",".join(map(str, Counter(core_haplo).values()))}')

    with open(snp_file, mode='w', encoding='utf-8') as out_snp:
        out_snp.write('\n'.join(primary_core))

    haplo = {}
    for i in range(population_size):
        haplo.setdefault(core_haplo[i], [])
        haplo[core_haplo[i]].append(samples[i])

    with open(haplo_file, mode='w', encoding='utf-8') as out_haplo:
        for hap, sample_list in haplo.items():
            out_haplo.write(f"{hap}: {','.join(sample_list)}\n")


def find_diff_index(string1, string2):
    """
    Given two strings, return the index list that the character is different in these positions.

    Parameters:
        string1: a str string
        string2: a str string

    Returns:
        index: a set, [int]
    """

    if len(string1) != len(string2):
        raise ValueError('These two strings should have the same length.')

    index = set()
    for i, (a, b) in enumerate(zip(string1, string2)):
        if a != b and a != '-' and b != '-':
            index.add(i)

    return index


def hamming_distance(string1, string2):
    """
    Calculate the Hamming distance between two strings.

    Parameters:
        string1: a str string
        string2: a str string

    Returns:
        distance, int.
    """

    if len(string1) != len(string2):
        raise ValueError('These two strings should have the same length.')

    return sum(a != b for a, b in zip(string1, string2))


def replace_missing(prefix):
    """
    Fill the missing SNP.

    Parameters:
        prefix: prefix of filename

    Returns:
        sample_pairs: a list, containing pairs of samples as tuple, [(sample1, sample2), ...]
        final_core: a list, final core SNP set
    """

    script = (
        plink,
        '--bfile',
        f'{args.out}/data',
        '--double-id',
        '--allow-extra-chr',
        '--extract',
        f'{args.out}/{prefix.s1}.snp',
        '--recode',
        '--out',
        f'{args.out}/{prefix.s2}',
    )
    subprocess.run(script, stdout=subprocess.PIPE)

    haplo2sample = {}
    haplo_with_missing = set()
    haplo_without_missing = set()
    with open(f'{args.out}/{prefix.s2}.ped', mode='r', encoding='utf-8') as f:
        for line in f:
            sample, *_, genotype = line.rstrip().split(maxsplit=6)
            genotype = re_pattern.space.sub('', genotype)
            haplo = re_pattern.cap_double_w.sub(lambda m: BI2SINGLE[m.group(1)], genotype)
            haplo2sample.setdefault(haplo, [])
            haplo2sample[haplo].append(sample)

            if re_pattern.missing_txt.search(haplo):
                haplo_with_missing.add(haplo)
            else:
                haplo_without_missing.add(haplo)

    # Recording samples that share the same haplotype.
    sample_pairs = deque()
    sample_pairs_add = deque()

    for miss_haplo, haplo in itertools.product(haplo_with_missing, haplo_without_missing):
        positions = record_positions(re_pattern.missing_txt, miss_haplo)
        miss_haplo_m = re_pattern.missing_txt.sub('', miss_haplo)
        haplo_m = ''.join(haplo[i] for i in range(len(haplo)) if i not in positions)
        if args.minimal == 1:
            if miss_haplo_m == haplo_m:
                sample1 = tuple(haplo2sample[miss_haplo])
                sample2 = tuple(haplo2sample[haplo])
                sample_pairs.append((sample1, sample2))
                logger.debug(f'{miss_haplo} ({",".join(sample1)})\t{haplo} ({",".join(sample2)})')
        elif args.minimal == 2:
            if miss_haplo_m == haplo_m:
                sample1 = tuple(haplo2sample[miss_haplo])
                sample2 = tuple(haplo2sample[haplo])
                sample_pairs.append((sample1, sample2))
                logger.debug(f'{miss_haplo} ({",".join(sample1)})\t{haplo} ({",".join(sample2)})')
            elif hamming_distance(miss_haplo_m, haplo_m) < 2:
                sample1 = tuple(haplo2sample[miss_haplo])
                sample2 = tuple(haplo2sample[haplo])
                sample_pairs_add.append((sample1, sample2))

    for haplo1, haplo2 in itertools.combinations(haplo_with_missing, 2):
        positions = record_positions(re_pattern.missing_txt, haplo1)
        positions.update(record_positions(re_pattern.missing_txt, haplo2))
        haplo1_m = ''.join(haplo1[i] for i in range(len(haplo1)) if i not in positions)
        haplo2_m = ''.join(haplo2[i] for i in range(len(haplo2)) if i not in positions)
        if args.minimal == 1:
            if haplo1_m == haplo2_m:
                sample1 = tuple(haplo2sample[haplo1])
                sample2 = tuple(haplo2sample[haplo2])
                sample_pairs.append((sample1, sample2))
                logger.debug(f'{haplo1} ({",".join(sample1)})\t{haplo2} ({",".join(sample2)})')
        elif args.minimal == 2:
            if haplo1_m == haplo2_m:
                sample1 = tuple(haplo2sample[haplo1])
                sample2 = tuple(haplo2sample[haplo2])
                sample_pairs.append((sample1, sample2))
                logger.debug(f'{haplo1} ({",".join(sample1)})\t{haplo2} ({",".join(sample2)})')
            elif hamming_distance(haplo1_m, haplo2_m) < 2:
                sample1 = tuple(haplo2sample[haplo1])
                sample2 = tuple(haplo2sample[haplo2])
                sample_pairs_add.append((sample1, sample2))

    if args.minimal == 2:
        for haplo1, haplo2 in itertools.combinations(haplo_without_missing, 2):
            if hamming_distance(haplo1, haplo2) < 2:
                sample1 = tuple(haplo2sample[haplo1])
                sample2 = tuple(haplo2sample[haplo2])
                sample_pairs_add.append((sample1, sample2))

    # There still exists samples that share the same haplotype, re-pick new SNPs to distinguish them.
    script = (
        plink,
        '--bfile',
        f'{args.out}/data',
        '--double-id',
        '--allow-extra-chr',
        '--exclude',
        f'{args.out}/{prefix.s1}.snp',
        '--recode',
        '--freq',
        '--out',
        f'{args.out}/{prefix.s3}',
    )
    subprocess.run(script, stdout=subprocess.PIPE)

    sample2geno = {}
    with open(f'{args.out}/{prefix.s3}.ped', mode='r', encoding='utf-8') as f:
        if args.minimal == 1:
            candidate_set = set(itertools.chain.from_iterable(itertools.chain.from_iterable(sample_pairs)))
        elif args.minimal == 2:
            candidate_set = set(
                itertools.chain.from_iterable(itertools.chain.from_iterable((*sample_pairs, *sample_pairs_add))))

        for line in f:
            sample, *_, genotype = line.strip().split(maxsplit=6)
            if sample in candidate_set:
                genotype = re_pattern.space.sub('', genotype)
                haplo = re_pattern.cap_double_w.sub(lambda m: BI2SINGLE[m.group(1)], genotype)
                sample2geno[sample] = haplo

    diff_sets_1 = deque()
    for sample1, sample2 in sample_pairs:
        diff_index = find_diff_index(sample2geno[sample1[0]], sample2geno[sample2[0]])
        if diff_index:
            diff_sets_1.append(diff_index)
        else:
            logger.warning(f'identical samples: ({",".join(sample1)}) -vs- ({",".join(sample2)})')

    if args.minimal == 2:
        diff_sets_2 = deque()
        for sample1, sample2 in sample_pairs_add:
            diff_index = find_diff_index(sample2geno[sample1[0]], sample2geno[sample2[0]])
            if diff_index:
                diff_sets_2.append(diff_index)
            else:
                logger.warning(f'only 1 SNP differ: ({",".join(sample1)}) -vs- ({",".join(sample2)})')

    # Read primary core SNPs.
    with open(f'{args.out}/{prefix.s1}.snp', mode='r', encoding='utf-8') as f:
        final_core = [i.strip() for i in f.readlines()]

    # If there is no more SNPs that could distinguish the sample pairs, return s1.snp as core SNP.
    if args.minimal == 1:
        try:
            next(itertools.chain.from_iterable(diff_sets_1))
        except StopIteration:
            return (sample_pairs, final_core)
    elif args.minimal == 2:
        try:
            next(itertools.chain.from_iterable((*diff_sets_1, *diff_sets_2)))
        except StopIteration:
            return (sample_pairs, final_core)

    # Fill missing by adding more SNPs.
    snp_index_to_id = {}
    snp_id_to_index = {}
    snp_maf = {}
    with open(f'{args.out}/{prefix.s3}.frq', mode='r', encoding='utf-8') as f:
        f.readline()  # .frq file has a title line, skip it
        for i, line in enumerate(f):
            _, snp_id, _, _, maf, _ = line.strip().split()
            snp_index_to_id[i] = snp_id
            snp_id_to_index[snp_id] = i
            snp_maf[snp_id] = float(maf)

    # Check if there is a perfect intersection set, if True then just return it.
    perfect_set = diff_sets_1[0].copy()
    for set_ in diff_sets_1:
        perfect_set.intersection_update(set_)

    num = len(final_core)
    selected = set()
    if perfect_set:
        # Fortunately, we got the perfect SNPs!
        logger.debug(f'Candidate SNP index: {perfect_set}')
        diff_snp = [snp_index_to_id[i] for i in perfect_set]
        selected_snp = max(diff_snp, key=snp_maf.get)
        final_core.append(selected_snp)
        num += 1
        logger.info(f'Selecting SNP #{num}: {selected_snp}')
    else:
        # Sadly, there is not a prefect intersection set, pick the SNPs which has the max frequency
        max_group = len(diff_sets_1)
        logger.debug(f'Number of sample pairs should be distinguished: {max_group}')

        temp_sets = diff_sets_1.copy()
        counter = 0
        while counter < max_group:
            frequency = Counter(itertools.chain.from_iterable(temp_sets))
            mx = max(Counter(frequency.values()))
            logger.debug(f'Top count of SNP in diff sets: {mx}')
            candidate = [k for k, v in frequency.items() if v == mx]
            logger.debug(f'Candidate SNP index: {candidate}')
            diff_snp = [snp_index_to_id[i] for i in candidate]
            selected_snp = max(diff_snp, key=snp_maf.get)
            final_core.append(selected_snp)
            num += 1
            logger.info(f'Selecting SNP #{num}: {selected_snp}')
            selected.add(snp_id_to_index[selected_snp])

            counter += mx
            logger.debug(f'Number of sample pairs have been distinguished: {counter} of {max_group}')
            if counter == max_group:
                break

            # Update temp_sets, remove index set that contains selected SNP
            temp_sets = tuple(s for s in temp_sets if snp_id_to_index[selected_snp] not in s)

    if args.minimal == 2:
        sample_pairs.extend(sample_pairs_add)
        diff_sets_1.extend(diff_sets_2)
        diff_sets = tuple(s - selected for s in diff_sets_1 if s - selected)

        max_group = len(diff_sets)
        logger.debug(f'Number of sample pairs should be distinguished: {max_group}')
        counter = 0
        while counter < max_group:
            frequency = Counter(itertools.chain.from_iterable(diff_sets))
            mx = max(Counter(frequency.values()))
            logger.debug(f'Top count of SNP in diff sets: {mx}')
            candidate = [k for k, v in frequency.items() if v == mx]
            logger.debug(f'Candidate SNP index: {candidate}')
            diff_snp = [snp_index_to_id[i] for i in candidate]
            selected_snp = max(diff_snp, key=snp_maf.get)
            final_core.append(selected_snp)
            num += 1
            logger.info(f'Selecting SNP #{num}: {selected_snp}')

            counter += mx
            logger.debug(f'Number of sample pairs have been distinguished: {counter} of {max_group}')
            if counter == max_group:
                break

            diff_sets = tuple(s for s in diff_sets if snp_id_to_index[selected_snp] not in s)

    return (sample_pairs, final_core)


def pic(freq):
    """
    Calculate PIC value using allele frequency.

    Parameters:
        freq: allele frequency

    Returns:
        SNP PIC value
    """

    # return 1 - freq**2 - (1 - freq)**2 - 2 * freq**2 * (1 - freq)**2
    return 2 * freq * (1 - freq) * (1 - freq * (1 - freq))


def check_result(prefix, sample_pairs):
    """
    Check if there still exists samples that has not been distinguished.

    Parameters:
        prefix: prefix of filename
        sample_pairs: pairs of samples that share the same haplotype before filling missing
    """

    sample2haplo = {}
    haplo2sample = {}
    with open(f'{args.out}/{prefix.s4}.ped', mode='r', encoding='utf-8') as f:
        for line in f:
            sample, *_, genotype = line.strip().split(maxsplit=6)
            # if sample in set(itertools.chain.from_iterable(sample_pairs)):
            genotype = re_pattern.space.sub('', genotype)
            haplo = re_pattern.cap_double_w.sub(lambda m: BI2SINGLE[m.group(1)], genotype)
            sample2haplo[sample] = haplo

            haplo2sample.setdefault(haplo, [])
            haplo2sample[haplo].append(sample)

    existed = set()
    for sample1, sample2 in sample_pairs:
        haplo1 = sample2haplo[sample1[0]]
        haplo2 = sample2haplo[sample2[0]]
        positions = record_positions(re_pattern.missing_txt, haplo1)
        positions.update(record_positions(re_pattern.missing_txt, haplo2))
        haplo1_m = ''.join(haplo1[i] for i in range(len(haplo1)) if i not in positions)
        haplo2_m = ''.join(haplo2[i] for i in range(len(haplo2)) if i not in positions)
        if haplo1_m == haplo2_m:
            existed.add(tuple(haplo2sample[haplo1]))
            existed.add(tuple(haplo2sample[haplo2]))
            logger.debug(f'{haplo1} {sample1}\t{haplo2} {sample2}')

    multiple = set(tuple(samples) for samples in haplo2sample.values() if len(samples) > 1)
    logger.debug(f'multiple samples share the same haplotype: {multiple}')
    logger.debug(f'cannot be distinguished: {existed}')
    if multiple - existed:
        logger.debug(f'{multiple - existed}')


def remove_files(prefix):
    for file in os.listdir(args.out):
        if re.match(f'{prefix.z}.s\d', file):
            os.remove(f'{args.out}/{file}')


def main(prefix, samples, snp2geno):
    logger.info(f'Selecting CoreSNP - {prefix.z}.')
    core_selection(prefix, samples, snp2geno)
    logger.info('Fix the missing introduced by major allele imputation.')
    sample_pairs, final_core = replace_missing(prefix)

    with open(f'{args.out}/{prefix.s4}.snp.id', mode='w', encoding='utf-8') as f:
        f.write('\n'.join(final_core))
    script = (
        plink,
        '--bfile',
        f'{args.out}/data',
        '--double-id',
        '--allow-extra-chr',
        '--extract',
        f'{args.out}/{prefix.s4}.snp.id',
        '--recode',
        '--freq',
        '--out',
        f'{args.out}/{prefix.s4}',
    )
    subprocess.run(script, stdout=subprocess.PIPE)

    haplo2sample = {}
    with open(f'{args.out}/{prefix.s4}.ped', mode='r', encoding='utf-8') as f:
        for line in f:
            sample, *_, genotype = line.strip().split(maxsplit=6)
            genotype = re_pattern.space.sub('', genotype)
            haplo = re_pattern.cap_double_w.sub(lambda m: BI2SINGLE[m.group(1)], genotype)
            haplo2sample.setdefault(haplo, [])
            haplo2sample[haplo].append(sample)

    with open(f'{args.out}/{prefix.z}.haplo', mode='w', encoding='utf-8') as f:
        for haplo, sample_list in haplo2sample.items():
            f.write(f"{haplo}: {','.join(sample_list)}\n")

    snp_pos = {}
    with open(f'{args.out}/{prefix.s4}.map', mode='r', encoding='utf-8') as f:
        for line in f:
            _, snp_id, _, position = line.strip().split()
            snp_pos[snp_id] = position

    out_snp = open(f'{args.out}/{prefix.z}.snp', mode='w', encoding='utf-8')
    out_snp.write('\t'.join(('SNP_ID', 'CHROM', 'POSITION', 'MAF', 'PIC\n')))
    with open(f'{args.out}/{prefix.s4}.frq', mode='r', encoding='utf-8') as f:
        f.readline()  # skip the title line
        for line in f:
            chrom, snp_id, _, _, maf, _ = line.split()
            pic_value = pic(float(maf))
            out_snp.write('\t'.join((snp_id, chrom, snp_pos[snp_id], maf, f'{pic_value}\n')))
    out_snp.close()

    # check_result(prefix, sample_pairs)
    logger.info(f'{prefix.z} done.\n')


if __name__ == '__main__':
    logger.info('Initializing.')
    vcf2bed(args.vcf)
    logger.info('Convert heterozygous SNP to missing, and impute by major allele.')
    (samples, genotype) = fill_by_major(f'{args.out}/data.vcf.gz')

    for i in range(args.count):
        prefix = FilePrefix(f'core_{i+1}')
        main(prefix, samples, genotype)
        remove_files(prefix)

# TODO
# 1. Because of missing filtration, if there still exists samples that share the same haplotype, use SNPs in the raw vcf.
# 2. In the 'core selection' step, remove the sample in next iteration if it is unique for one haplotype.
# 3. Adding a new module that selecting core SNP using 'diff sets' method.
# 4. Generate some graphs such as IBS/MAF/missing distribution, visualization of selecting, SNP on chromosome etc.
# 5. use multiprocessing on 'core selection'.
# 6. A GUI version for Windows users.

#
#
#                                                               d
#      ddddddddddddddddddddddddddddddd                           dddd
#      ddddddddddddddddddddddddddddddd                             d@
#                    ddd                                         *dddddddWWdd
#                    ddd                                  d dd@ddd   dd     ddd
#                    ddd                                ddd   ddd#dd @dddddddddd
#                    ddd                                         ddddddd*
#       ddddddddddddddddddddddddddddd                       dddddddd
#                    ddd                                     &ddddddd   ddd
#                    ddd                                    qdd@  ddddddddd
#                    ddd                                   d&     dd
#                    ddd                                   *dddbWddddddddddddd@
#                    ddd                             kddddddddddddd         B&dh
#    ddddddddddddddddddddddddddddddddddd                      ddd   8dddddp
#    ddddddddddddddddddddddddddddddddddd                   ddddd       ddddd
#                                                                         dd
#                    ddd
#      ddddddddddddddddddddddddddddddd                                 dd
#                   ddd                                               ddddh
#                  ddd                                       ddd              d
#        ddddddddddddddddddddddddddd                        ddd   qddddddddddddo
#                ddd                                        dd  @do   @dddddk
#     ddddddddddddddddddddddddddddddddd                    ddd  dd  ddd   ddd
#            dddd          #ddd                        pdddddddddk   ddddddd
#         &dddd               dddd                  dddd* dd  %dd   hdddddM   Bdd
#      dddddddddddddddddddddddddddddddM                  ddddddd dd     @ddd  ddddB
#    dddd   ddd               ddd   #ddd                    ddddd@  ddddddd
#           ddddddddddddddddddddd                       &ddd   Mdd     ddd
#           ddd               ddd                                      odd
#           ddd               ddd                                      ddd
#           ddddddddddddddddddddd                                   dddddd
#
#            ddd
#            ddd      kddddddddddddddddd                 ddd     ddddddddddd
#      dddddddddddddd#     @dd       ddd                 ddddh   dddddddd
#            ddd           ddd       ddd                            pddd
#            ddd          ddd        ddo                       dddddddd ddddd
#            ddd       dddd#  ddddddddb                ddddd      ddd  #ddd
#     dddddddddddddddd  d                          dddddddd      ddddddddddddddd
#             ddd       dddddddddddddddd                dd   ddddM
#        dd   ddd       ddd          ddd                dk   d  WdddddddddddW
#       bdd   dddddddd  ddd          ddd                dq*dd   dddd    #dddd
#       ddd   ddd       ddd          ddd                dddk     dddddddddd
#       dddd  ddd       dddddddddddddddd                           dddddddd
#      %dd dddddd       ddd          ddd
#     #dd&  addddd*
#     ddb       dddddddddddddddddddddddddM
#    h
