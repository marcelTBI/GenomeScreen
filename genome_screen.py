import numpy as np
import pandas as pd
import argparse
import os
from datetime import datetime

# import R
import rpy2.robjects
import rpy2.robjects.numpy2ri
import rpy2.robjects.pandas2ri


all_chromosomes = list(range(24))
all_autosomes = list(range(22))

# create chromosome list:
chromosome_list = ['chr%d' % i for i in range(1, 23)] + ['chrX', 'chrY']
chrom_nums = {k: v for v, k in enumerate(chromosome_list)}
to_int = np.vectorize(lambda x: int(chrom_nums[x]))


def load_arguments():
    """
    Loads and parses arguments.
    :return: argparse arguments
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    required = parser.add_argument_group('Required')
    required.add_argument('-m', '--mean', type=nonempty_file, required=True, help="Mean count of bins. Numpy file - either .npy or .npz with 'means' entry.")
    required.add_argument('-b', '--bins', type=nonempty_file, required=True, help="Bin count of a sample. Numpy file - either .npy or .npz with 'cnv' entry.")

    defaults = parser.add_argument_group('Default')
    defaults.add_argument('-w', '--window', type=positive_nonzero_int, default=20000, help="Window size for a bin. Default:20000")
    defaults.add_argument('-r', '--rscript', type=nonempty_file, default=None, help="Path to the R-script.")

    output = parser.add_argument_group('Output')
    output.add_argument('-o', '--output', type=str, default=None, help="Directory path for output files.")
    output.add_argument('-n', '--name', type=str, default=None, help="Name prefix of the output file. Use --bins if not provided.")

    misc = parser.add_argument_group('Miscellaneous')
    misc.add_argument('-x', '--x-count', type=int, default=None, help="Number of X chromosomes 1 for a male sample, 2 for a female sample. Used only for gonosomal chromosomes.")
    misc.add_argument('-l', '--min-length', type=int, default=5000000,
                      help="Minimal effective length of detection (in bases) to report. Supply -1 for infinity, 0 for no restriction. Default=5,000,000.")

    filtration = parser.add_argument_group('Filtration')
    filtration.add_argument('-k', '--filter-on', action='store_true', help="Do filtering of bins to predict CNV on.")
    filtration.add_argument('--minimal-mean', type=float, default=3.0,
                            help="Minimal mean number of reads per bin. Average is around 7 (7.7 for a healthy bin), use numbers in range (0.0 - 6.0). Default=3.0")
    filtration.add_argument('--maximal-mean', type=float, default=None,
                            help="Maximal mean number of reads per bin. Average is around 7 (7.7 for a healthy bin), use numbers greater than 15.0. Default=None(infinity)")
    filtration.add_argument('--maximal-variance', type=float, default=1.5,
                            help="Maximal variance allowed in a bin. Usual variance of a healthy bin is around '0.85', use numbers greater than '1.1'. Default=1.5")
    filtration.add_argument('--filter-file', type=nonempty_file, default=None,
                            help="File, from where load bad bins (gaps). This option overrides other filtration arguments.")

    args = parser.parse_args()

    # assign defaults:
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if args.rscript is None:
        args.rscript = '%s/cbs_lib.R' % script_dir
    if args.output is None:
        args.output = '.'
    if args.maximal_mean is None:
        args.maximal_mean = np.inf
    if args.name is None:
        name = args.bins.split('/')[-1]
        args.name = name[:name.rfind('.')]
    if args.x_count == 0:
        args.x_count = None
    if args.min_length < 0:
        args.min_length = np.inf
    if args.min_length == 0:
        args.min_length = None

    return args


def positive_float(value, max_limit=None):
    """
    Represents positive float
    :param value: string value to estimate
    :param max_limit: maximal allowed value, skip validation if None
    :return: integer value, if can be converted into positive int else ArgumentTypeError
    """
    try:
        float_value = float(value)
    except ValueError:
        error = "Value %s is not float" % value
        # report.log_str(error, priority=logging.ERROR)
        raise argparse.ArgumentTypeError(error)
    if float_value <= 0:
        error = "Value %s is not positive float" % value
        # report.log_str(error, priority=logging.ERROR)
        raise argparse.ArgumentTypeError(error)
    if max_limit and max_limit < float_value:
        error = "Value %s must be lower or equal to %s" % (value, max_limit)
        # report.log_str(error, priority=logging.ERROR)
        raise argparse.ArgumentTypeError(error)
    return float_value


def positive_int(value, max_limit=None):
    """
    Represents positive decimal number, 0 included
    :param value: string value to estimate
    :param max_limit: maximal allowed value, skip validation if None
    :return: integer value, if can be converted into positive int else ArgumentTypeError
    """
    try:
        int_value = int(value)
    except ValueError:
        error = "Value %s is not integer" % value
        # report.log_str(error, priority=logging.ERROR)
        raise argparse.ArgumentTypeError(error)
    if int_value < 0:
        error = "Value %s is not positive integer" % value
        # report.log_str(error, priority=logging.ERROR)
        raise argparse.ArgumentTypeError(error)
    if max_limit and max_limit < int_value:
        error = "Value %s must be lower or equal to %s" % (value, max_limit)
        # report.log_str(error, priority=logging.ERROR)
        raise argparse.ArgumentTypeError(error)
    return int_value


def positive_nonzero_int(value):
    """
    Represents positive decimal number, 0 excluded
    :param value: string value to estimate
    :return: integer value, if can be converted into positive int else ArgumentTypeError
    """
    int_value = positive_int(value)
    if int_value == 0:
        error = "Value %s cannot be 0" % value
        # report.log_str(error, priority=logging.ERROR)
        raise argparse.ArgumentTypeError(error)
    return int_value


def nonempty_file(file_path):
    """
    Checks if the filename in input is non-empty file.
    :param file_path: str - filename to check
    :return: str - filename
    """
    if not os.path.exists(file_path):
        error = "File %s does not exist" % file_path
        # report.log_str(error, priority=logging.ERROR)
        raise argparse.ArgumentTypeError(error)
    if os.path.getsize(file_path) == 0:
        error = "File %s is empty" % file_path
        # report.log_str(error, priority=logging.ERROR)
        raise argparse.ArgumentTypeError(error)
    return file_path


def generate_files(cbs, prefix):
    """
    Generate output files for Werner's (circular) visualization
    :param cbs: pd.DataFrame - CBS smoothing of bins
    :param prefix: str - path and prefix for generated files
    :return: None
    """
    # save cbs line
    cbs_new = cbs.copy()
    cbs_new['start'] = cbs_new['start'].map(lambda x: '%10d' % x)
    cbs_new['end'] = cbs_new['end'].map(lambda x: '%10d' % x)
    cbs_new['length'] = cbs_new['length'].map(lambda x: '%5d' % x)
    cbs_new['effective_length'] = cbs_new['effective_length'].map(lambda x: '%5d' % x)
    filename_cbs = prefix + "_cbs.txt"
    cbs_new.to_csv(filename_cbs, sep='\t', float_format='%7.4f')


def get_significance(mean_bins, cbs, ind_good=None, window=20000, x_count=None, min_length=None):
    """
    Assign significance levels:
    5 significance levels: ok - green, intermediate - yellow, short cnv - orange, cnv - red, cnv(mother) - crimson
    :param mean_bins: ndarray - mean bin count (from many training samples)
    :param cbs: pd.DataFrame - CBS smoothing of bins
    :param ind_good: ndarray(bool) - numpy array indices with good bins
    :param window: int - size of bin
    :param x_count: int - number of X chromosomes for fetus
    :param min_length: int - minimal length for detection
    :return: pd.DataFrame, (float, float, float, float, float, float) - CBS smoothing of bins with significance info, expected red and magenta levels for autosomal, X and Y chromosome
    """
    mean_bins_autosomes = mean_bins[mean_bins['chromosome'] < 22]
    assert len(mean_bins_autosomes) == 144061  # 20000 bins on all autosomes

    if ind_good is None:
        ind_good = np.ones_like(mean_bins['bins_PCA'], dtype=bool)

    ind_good_autosomes = ind_good[:len(mean_bins_autosomes)]
    avg_bin = float(np.sum(mean_bins_autosomes['bins_PCA'][ind_good_autosomes], dtype=np.float128) / len(mean_bins_autosomes['bins_PCA'][ind_good_autosomes]))

    # significance levels
    yellow_sign_real = 0.25 * avg_bin / 2.0
    red_sign_real = 0.5 * avg_bin / 2.0
    crimson_sign_real = 1.0 * avg_bin / 2.0

    # significance for X and Y:
    yellow_sign_real_x = None
    red_sign_real_x = None
    crimson_sign_real_x = None
    avg_bin_x = None
    if x_count is not None and len(mean_bins) == len(ind_good):
        ind_x = (mean_bins['chromosome'] == 22) & ind_good
        avg_bin_x = float(np.sum(mean_bins['bins_PCA'][ind_x], dtype=np.float128) / len(mean_bins['bins_PCA'][ind_x]))
        yellow_sign_real_x = 0.25 * avg_bin_x / float(x_count)
        red_sign_real_x = 0.5 * avg_bin_x / float(x_count)
        crimson_sign_real_x = 1.0 * avg_bin_x / float(x_count)
    yellow_sign_real_y = None
    red_sign_real_y = None
    crimson_sign_real_y = None
    avg_bin_y = None
    if x_count is not None and len(mean_bins) == len(ind_good) and x_count == 1:
        ind_y = (mean_bins['chromosome'] == 23) & ind_good
        avg_bin_y = float(np.sum(mean_bins['bins_PCA'][ind_y], dtype=np.float128) / len(mean_bins['bins_PCA'][ind_y]))
        yellow_sign_real_y = 0.25 * avg_bin_y
        red_sign_real_y = 0.5 * avg_bin_y
        crimson_sign_real_y = 1.0 * avg_bin_y

    # significance levels with height_significance
    yellow_sign = max(0.1, 0.5 * yellow_sign_real)
    orange_sign = max(0.15, 0.5 * red_sign_real)
    red_sign = max(0.2, 0.5 * red_sign_real)
    crimson_sign = max(1.0, 0.75 * crimson_sign_real)
    # thus magenta (75% - 100%), red (25% - 75%), yellow (12.5% - 25%)

    # significance levels with height_significance X chromosome
    yellow_sign_x = None
    orange_sign_x = None
    red_sign_x = None
    crimson_sign_x = None
    if yellow_sign_real_x is not None:
        yellow_sign_x = max(0.1, 0.5 * yellow_sign_real_x)
        orange_sign_x = max(0.15, 0.5 * red_sign_real_x)
        red_sign_x = max(0.2, 0.5 * red_sign_real_x)
        crimson_sign_x = max(1.0, 0.75 * crimson_sign_real_x)

    # significance levels for Y:
    yellow_sign_y = None
    orange_sign_y = None
    red_sign_y = None
    crimson_sign_y = None
    if yellow_sign_real_y is not None:
        yellow_sign_y = max(0.1, 0.5 * yellow_sign_real_y)
        orange_sign_y = max(0.15, 0.5 * red_sign_real_y)
        red_sign_y = max(0.2, 0.5 * red_sign_real_y)
        crimson_sign_y = max(1.0, 0.75 * crimson_sign_real_y)

    # define signs arrays:
    signs = np.array([0, yellow_sign, orange_sign, red_sign, crimson_sign])
    signs_x = None
    if yellow_sign_x is not None:
        signs_x = np.array([0, yellow_sign_x, orange_sign_x, red_sign_x, crimson_sign_x])
    signs_y = None
    if yellow_sign_y is not None:
        signs_y = np.array([0, yellow_sign_y, orange_sign_y, red_sign_y, crimson_sign_y])

    # print significance levels
    print("Significance levels: %5.3f (%5.3f - yellow) %5.3f (%5.3f - orange) %5.3f (%5.3f - red) %5.3f (%5.3f - magenta) (%5.3f avg_bin)" % (
        yellow_sign, yellow_sign_real, orange_sign, red_sign_real, red_sign, red_sign_real, crimson_sign, crimson_sign_real, avg_bin))
    if yellow_sign_x is not None:
        print("Significance levels X: %5.3f yellow %5.3f orange %5.3f red %5.3f magenta (%5.3f avg_bin)" % (
            yellow_sign_x, orange_sign_x, red_sign_x, crimson_sign_x, avg_bin_x))
    if yellow_sign_y is not None:
        print("Significance levels Y: %5.3f yellow %5.3f orange %5.3f red %5.3f magenta (%5.3f avg_bin)" % (
            yellow_sign_y, orange_sign_y, red_sign_y, crimson_sign_y, avg_bin_y))

    def assign_significance(chromosome, level, length, min_length=None):
        """
        Assign significance level of single level and length.
        :param chromosome: int - chromosome number
        :param level: float - level of cbs smoothing
        :param length: int - length of interval in bins
        :param min_length: int - minimal length for detection
        :return: int<-1, 3> - significance level
        """
        used_signs = signs
        if chromosome == 22:
            used_signs = signs_x
        if chromosome == 23:
            used_signs = signs_y
        if used_signs is None:
            if abs(level) > 20:  # return orange if level is too much
                return 2
            return -1

        if min_length is None:
            min_length = 10 * 20000

        min_lengths = np.array([0, 2, 2, min_length / 20000, 10])

        for i in reversed(range(len(used_signs))):
            if abs(level) >= used_signs[i] and length >= min_lengths[i]:
                # if length >= 500 and i in [1, 2]:
                #    return i + 1
                return i
        return 0

    # calculate effective length
    effective_length = np.array(cbs['length'], dtype=int)
    for i, row in cbs.iterrows():
        ind_chr = mean_bins_autosomes['chromosome'] == row['chromosome']
        start_bin = int(row['start'] // window)
        if row['chromosome'] < 22:
            effective_length[i] = np.count_nonzero(ind_good_autosomes[ind_chr][start_bin: start_bin + int(row['length'])])
        if row['chromosome'] == 22:
            ind_chrX = mean_bins['chromosome'] == row['chromosome']
            effective_length[i] = np.count_nonzero(ind_good[ind_chrX][start_bin: start_bin + int(row['length'])])
        if row['chromosome'] == 23:
            ind_chrY = mean_bins['chromosome'] == row['chromosome']
            effective_length[i] = np.count_nonzero(ind_good[ind_chrY][start_bin: start_bin + int(row['length'])])
    cbs['effective_length'] = pd.Series(effective_length, index=cbs.index)

    # assign significance
    assign_significance_vec = np.vectorize(assign_significance)
    significant = assign_significance_vec(cbs["chromosome"], cbs["level"], cbs["effective_length"], min_length)

    colors = np.array(["cadetblue", "green", "yellow", "orange", "red", "magenta"])  # cadetblue is for X and Y chromosomes, when no significance is applied

    cbs['significant'] = pd.Series(significant, index=cbs.index)
    cbs['color'] = pd.Series(colors[significant + 1], index=cbs.index)

    def expected_ff(level, avg_bin, avg_bin_x, avg_bin_y, chromosome, x_count):
        """
        Calculate fraction that supposedly generated the detection.
        :param level: float - level of detection
        :param avg_bin: float - average bin count for autosomal chromosomes
        :param avg_bin_x: float - average bin count for X chromosome
        :param avg_bin_y: float - average bin count for Y chromosome
        :param chromosome: int - chromosome number
        :param x_count: int - count of X chromosomes
        :return: float - fraction that generated this detection
        """
        # gonosomal:
        if chromosome == 22:
            if avg_bin_x is None or x_count == 0:
                return 0.0
            else:
                if x_count == 2:
                    return np.abs(level / avg_bin_x * 2)
                elif x_count == 1:
                    return np.abs(level / avg_bin_x * 1)
                else:
                    return 0.0
        elif chromosome == 23:
            if avg_bin_y is None:
                return 0.0
            else:
                return np.abs(level / avg_bin_y)

        # autosomal:
        return np.abs(level / avg_bin * 2)

    # find expected ff for each detection:
    expected_ffs = [expected_ff(cbs['level'][i], avg_bin, avg_bin_x, avg_bin_y, cbs['chromosome'][i], x_count) for i in cbs.index]
    cbs['expected_ff'] = pd.Series(expected_ffs, index=cbs.index)

    return cbs, (red_sign_real, crimson_sign_real, red_sign_real_x, crimson_sign_real_x, red_sign_real_y, crimson_sign_real_y)


def extrem(arr):
    """
    Get extremal value.
    :param arr: list(float) - values
    :return: float - extremal value
    """
    if np.sum(arr, dtype=float) < 0:
        return min(arr)
    else:
        return max(arr)


def load_tsv(tsv_path):
    """
    Loads tab-separated table from the provided filepath.
    :param tsv_path: str - filepath to the table
    :return: None/pandas.DataFrame - loaded table or None
    """
    if not os.path.exists(tsv_path):
        print("Path %s does not exist" % tsv_path)
        return None
    try:
        return pd.read_csv(tsv_path, sep='\t')
    except IOError:
        pass

    print("Path %s cannot be loaded" % tsv_path)
    return None


def filtration(args, mean_stat, var_stat, chromosomes, excludey=False, excludey_variance=True):
    """
    Filtration of bins according to filters from the filtration tab in arguments. Either according to a file, or according to params.
    :param args: argparse arguments
    :param mean_stat: ndarray - numpy array with mean statistics of bins
    :param var_stat: ndarray - numpy array with variance statistics of bins
    :param chromosomes: ndarray - numpy array with chromosomes
    :param excludey: bool - exclude Y chromosome from filtration
    :param excludey_variance: bool - exclude variance on Y chromosome from filtration
    :return: pandas.DataFrame - gaps = bad bins in a nice DataFrame
    """
    assert var_stat is None or len(var_stat) == len(mean_stat), "len(var_stat) == %d, len(mean_stat) == %d, unequal!" % (len(var_stat), len(mean_stat))

    # filter according to params or file
    if args.filter_file is None:
        print('Filtering (max_var=%.1f min_mean=%.1f max_mean=%.1f) ' % (args.maximal_variance, args.minimal_mean, args.maximal_mean))
        # filter bins (remaining are True in ind):
        ind = (mean_stat > args.minimal_mean) & (mean_stat < args.maximal_mean)
        if var_stat is not None and chromosomes is not None:
            if excludey_variance:
                ind = ((var_stat < args.maximal_variance) | (chromosomes == 23)) & ind
            else:
                ind = (var_stat < args.maximal_variance) & ind

        # add Y chromosome if needed:
        if excludey:
            ind = (chromosomes == 23) | ind

    else:
        print('Filtering (according to gaps in %s) ' % args.filter_file)
        # load filtered bins
        gaps_dict = pd.read_csv(args.filter_file, index_col=None, header=0, sep='\t')

        # convert to indices
        ind = np.ones_like(mean_stat, dtype=bool)
        for i, row in gaps_dict.iterrows():
            start = int(row['full_start'] // args.window)
            # if not ind[start:start + row['length']].all():
            #    print(row, ind[start:start + row['length']])
            ind[start:start + int(row['length'])] = False

    print('Bins passed filtering %d/%d (%.1f%%)' % (np.count_nonzero(ind), len(ind), np.count_nonzero(ind) / len(ind) * 100.0))
    print('Bins filtered due to minimal mean    :', len(ind) - np.count_nonzero(mean_stat > args.minimal_mean))
    print('Bins filtered due to maximal mean    :', len(ind) - np.count_nonzero(mean_stat < args.maximal_mean))
    print('Bins filtered due to maximal variance:', len(ind) - np.count_nonzero(var_stat < args.maximal_variance))

    return ind


def normalize(args, bins_data, mean_stat, ind_good=None):
    """
    Normalize bins (subtract means) and save the file for R script.
    :param args: argparse arguments
    :param bins_data: ndarray - numpy array with read count per bins
    :param mean_stat: ndarray - numpy array with mean statistics of bins
    :param ind_good: ndarray(bool) - numpy array indices with good bins
    :return: pd.DataFrame, ndarray - bad bins (gaps), normalized bins
    """
    # modify the R input (remove some regions)
    bins_npy_mod = args.output + args.name + '_norm.npz'
    gaps_mod = args.output + args.name + '_gaps.txt'

    if ind_good is None:
        ind_good = np.ones_like(mean_stat, dtype=bool)

    assert len(ind_good) == len(mean_stat)

    # extract sums for normalization
    ind_autosomes = ind_good & (mean_stat['chromosome'] < 22)
    ind_x = ind_good & (mean_stat['chromosome'] == 22)
    ind_y = ind_good & (mean_stat['chromosome'] == 23)
    sum_mean_autosomes = np.sum(mean_stat['bins_PCA'][ind_autosomes], dtype=np.float128)
    sum_mean_x = np.sum(mean_stat['bins_PCA'][ind_x], dtype=np.float128)
    sum_mean_y = np.sum(mean_stat['bins_PCA'][ind_y], dtype=np.float128)
    bins_norm = bins_data[:len(ind_good)]
    bins_norm_PCA = bins_norm['bins_PCA'].astype(np.float128)
    bins_norm_PCA[~ind_good] = np.nan

    # normalize autosomes:
    print("Normalizing autosomes with mean of %10.2f" % sum_mean_autosomes)
    bins_norm_PCA[ind_autosomes] /= np.sum(bins_norm_PCA[ind_autosomes], dtype=np.float128) / sum_mean_autosomes

    # normalize gonosomes
    if sum_mean_x > 0:
        print("Normalizing X chrom.  with mean of %10.2f" % sum_mean_x)
        bins_norm_PCA[ind_x] /= np.sum(bins_norm_PCA[ind_x], dtype=np.float128) / sum_mean_x
    if sum_mean_y > 0:
        print("Normalizing Y chrom.  with mean of %10.2f" % sum_mean_y)
        sum_y = np.sum(bins_norm_PCA[ind_y], dtype=np.float128)
        # fix insufficient mapping (girls only)
        if sum_y < 1.0:
            sum_y = 1.0
        bins_norm_PCA[ind_y] *= sum_mean_y
        bins_norm_PCA[ind_y] /= sum_y

    # subtract means:
    bins_norm_PCA[ind_good] -= mean_stat['bins_PCA'][ind_good]

    # create whole normalization numpy data
    bins_norm['bins_loess'][~ind_good] = np.nan
    bins_norm['bins_PCA'] = bins_norm_PCA

    # save gaps
    save_bins(bins_norm, args.window, save_gaps=gaps_mod)

    # load gaps
    gaps_dict = pd.read_csv(gaps_mod, index_col=None, header=0, sep='\t')

    return gaps_dict, bins_norm


def save_bins(binarray, window, save_gaps=None):
    """
    Save bins and normalized bins to a numpy file together with their chromosome numbers.
    :param binarray: structured ndarray - fragment counts per bin for a single sample
    :param window: int - window size
    :param save_gaps: str - filename where to save the gaps
    :return: None
    """
    # save gaps:
    if save_gaps is not None:
        with open(save_gaps, 'w') as f:
            print('chromosome\tstart\tlength\tend\tfull_start', file=f)
            pos = 0
            last_ch = 0
            to_save = None
            full_pos = 0
            for ch, l, p in zip(binarray['chromosome'], binarray['bins_loess'], binarray['bins_PCA']):
                if last_ch != ch:
                    if to_save is not None:
                        print("%d\t%d\t%d\t%d\t%d" % (to_save[0], to_save[1], (pos - to_save[1]) // window, pos, to_save[2]), file=f)
                        to_save = None
                    pos = 0
                if np.isnan(l) or np.isnan(p):
                    if to_save is None:
                        to_save = [ch, pos, full_pos]
                    last_ch = ch
                    pos += window
                    full_pos += window
                    continue
                elif to_save is not None:
                    print("%d\t%d\t%d\t%d\t%d" % (to_save[0], to_save[1], (pos - to_save[1]) // window, pos, to_save[2]), file=f)
                    to_save = None
                last_ch = ch
                pos += window
                full_pos += window
            if to_save is not None:
                print("%d\t%d\t%d\t%d\t%d" % (to_save[0], to_save[1], (pos - to_save[1]) // window, pos, to_save[2]), file=f)


def gen_pos(chromosomes, window):
    """
    Generate positions for chromsoomes
    :param chromosomes: ndarray - chromosome numbers
    :param window: int - size of a bin
    :return: ndarray - positions
    """
    # for each chromosome generate postitions
    pos_arrays = []
    for chromosome in all_chromosomes:
        ind = chromosomes == chromosome
        x = np.arange(np.count_nonzero(ind)) * window
        pos_arrays.append(x)

    # return concatenated positions
    return np.hstack(pos_arrays)


def call_cbs_r(bins_norm, rscript, window):
    """
    Call R subroutine for computing Circular Binary Segmentation.
    :param bins_norm: structured ndarray - fragment counts per bin for a single sample
    :param rscript: str - path to the R library
    :param window: int - size of a bin
    :return: pd.DataFrame - dataframe with CBS results
    """
    # get valid bins
    ind_valid = ~np.isnan(bins_norm['bins_PCA'])

    # compute positions
    positions = gen_pos(bins_norm['chromosome'], window)

    # load the library
    rpy2.robjects.r.source(rscript)

    # convert data to R format
    cbs_data = rpy2.robjects.numpy2ri.py2ri(bins_norm['bins_PCA'][ind_valid])
    cbs_chromosome = rpy2.robjects.numpy2ri.py2ri(bins_norm['chromosome'][ind_valid])
    cbs_positions = rpy2.robjects.numpy2ri.py2ri(positions[ind_valid])

    # call the CBS function (this line takes ~11secs out of 12.5secs)
    ret = rpy2.robjects.r.CBS(cbs_data, cbs_chromosome, cbs_positions)

    # covert the result to dataframe and fix it
    cbs_r = rpy2.robjects.pandas2ri.ri2py(ret)
    cbs_r = cbs_r.drop(['ID'], axis=1)
    cbs_r.columns = ["chromosome", "start", "end", "length", "level"]
    cbs_r = cbs_r[["start", "length", "end", "chromosome", "level"]]
    cbs_r['length'] = np.divide(cbs_r['end'] - cbs_r['start'], 20000).astype(int)
    cbs_r.set_index(np.arange(len(cbs_r.index)), inplace=True)

    # return the dataframe
    return cbs_r


"""
Actual code starts here:
"""
if __name__ == "__main__":
    # print time of the start:
    start_time = datetime.now()
    print("Launching CNVCaller:     %s" % start_time.isoformat())

    # load arguments
    args = load_arguments()

    # extract sample name
    if not args.output.endswith("/"):
        args.output += "/"
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    # load data and normalize them
    print("Loading data             %s" % datetime.now().isoformat())
    bins_data = np.load(args.bins)
    if args.bins.endswith('.npz'):
        bins_data = bins_data['cnv']
    mean_data = np.load(args.mean)
    if args.mean.endswith('.npz'):
        mean_data = mean_data['means']

    # load statistics of bins - first try to extract it from mean file, then from stats file
    var_stat = mean_data['var_PCA']

    # filtration
    ind_filtered = None
    if args.filter_on:
        print("Filtering data           %s" % datetime.now().isoformat())
        ind_filtered = filtration(args, mean_data['bins_PCA'], var_stat, mean_data['chromosome'])

    # normalization
    print("Normalization data       %s" % datetime.now().isoformat())
    gaps_dict, bins_norm = normalize(args, bins_data, mean_data, ind_filtered)

    # call the CBS in R
    print("Calling CBS              %s" % datetime.now().isoformat())
    cbs_r = call_cbs_r(bins_norm, args.rscript, args.window)

    # filter cbs:
    print("Calculating significance %s" % datetime.now().isoformat())
    cbs, expected_levels = get_significance(mean_data, cbs_r, ind_filtered, x_count=args.x_count,
                                            min_length=args.min_length)

    # generate additional files
    print("Generating output into: %s %s" % (args.output + args.name, datetime.now().isoformat()))
    generate_files(cbs, args.output + args.name)

    # exit:
    end_time = datetime.now()
    print("Exiting CNVCaller:       %s" % end_time.isoformat())
    print("Time of run      :      ", (end_time - start_time))
