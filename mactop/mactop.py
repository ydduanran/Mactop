import cooler
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
import iced
from scipy import stats
import networkx as nx
import markov_clustering as mc
pd.options.mode.chained_assignment = None
import pickle


def expand_matrix(orign_mat, probability_mat, power):
    expand_mat = orign_mat
    for i in range(power):
        expand_mat = expand_mat @ expand_mat
    return expand_mat

def inflate(expand_mat, inflation):
    inflate_mat = expand_mat
    inflate_mat = np.power(inflate_mat, inflation)
    inflate_column_sum = np.sum(inflate_mat, axis=0)
    inflate_mat_result = inflate_mat / inflate_column_sum
    return inflate_mat_result

def markov_cluster(adjacency_mat, inflation, power):
    column_sum = np.sum(adjacency_mat, axis=0)
    M1 = adjacency_mat / column_sum
    flag = True
    diedai = 0
    while flag:
        diedai += 1
        M2 = expand_matrix(M1, M1, power)
        M1 = inflate(M2, inflation)
        if (M1 == M2).all():
            flag = False
        if diedai == 100:
            print('do convergence')
            flag = False
    return M1

def result_split(inflation_mat):
    coo = coo_matrix(inflation_mat)
    re = np.unique(coo.row, return_counts=True)
    return re

def result_split2(inflation_mat):
    coo = coo_matrix(inflation_mat)
    re = np.unique(coo.row, return_counts=True)
    return re[1]

def manage_clustering(mat, inflation, power, cindex):
    mat[mat <= 0] = 0
    for aa in range(np.shape(mat)[0]):
        if mat[aa, aa] == 0:
            mat[aa, aa] = 1
    result = markov_cluster(mat, inflation, power)
    return result_split2(result)

def manage_clustering_for_clique(mat, inflation, power=1):
    mat[mat <= 0] = 0
    for aa in range(np.shape(mat)[0]):
        if mat[aa, aa] == 0:
            mat[aa, aa] = 1
    result = markov_cluster(mat, inflation, power)
    return result

def split_by_window(matrix, power, inflation, window_size, variance):
    split_begin = 0
    split_step = window_size
    all_result = []
    mat_length = int(np.shape(matrix)[0])
    clustingindex = 0
    while (split_begin + split_step) < mat_length:
        split_matrix = matrix[split_begin:split_begin + split_step, split_begin:split_begin + split_step]
        random_mat = Gaussian_matrix(split_matrix, variance)
        first_split = manage_clustering(random_mat, inflation, power, clustingindex)
        clustingindex += 1
        first_split = first_split[:-1]
        if len(first_split) == 0:
            first_split = [window_size]
        all_result.append(first_split)
        split_begin = split_begin + sum(first_split)

    split_matrix = matrix[split_begin:mat_length, split_begin:mat_length]
    random_mat = Gaussian_matrix(split_matrix, variance)
    first_split = manage_clustering(random_mat, inflation, power, clustingindex)
    all_result.append(first_split)
    filter_result = []
    filter_result.append(0)
    begin = 0
    for aa in all_result:
        for bb in aa:
            begin += bb
            filter_result.append(begin)
    return filter_result

def get_consistent_matrix(results, mat_length):
    zeores_mat = np.zeros(shape=(mat_length, mat_length))
    for re in results:
        relist = re[:-1]
        begin_bin = 0
        for floatvalue in relist:
            value = int(floatvalue)
            zeores_mat[begin_bin:value, begin_bin:value] += 1
            begin_bin = value
    return zeores_mat

def Gaussian_matrix(hicmatrix, variance):
    hic_matrix = hicmatrix
    beta = variance
    beta_matrix = np.random.uniform(low=-beta, high=beta, size=np.shape(hic_matrix))
    result = np.floor(hic_matrix + np.multiply(beta_matrix, hic_matrix))
    result = result.astype(np.int)
    triuma = np.triu(result, -1)
    final = triuma + triuma.T
    return final

def split_by_window(matrix, power, inflation, window_size, variance):
    split_begin = 0
    split_step = window_size
    all_result = []
    mat_length = int(np.shape(matrix)[0])
    clustingindex = 0
    while (split_begin + split_step) < mat_length:
        split_matrix = matrix[split_begin:split_begin + split_step, split_begin:split_begin + split_step]
        random_mat = Gaussian_matrix(split_matrix, variance)
        first_split = manage_clustering(random_mat, inflation, power, clustingindex)
        clustingindex += 1
        first_split = first_split[:-1]
        if len(first_split) == 0:
            first_split = [window_size]
        all_result.append(first_split)
        split_begin = split_begin + sum(first_split)

    split_matrix = matrix[split_begin:mat_length, split_begin:mat_length]
    random_mat = Gaussian_matrix(split_matrix, variance)
    first_split = manage_clustering(random_mat, inflation, power, clustingindex)
    all_result.append(first_split)
    filter_result = []
    filter_result.append(0)
    begin = 0
    for aa in all_result:
        for bb in aa:
            begin += bb
            filter_result.append(begin)
    return filter_result

def consistent_boundary_score(matrix):
    mat = matrix
    window_size = 5
    mat_length = np.shape(mat)[0]
    relist = []
    for bin_index in range(mat_length):
        if bin_index < window_size:
            # print(bin_index)
            relist.append(0)
        if bin_index >= window_size and bin_index < mat_length - window_size:
            # print(bin_index, mat[bin_index][bin_index])
            # print(mat[bin_index-window_size:bin_index, bin_index+1:bin_index+window_size+1])
            sum = np.sum(mat[bin_index - window_size:bin_index, bin_index + 1:bin_index + window_size + 1])
            relist.append(sum)
        if bin_index > mat_length - window_size:
            relist.append(0)

        bin_index += 1

    min_list = []
    value_list = []
    for index in range(1, len(relist) - 1):
        if relist[index] < 2500:
            if relist[index] < relist[index - 1] and relist[index] < relist[index + 1]:
                min_list.append(index + 1)
                value_list.append(relist[index])
            if relist[index] < relist[index - 1] and relist[index] == relist[index + 1]:
                min_list.append(index + 1)
                value_list.append(relist[index])

    return relist, min_list, value_list

def mactop_for_single_chrom(input_mat, mat_len, inflation=1.7, window_size=50, variance=0.2):
    power = 1
    hic_mat = input_mat
    print('begin cluster')
    results = []
    for index in range(100):
        result = split_by_window(hic_mat, power, inflation, window_size, variance)
        results.append(result)
    print('get consistent matrix')
    consistentmat = get_consistent_matrix(results, mat_len)
    relist, min_list, value_list = consistent_boundary_score(consistentmat)
    return relist, min_list, value_list, consistentmat

def get_tad_list(hic_result):
    tad_list = []
    begin = 0
    for vv in hic_result:
        end = int(vv)
        tad_list.append(np.array([begin, end]))
        begin = end
    return tad_list

def significance(mat, tad_1, tad_2):
    b_a = tad_1[0]
    b_b = tad_1[1]
    b_c = tad_2[0]
    b_d = tad_2[1]
    fisrt_digonal = b_c - b_b
    secend_digonal = b_d - b_a
    back_ground_dis = np.array([])
    for dia_indx in range(fisrt_digonal, secend_digonal + 1):
        diag_vv = np.diagonal(mat, offset=dia_indx)
        back_ground_dis = np.concatenate((back_ground_dis, diag_vv))

    true_value_dis = mat[b_c:b_d + 1, b_a:b_b + 1]
    true_value_dis = true_value_dis.flatten()
    t_statistic, p_value = stats.ttest_ind(true_value_dis, back_ground_dis)
    return t_statistic, p_value

def significance2(mat, tad_1, tad_2):
    b_a = tad_1[0]
    b_b = tad_1[1]
    b_c = tad_2[0]
    b_d = tad_2[1]
    fisrt_digonal = b_c - b_b
    secend_digonal = b_d - b_a
    back_ground_dis = np.array([])
    for dia_indx in range(fisrt_digonal, secend_digonal + 1):
        diag_vv = np.diagonal(mat, offset=dia_indx)
        back_ground_dis = np.concatenate((back_ground_dis, diag_vv))

    true_value_dis = mat[b_c:b_d + 1, b_a:b_b + 1]
    true_value_dis = true_value_dis.flatten()
    t_statistic, p_value = stats.ttest_ind(true_value_dis, back_ground_dis)
    return t_statistic, p_value

def nearest_naber(c_list, knumber):
    naber_list = []
    c_list = np.array(c_list)
    c_list_index = np.where(c_list > 0)[0]
    c_list_value = c_list[c_list_index]
    if len(c_list_index) <= knumber:
        return c_list_index
    else:
        max_knumber_index = np.argpartition(c_list_value, -knumber)[-knumber:]
        max_knumber = c_list_index[max_knumber_index]
        return max_knumber

def cluster_to_heatmap(concatemers_data, begin_bin, end_bin, concatemers_cluster_id, resulotion):
    concatemers_copy = concatemers_data.copy()
    concatemers_copy['align2_start'] = concatemers_copy['align2_start'] // resulotion
    concatemers_copy['align1_start'] = concatemers_copy['align1_start'] // resulotion
    concatemers_copy['align2_start'] = concatemers_copy['align2_start'] - begin_bin
    concatemers_copy['align1_start'] = concatemers_copy['align1_start'] - begin_bin

    mat_length = end_bin - begin_bin
    mat = np.zeros((mat_length + 1, mat_length + 1))

    read_idx_list = concatemers_cluster_id
    for read_idx in read_idx_list:
        data_read_idx = concatemers_copy[concatemers_copy['read_idx'] == read_idx]
        bins_set = set(data_read_idx['align1_start'].tolist() + data_read_idx['align2_start'].tolist())
        bins_list = list(bins_set)
        coordinates = []
        for i in range(len(bins_list)):
            for j in range(i, len(bins_list)):
                if i != j:
                    coordinates.append((bins_list[i], bins_list[j]))
        same_align = data_read_idx[data_read_idx['align1_start'] == data_read_idx['align2_start']]
        same_align_list = list(same_align['align1_start'])
        for single_align in same_align_list:
            coordinates.append((single_align, single_align))
        for coordinate in coordinates:
            mat[coordinate[0], coordinate[1]] += 1
            mat[coordinate[1], coordinate[0]] += 1
    return mat

def original_matrix_to_ice_matrix(hic_matrix):
    b = hic_matrix[:]
    filter_low_counts_perc = 0.02
    filter_high_counts_perc = 0
    mat = b.astype('float')
    if filter_low_counts_perc != 0:
        mat = iced.filter.filter_low_counts(mat, percentage=filter_low_counts_perc, copy=False, sparsity=False,
                                            verbose=False)
    if filter_high_counts_perc != 0:
        mat = iced.filter.filter_high_counts(mat, percentage=filter_high_counts_perc, copy=False)

    norm_mat = iced.normalization.ICE_normalization(mat, SS=None, max_iter=100, eps=0.0001, copy=True, norm='l1',
                                                    verbose=0,
                                                    output_bias=False, total_counts=None, counts_profile=None)

    return norm_mat

def parse_chr_num_from_str(chr_num_str):
    parsed_list = []
    items = chr_num_str.split(',')
    for item in items:
        if '-' in item:
            start, end = item.split('-')
            parsed_list.extend(range(int(start), int(end) + 1))
        else:
            parsed_list.append(int(item))
    return parsed_list

def get_boundary_type(value_list, min_list,stable_cut_of):
    value_cut_off_90 = np.percentile(value_list, 90)
    value_cut_off = np.percentile(value_list, stable_cut_of)
    count_index = []
    result_list = []
    for index in range(len(value_list)):
        if value_list[index] < value_cut_off_90:
            count_index.append(value_list[index])
            result_list.append(min_list[index])
    begin_list = []
    begin_type_list = []
    end_list = []
    end_type_list = []

    tad_list = get_tad_list(result_list)
    for tad_index in range(len(tad_list) - 1):
        tad_begin = tad_index
        tad_end = tad_index + 1
        begin_list.append(tad_list[tad_index][0])
        end_list.append(tad_list[tad_index][1])

        tad_begin_count = count_index[tad_begin]
        tad_end_count = count_index[tad_end]

        if tad_begin_count < value_cut_off:
            begin_type_list.append('stable')
        else:
            begin_type_list.append('dynamic')

        if tad_end_count < value_cut_off:
            end_type_list.append('stable')
        else:
            end_type_list.append('dynamic')

    return begin_list, begin_type_list, end_list, end_type_list, result_list, count_index

def get_community_list(tad_df, dense_mat, inflation = 5):
    bbls = tad_df['begin'].tolist()
    eels = tad_df['end'].tolist()
    tad_list = []
    for index in range(len(bbls)):
        tad_list.append([bbls[index], eels[index]])
    diag_mat = np.copy(dense_mat)
    node_one_list = []
    node_two_list = []
    count_list = []

    for index in range(len(tad_list)):
        for index2 in range(index + 1, len(tad_list)):
            tad_one = tad_list[index]
            tad_two = tad_list[index2]
            tad_one_begin = tad_one[0]
            tad_one_end = tad_one[1]
            tad_two_begin = tad_two[0]
            tad_two_end = tad_two[1]
            ot_conut = diag_mat[tad_two_begin:tad_two_end + 1, tad_one_begin:tad_one_end + 1]
            weight_ori = np.sum(ot_conut) / ((tad_two_end - tad_two_begin) * (tad_one_end - tad_one_begin))
            t_statistic, p_value = significance(diag_mat, tad_one, tad_two)
            if t_statistic > 0 and p_value < 0.001:
                node_one_list.append(index)
                node_two_list.append(index2)
                count_list.append(weight_ori)
    df = pd.DataFrame({'node_one': node_one_list, 'node_two': node_two_list, 'weight': count_list})

    G = nx.from_pandas_edgelist(df, 'node_one', 'node_two', edge_attr='weight')
    adj_mat = nx.to_numpy_matrix(G)
    ddre = manage_clustering_for_clique(adj_mat,inflation)
    re = result_split(ddre)
    community_for_each_tad = []
    for index in range(len(tad_list)):
        for lines in re[0]:
            position = ddre[lines]
            position = np.where(position != 0)
            if index in position[1]:
                community_for_each_tad.append(position[1])
    return community_for_each_tad

def get_mat_from_muti_reads(ori_data, chr_num, resolution):
    ori_data_not_change = ori_data.copy()
    select_chr_reads = ori_data_not_change[(ori_data_not_change['align1_chrom'] == 'chr{0}'.format(chr_num)) & (
            ori_data_not_change['align2_chrom'] == 'chr{0}'.format(chr_num))]
    bin_range_read_new = select_chr_reads[
        ['read_idx', 'align1_start', 'align1_end', 'align2_start', 'align2_end']]
    bin_range_read_new['align1_length'] = bin_range_read_new['align1_end'] - bin_range_read_new['align1_start']
    bin_range_read_new['align2_length'] = bin_range_read_new['align2_end'] - bin_range_read_new['align2_start']

    bin_range_read_change = bin_range_read_new.copy()
    concatemers_adj_reselution = resolution
    bin_range_read_change['align1_start'] = bin_range_read_change['align1_start'] // concatemers_adj_reselution
    bin_range_read_change['align2_start'] = bin_range_read_change['align2_start'] // concatemers_adj_reselution
    bin_range_read_change['align1_end'] = bin_range_read_change['align1_end'] // concatemers_adj_reselution
    bin_range_read_change['align2_end'] = bin_range_read_change['align2_end'] // concatemers_adj_reselution

    concatemers_copy = bin_range_read_change.copy()
    number = concatemers_copy[['align1_start', 'align2_start']].min(axis=1).min()
    concatemers_copy['align1_start'] = concatemers_copy['align1_start'] - number
    concatemers_copy['align2_start'] = concatemers_copy['align2_start'] - number

    max_number = concatemers_copy[['align1_start', 'align2_start']].max(axis=1).max().max()
    jacard_list_matrix = np.zeros((max_number + 1, max_number + 1))
    for index in range(len(concatemers_copy)):
        align1_start = concatemers_copy.iloc[index]['align1_start']
        align2_start = concatemers_copy.iloc[index]['align2_start']
        jacard_list_matrix[align1_start, align2_start] += 1
        jacard_list_matrix[align2_start, align1_start] += 1
    return jacard_list_matrix

def call_tad(dense_mat, inflation=1.7, window_size=50, variance=0.2,stable_cut_of = 30):
    tad_df = pd.DataFrame(columns=['begin', 'begin_type', 'end', 'end_type'])
    relist, min_list, value_list, consistentmat = mactop_for_single_chrom(dense_mat,
                                                                          np.shape(dense_mat)[0],
                                                                          inflation=inflation,
                                                                          window_size=window_size,
                                                                          variance=variance)
    begin_list, begin_type_list, end_list, end_type_list, result_list, count_index = get_boundary_type(value_list,
                                                                                                       min_list,stable_cut_of)
    tad_df['begin'] = begin_list
    tad_df['begin_type'] = begin_type_list
    tad_df['end'] = end_list
    tad_df['end_type'] = end_type_list

    return tad_df,consistentmat

def call_community(tad_df, dense_mat, inflation = 5):
    community_for_each_tad = get_community_list(tad_df, dense_mat, inflation)
    return community_for_each_tad

def call_chromunity(select_chr_reads, resolution, tad_df, inflation=1.7):
    bbls = tad_df['begin'].tolist()
    eels = tad_df['end'].tolist()
    tad_list = []
    for index in range(len(bbls)):
        tad_list.append([bbls[index], eels[index]])
    tad_chromunity = {}
    for tad_indexs in range(len(tad_list)):
        begin = tad_list[tad_indexs][0]
        end = tad_list[tad_indexs][1]
        bin_range_begin = begin
        bin_range_end = end

        data_chr_unchange = select_chr_reads.copy()
        bin_range_read = data_chr_unchange[(data_chr_unchange['align1_start'] > bin_range_begin * resolution) & (
                data_chr_unchange['align1_start'] < bin_range_end * resolution) & (data_chr_unchange[
                                                                                       'align2_start'] > bin_range_begin * resolution) & (
                                                   data_chr_unchange[
                                                       'align2_start'] < bin_range_end * resolution)]

        bin_range_read = bin_range_read[
            bin_range_read['read_idx'].map(bin_range_read['read_idx'].value_counts()) > 1]

        bin_range_read_new = bin_range_read[
            ['read_idx', 'align1_start', 'align1_end', 'align2_start', 'align2_end']]
        bin_range_read_new['align1_length'] = bin_range_read_new['align1_end'] - bin_range_read_new['align1_start']
        bin_range_read_new['align2_length'] = bin_range_read_new['align2_end'] - bin_range_read_new['align2_start']

        bin_range_read_change = bin_range_read_new.copy()
        concatemers_adj_resolution = resolution
        bin_range_read_change['align1_start'] = bin_range_read_change['align1_start'] // concatemers_adj_resolution
        bin_range_read_change['align2_start'] = bin_range_read_change['align2_start'] // concatemers_adj_resolution
        bin_range_read_change['align1_end'] = bin_range_read_change['align1_end'] // concatemers_adj_resolution
        bin_range_read_change['align2_end'] = bin_range_read_change['align2_end'] // concatemers_adj_resolution

        concatemers_copy = bin_range_read_change.copy()
        number = concatemers_copy[['align1_start', 'align2_start']].min(axis=1).min()
        concatemers_copy['align1_start'] = concatemers_copy['align1_start'] - number
        concatemers_copy['align2_start'] = concatemers_copy['align2_start'] - number

        # print(concatemers_copy[['align1_start','align2_start']])
        if concatemers_copy[['align1_start', 'align2_start']].empty:
            continue
        max_number = concatemers_copy[['align1_start', 'align2_start']].max(axis=1).max().max()

        concatemer_list = []
        read_ids = concatemers_copy['read_idx'].unique().tolist()
        for read_idx in read_ids:
            data_read_idx = concatemers_copy[concatemers_copy['read_idx'] == read_idx]
            bins_set = set(data_read_idx['align1_start'].tolist() + data_read_idx['align2_start'].tolist())
            concatemer_list.append(bins_set)

        jacard_list_matrix = []
        for i in range(len(concatemer_list)):
            jaccard_list = []
            for j in range(len(concatemer_list)):
                if i != j:
                    union = concatemer_list[i] | concatemer_list[j]
                    intersection = concatemer_list[i] & concatemer_list[j]
                    javale = len(intersection) / len(union)
                    if javale > 0.2:
                        jaccard_list.append(javale)
                    else:
                        jaccard_list.append(0)
                else:
                    jaccard_list.append(0)
            jacard_list_matrix.append(jaccard_list)

        jacard_list_matrix = np.array(jacard_list_matrix)

        adj_mat = np.zeros(np.shape(jacard_list_matrix))

        for i in range(len(jacard_list_matrix)):
            for j in range(i + 1, len(jacard_list_matrix)):
                first_naber = nearest_naber(jacard_list_matrix[i], 25)
                secend_naber = nearest_naber(jacard_list_matrix[j], 25)
                intersection = set(first_naber) & set(secend_naber)
                if len(intersection) > 3:
                    adj_mat[i, j] = len(intersection)
                    adj_mat[j, i] = len(intersection)
        for i in range(len(adj_mat)):
            adj_mat[i, i] = 1

        mc_copy = adj_mat.copy()
        adj_G = nx.from_numpy_matrix(mc_copy)
        # avg_clustering = nx.average_clustering(adj_G)
        # if avg_clustering > 0.7:
        #     up_cut_off.append(tad_indexs)
        # else:
        #     down_cut_off.append(tad_indexs)
        # avg_cluster_list.append(avg_clustering)
        clusters = mc.get_clusters(mc.run_mcl(mc_copy, inflation))
        cris = []
        cluster_length_list = []
        for cc in clusters:
            cluster_length_list.append(len(cc))
            if len(cc) > 20:
                cluster_read_idx = []
                for c in cc:
                    cluster_read_idx.append(read_ids[c])
                cris.append(cluster_read_idx)

        chromunitys = []
        for i in range(len(cris)):
            finmat = cluster_to_heatmap(bin_range_read_new, bin_range_begin, bin_range_end, cris[i],
                                        resolution)
            chromunitys.append(finmat)
        tad_chromunity[tad_indexs] = chromunitys
    return tad_chromunity