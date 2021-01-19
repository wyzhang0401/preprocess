import datetime
import re
import pandas as pd
import os
import numpy as np
import gzip
import multiprocessing
import time

# 定义全局参数
effectiveSymbols = "ACGT"
countOfEffectSymbols = len(effectiveSymbols)
# 这样取数组中的最后一个值，在下面的代码中，在值的数组中末尾添加平均值
n_index = -1
quaIndex = [0, -1, 1, -1, -1, -1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1]

'''
nucleotide 表示一元，二元，或三元核苷酸的串
length表示的是1,2,3，分别表示一元，二元或者三元核苷酸
功能：将4进制的数转换为10进制的数，10进制表示对应的串在值的数组中的下标
'''


def gen_index(nucleotide, length):
    index = 0
    if 'N' in nucleotide:
        index = n_index
    else:
        for i in range(length):
            diff = ord(nucleotide[i]) - ord('A')
            index = index * countOfEffectSymbols + quaIndex[diff]
    # print(nucleotide, index)
    return index


'''
读取分染色体存储的FASTA文件，一个文件一个序列，因此一个文件中只有一个>
genString参数表示的一个基因序列
返回的是将序列id和sequence分开的字典
'''


def read_single_fasta(genstring):
    id_seq = {}
    index = genstring.find("\n")
    id_seq[genstring[1:index]] = re.sub('\s+', "", genstring[index + 1:]).strip()
    return id_seq


'''
读取理化特性，xlsx文件
fname 文件名
'''


def read_property(fname):
    properties = pd.read_excel(fname)
    return properties


'''
取数组的平均值，并且保留3位小数
'''


def avg(array):
    mean = np.mean(array)
    mean_around = np.around(mean, decimals=3)
    return mean_around


'''
将基因组序列转换为值的序列, 边转换边保存到wiggle文件
参数：保存的文件名fout，理化特性值的数组pro_value，基因组序列gen_string，表示k值的参数k，染色体编号chrom
'''


def sequence_convert_index(gen_string, k):
    gen_str = gen_string.upper()
    length = len(gen_str) - k + 1
    print(length)
    value_index = np.zeros(length)
    for i in range(length):
        value_index[i] = gen_index(gen_str[i: i + k], k)
    return value_index


def save_wiggle(fout, value_array, chrom):
    # 在理化特性值的数组最后添加了一个值：该值是整个数组的平均值，如果kmers中出现N，index表示为-1，值的数组的倒数第一个值
    # value_array = np.array(np.append(pro_value, avg(pro_value)))[value_index.astype("int64")]
    # print(np.append(pro_value, avg(pro_value)))
    title = "fixedStep chrom=" + chrom + " start=1 step=1\n"
    title_bytes = title.encode(encoding="utf-8")
    value_str = "\n".join(str(i) for i in value_array)
    with gzip.open(fout, "w") as outfile:
        outfile.write(title_bytes)
        outfile.write(value_str.encode(encoding="utf-8"))
    outfile.close()


'''
对一个基因组序列（一个染色体文件）进行值的转换并且保存到wig.gz文件
参数：染色体文件路径 chrom_file, 理化特性值的集合property_values
'''


def get_values(chrom_file, property_values):
    print(multiprocessing.current_process().name + '-' + chrom_file)
    if os.path.exists(chrom_file):
        # 读fa.gz压缩文件
        file = gzip.GzipFile(chrom_file, "r")
        # decode bytes转换为string类型
        gen_string = file.read().decode()
        id_seq = read_single_fasta(gen_string)
        # 对读取到的每个序列（这里实际只有一个序列）
        for id in id_seq:
            seq_index_1 = sequence_convert_index(id_seq[id], 1)
            seq_index_2 = sequence_convert_index(id_seq[id], 2)
            seq_index_3 = sequence_convert_index(id_seq[id], 3)
            seq_index = np.array([seq_index_1, seq_index_2, seq_index_3], dtype=object)
            print(len(seq_index[0]), len(seq_index[1]), len(seq_index[2]))
            for i in range(len(property_values)):
                for j in range(len(property_values[i]["values"])):
                    fout_path = os.path.dirname(__file__) + "/values/" + os.path.basename(
                        os.path.dirname(chrom_file)) + "/" + property_values[i]["orista"] + "/" + property_values[i][
                                    "name"] + "/"
                    fout_name = id + "_" + property_values[i]["values"][j][0] + "_" + property_values[i]["values"][j][
                        1] + "_" + property_values[i]["orista"] + ".wig.gz"
                    value_index = seq_index[property_values[i]["k"] - 1]
                    pro_value = property_values[i]["values"][j][4:]
                    # 在理化特性值的数组最后添加了一个值：该值是整个数组的平均值，如果kmers中出现N，index表示为-1，值的数组的倒数第一个值
                    value_array = np.array(np.append(pro_value, avg(pro_value)))[value_index.astype("int64")]
                    save_wiggle(fout_path + fout_name, value_array, id)
                    time.sleep(2)
                time.sleep(1)
        file.close()
    # print(multiprocessing.current_process().name + '-' + chrom_file)


'''
预处理程序的主要函数
'''


if __name__ == '__main__':
    print("preprocess start... ...")

    # 创建processes个进程
    pool = multiprocessing.Pool(processes=1)
    start = datetime.datetime.now()
    base_dir = os.path.dirname(__file__)
    ori_base_path = os.path.join(base_dir, "property-values", "original")
    sta_base_path = os.path.join(base_dir, "property-values", "standard")

    monoDNAOri = read_property(os.path.join(ori_base_path, "MonoDNA-original.xlsx")).values
    monoDNASta = read_property(os.path.join(sta_base_path, "MonoDNA-standardized.xlsx")).values
    diDNAOri = read_property(os.path.join(ori_base_path, "DiDNA-original.xlsx")).values
    diDNASta = read_property(os.path.join(sta_base_path, "DiDNA-standardized.xlsx")).values
    triDNAOri = read_property(os.path.join(ori_base_path, "TriDNA-original.xlsx")).values
    triDNASta = read_property(os.path.join(sta_base_path, "TriDNA-standardized.xlsx")).values

    # diRNAOri暂时先不考虑
    # diRNASta暂时先不考虑

    # 理化特性值
    property_values = [{"k": 1, "name": "monoDNAOri", "orista": "ori", "values": monoDNAOri},
                       {"k": 2, "name": "diDNAOri", "orista": "ori", "values": diDNAOri},
                       {"k": 3, "name": "triDNAOri", "orista": "ori", "values": triDNAOri}]
    # property_values = [{"k": 1, "name": "monoDNASta", "orista": "sta", "values": monoDNASta},
    #                    {"k": 2, "name": "diDNASta", "orista": "sta", "values": diDNASta},
    #                    {"k": 3, "name": "triDNASta", "orista": "sta", "values": triDNASta}]
    # property_values = [{"k": 1, "name": "monoDNAOri", "orista": "ori", "values": monoDNAOri},
    #                    {"k": 2, "name": "diDNAOri", "orista": "ori", "values": diDNAOri},
    #                    {"k": 3, "name": "triDNAOri", "orista": "ori", "values": triDNAOri},
    #                    {"k": 1, "name": "monoDNASta", "orista": "sta", "values": monoDNASta},
    #                    {"k": 2, "name": "diDNASta", "orista": "sta", "values": diDNASta},
    #                    {"k": 3, "name": "triDNASta", "orista": "sta", "values": triDNASta}]
    # 基因组数据：分染色体处理数据
    hg38_path = os.path.join(base_dir, "genomes", "hg38")
    mm39_path = os.path.join(base_dir, "genomes", "mm39")
    saccer3_path = os.path.join(base_dir, "genomes", "saccer3")

    genome_path = mm39_path
    for filename in os.listdir(genome_path):
        gz_file_path = os.path.join(genome_path, filename)
        print(gz_file_path)
        # 去循环
        # get_values(gz_file_path, property_values)
        # 多进程对每个染色体文件进行处理
        pool.apply_async(get_values, (gz_file_path, property_values,))

    pool.close()  # 关闭进程池，表示不能在往进程池中添加进程
    pool.join()  # 等待进程池中的所有进程执行完毕，必须在close()之后调用
    print("Sub-process(es) done.")
    finish = datetime.datetime.now()
    print("执行花费的时间为：\t{}".format(finish - start))
    # get_values("F:\PycharmProjects\PreProcess\genomes\saccer3\chrI.fa.gz", property_values)

