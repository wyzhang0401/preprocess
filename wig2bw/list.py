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


def get_values(chrom_file):
    print(multiprocessing.current_process().name + '-' + chrom_file)
    print(gz_file_path)


if __name__ == '__main__':
    print("preprocess start... ...")

    # 创建processes个进程
    pool = multiprocessing.Pool(processes=4)
    start = datetime.datetime.now()
    base_dir = os.path.dirname(__file__)

    # 基因组数据：分染色体处理数据
    saccer3_path = os.path.join(base_dir, "genomes", "saccer3")

    # 重新运行，生成文件之前一定要先把原来的文件删除，因为写文件使用的是追加
    genome_path = saccer3_path
    for filename in os.listdir(genome_path):
        gz_file_path = os.path.join(genome_path, filename)
        # 去循环
        # 多进程对每个染色体文件进行处理
        pool.apply_async(get_values, (gz_file_path,))

    pool.close()  # 关闭进程池，表示不能在往进程池中添加进程
    pool.join()  # 等待进程池中的所有进程执行完毕，必须在close()之后调用
    print("Sub-process(es) done.")
    finish = datetime.datetime.now()
    print("执行花费的时间为：\t{}".format(finish - start))


