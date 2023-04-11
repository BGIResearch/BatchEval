#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2023/4/11 11:13
# @Author  : zhangchao
# @File    : demo.py
# @Software: PyCharm
# @Email   : zhangchao5@genomics.cn
import os
import sys
import argparse
import scanpy as sc
from warnings import filterwarnings

from BatchQC import batch_qc

sys.path.append(os.getcwd())

sys.setrecursionlimit(1000000)
filterwarnings(action="ignore")


def main(args):
    data = [sc.read_h5ad(os.path.join(args.data_path, file)) for file in os.listdir(args.data_path) if file.endswith("h5ad")]
    print("UserWarning: The input gene expression matrix must be the original gene capture, "
          "otherwise the final quality control accuracy will be affected.")
    batch_qc(*data,
             qc_mode=args.mode,
             adjust_method=args.adjust_method,
             norm_log=args.norm_log,
             is_scale=args.is_scale,
             n_pcs=args.n_pcs,
             n_neighbors=args.n_neighbors,
             batch_key=args.batch_key,
             position_key=args.position_key,
             condition=args.condition,
             use_rep=args.use_rep,
             count_key=args.count_key,
             celltype_key=args.celltype_key,
             report_path=args.save_path,
             gpu=args.gpu)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="BatchQC Pipeline Running...")
    parser.add_argument("--data_path", default="./demo_data", type=str, help="dataset path")
    parser.add_argument("--mode", default="adjust", type=str, help="select running mode, only support 'raw' and 'adjust'")
    parser.add_argument("--adjust_method", default="harmony", type=str,
                        help="select adjust mode, it works for 'mode = adjust'")
    parser.add_argument("--norm_log", default=True, type=bool,
                        help="Whether to preprocess data. 'sc.pp.normalization()', 'sc.pp.log1p()'")
    parser.add_argument("--is_scale", default=True, type=bool,
                        help="Whether to preprocess data. 'sc.pp.scale()'")
    parser.add_argument("--n_pcs", default=50, type=int,
                        help="Number of principal components retained in PCA. default, 50.")
    parser.add_argument("--n_neighbors", default=30, type=int,
                        help="Calculate the nearest neighbors of a local area. default, 100.")
    parser.add_argument("--batch_key", default="batch", type=str, help="Label the data batches.")
    parser.add_argument("--use_rep", default="X_pca", type=str,
                        help="the embedding use to removal batch effects, it works for 'mode = adjust'")
    parser.add_argument("--position_key", default="X_umap", type=str,
                        help="Compute the coordinate space of the nearest neighbor.")
    parser.add_argument("--condition", default=None, type=str,
                        help="Label the experimental conditions. By default, the experimental conditions for each data are different.")
    parser.add_argument("--count_key", default="log1p_total_counts", type=str)
    parser.add_argument("--celltype_key", default=None, type=str)
    parser.add_argument("--gpu", default="0", type=str)
    parser.add_argument("--save_path", default="./output", type=str)
    args = parser.parse_args()

    main(args)
