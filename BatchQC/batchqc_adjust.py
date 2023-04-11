#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2023/4/11 10:51
# @Author  : zhangchao
# @File    : batchqc_adjust.py
# @Software: PyCharm
# @Email   : zhangchao5@genomics.cn
import os
import os.path as osp
import pkgutil
import pwd
import time
import scanpy as sc
import pandas as pd
from collections import defaultdict
from typing import Union
from anndata import AnnData
from lxml import etree

from test import *
from module import domain_variance_score
from utils import pca_lowrank, embed_text, embed_tabel, embed_table_imgs


def batchqc_adjust(
        merge_data: AnnData,
        adjust_method: str = "harmony",
        norm_log: bool = False,
        is_scale: bool = False,
        n_pcs: int = 50,
        use_rep: str = "X_pca",
        n_neighbors: int = 100,
        batch_key: str = "batch",
        position_key: str = "X_umap",
        count_key: str = "total_counts",
        celltype_key: Union[str, None] = None,
        report_path: str = "./",
        gpu: Union[str, int] = 0) -> dict:
    """BatchQC Adjust Dataset Pipeline

    Parameters
    -----------------
    merge_data: `AnnData`
        Data matrix with rows for cells and columns for genes.
    adjust_method: `str`
        batch effect removal method, default 'Harmony'
    norm_log: 'bool'
        Whether to preprocess data. 'sc.pp.normalization()', 'sc.pp.log1p()'
    is_scale: 'bool'
         Whether to preprocess data. 'sc.pp.scale()'
    n_pcs: 'int'
        Number of principal components retained in PCA. default, 50.
    use_rep: 'str'
    n_neighbors: 'int'
        Calculate the nearest neighbors of a local area. default, 100.
    batch_key: 'str'
        Label the data batches.
    position_key: 'str'
        Compute the coordinate space of the nearest neighbor.
    count_key: 'str'
    celltype_key: 'str'
    report_path: 'str'
    gpu: 'str', 'int'

    Returns
    -----------------
    """
    output_dict = defaultdict()

    n_batch = merge_data.obs[batch_key].cat.categories.size
    if merge_data.raw is None:
        merge_data.raw = merge_data
    else:
        merge_data = merge_data.raw.to_adata()
    if count_key not in merge_data.obs_keys():
        sc.pp.calculate_qc_metrics(merge_data, inplace=True)
    if norm_log:
        sc.pp.normalize_total(merge_data)
        sc.pp.log1p(merge_data)
    if is_scale:
        sc.pp.scale(merge_data, zero_center=False, max_value=10)

    if use_rep not in merge_data.obsm_keys():
        pca_lowrank(merge_data, use_rep=None, n_component=n_pcs)

    if adjust_method == "harmony":
        output_dict["adjust_method"] = adjust_method
        sc.external.pp.harmony_integrate(merge_data, key=batch_key, basis=use_rep)
        output_dict[f"dataset_{adjust_method}"] = merge_data.copy()
        sc.pp.neighbors(merge_data, use_rep="X_pca_harmony")
        sc.tl.umap(merge_data)

        domain_df = domain_variance_score(merge_data,
                                          input_dims=merge_data.obsm["X_pca_harmony"].shape[1],
                                          n_batch=n_batch,
                                          use_rep="X_pca_harmony",
                                          batch_key=batch_key,
                                          batch_size=4096,
                                          gpu=gpu,
                                          save_path=report_path)
        output_dict["table"] = {"domain": domain_df}

        metric_dict = metric_score(merge_data, n_neighbor=n_neighbors, batch_key=batch_key, metric_pos=position_key,
                                   celltype_key=celltype_key)
        output_dict["table"].update(metric_dict)
        pcr_raw = pca_regression(merge_data, n_pcs=n_pcs, batch_key=batch_key, embed_key=None)
        pcr_adjust = pca_regression(merge_data, n_pcs=n_pcs, batch_key=batch_key, embed_key="X_pca_harmony")

        pcr_score = (pcr_raw - pcr_adjust) / pcr_raw
        pcr_df = pd.DataFrame(data={"pcr raw": pcr_raw,
                                    f"pcr {adjust_method}": pcr_adjust,
                                    "pcr score": pcr_score if pcr_score > 0 else 0},
                              index=["pca_regression"])
        output_dict["table"]["pcr_df"] = pcr_df

        heatmap_gene_src = sample_heatmap(merge_data, feat_key="X_pca_harmony", metric="correlation",
                                          batch_key=batch_key)
        umap_batch_src = umap_plot(merge_data, visualize_key=batch_key)
        if celltype_key is not None:
            umap_type_src = umap_plot(merge_data, visualize_key=celltype_key)
            output_dict[f"imgs"]["umap_type"] = umap_type_src
        joint_srcs = joint_plot(merge_data, batch_key=batch_key, use_rep="X_pca_harmony")
        output_dict[f"imgs"] = {"heatmap": heatmap_gene_src, "umap_batch": umap_batch_src, "joint": joint_srcs}

        batch_score = output_dict["table"]["domain"].iloc[0]["Accept Rate"] * \
                      output_dict["table"]["kbet_df"].iloc[0]["Accept Rate"] * \
                      output_dict["table"]["lisi_df"].iloc[0]["LISI Mean"] / n_batch

        is_scale = False
        if output_dict["table"]["kbet_df"].iloc[0]["95% P Value"] < 0.05:
            threshold = 0.05
        elif output_dict["table"]["kbet_df"].iloc[0]["95% P Value"] < 0.1:
            is_scale = True
            threshold = 0.1
        else:
            threshold = 0.1

        if batch_score < threshold:
            conclusion = "Need to do batch effect removal."
        elif is_scale:
            conclusion = "The data difference is small. Recommend approach: 'z-score'."
        else:
            conclusion = "Don't have to do batch effects removal."

        summary_df = pd.DataFrame(
            data={"score": f"{batch_score:.4f}", "adaptive threshold": threshold, "conclusion": conclusion},
            index=["result"])

        output_dict["table"]["summary"] = summary_df

        generate_report(output_dict, report_path)

    return output_dict


def generate_report(data_dict: dict,
                    save_path: str) -> None:
    """Generate BatchQC Report

    Parameters
    -----------------
    data_dict: 'dict'
    save_path: 'str'
    """
    html = etree.HTML(pkgutil.get_data("template", "report_template_adjust.html").decode())

    # -------- set username & run time --------
    embed_text(html, pos="h4", name="username", text=f"Report By: {pwd.getpwuid(os.getuid())[0]}")
    embed_text(html, pos="h5", name="runtime",
               text=f"Report Time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")

    # -------- insert table --------
    embed_text(html, pos="h3", name="adj-name", text=f"Adjust Method: {data_dict['adjust_method']}")
    embed_tabel(data_dict["table"]['domain'], html, pos="h4", name="domain")
    embed_tabel(data_dict["table"]['kbet_df'], html, pos="h4", name="kbet")
    embed_tabel(data_dict["table"]['lisi_df'], html, pos="h4", name="lisi")
    embed_tabel(data_dict["table"]['ksim_df'], html, pos="h4", name="ksim")
    embed_tabel(data_dict["table"]["pcr_df"], html, pos="h4", name="pcr")
    embed_tabel(data_dict["table"]['summary'], html, pos="h4", name="conclusion")

    # -------- insert images --------
    src_dict = {
        "HeatMap": data_dict["imgs"]["heatmap"],
        "Joint": data_dict["imgs"]["joint"],
        "UMAP-Batch": data_dict["imgs"]["umap_batch"]
    }
    if "umap_type" in data_dict["imgs"].keys():
        src_dict.update({"UMAP-Type": data_dict["imgs"]["umap_type"]})
    embed_table_imgs(buffer_dict=src_dict, tree=html, pos="h4", class_name="compare")

    tree = etree.ElementTree(html)
    os.makedirs(save_path, exist_ok=True)
    tree.write(osp.join(save_path, "BatchQC_report_adjust.html"))
    print(
        f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())} The 'BatchQC_report_adjust.html' has been saved to {save_path}")
