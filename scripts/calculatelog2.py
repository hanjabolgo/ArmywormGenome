#!/usr/bin/env python
# -*- coding: utf-8 -*-
import math
import csv
import numpy as np
import pandas as pd
file = 'M2F_cov_raw.txt'
M2F_cov_df = pd.read_csv(file,sep='\t',names=["Chr_id","Chr_sta","Chr_end","Chr_idM","Chr_staM","Chr_endM","Chr_depM","Chr_idF","Chr_staF","Chr_endF","Chr_depF"])
M2F_med = M2F_cov_df['Chr_depM'].median()
F2F_med = M2F_cov_df['Chr_depF'].median()
M2F_cov_df['Chr_depM_nor'] = M2F_cov_df.apply(lambda x: x["Chr_depM"]/M2F_med, axis=1)
M2F_cov_df['Chr_depF_nor'] = M2F_cov_df.apply(lambda x: x["Chr_depF"]/F2F_med, axis=1)
M2F_cov_df['Chr_M2FnorLog2'] = M2F_cov_df.apply(lambda x: math.log2(x['Chr_depM_nor']/(x['Chr_depF_nor'])), axis=1)
M2F_cov_df_subset = M2F_cov_df[["Chr_id","Chr_sta","Chr_end","Chr_M2FnorLog2"]]
M2F_cov_df_subset.to_csv("M2F_cov_nor_log2.txt",sep = "\t")
