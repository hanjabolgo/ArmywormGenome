##!/usr/bin/env python
# -*- coding: utf-8 -*-
# Auther： amane
# URL： https://www.zhouxiaozhao.cn/2020/11/26/2020-11-26-annotion(4)/

import sys
import re
fout = open(sys.argv[3],'w')

ref_dict={}
with open(sys.argv[1]) as gene_a:
        for line in gene_a:
                line_s = line.strip().split('\t')
                info = re.split('=|;',line_s[8])
                ID = info[1]
                ref_set = []
                for n in range(0,8):
                        ref_set.append(line_s[n])
                ref_dict.setdefault(ID,[]).append(ref_set)
                ref_dict.setdefault(ID,[]).append(info[3])

with open(sys.argv[2]) as mrna:
        for eachline in mrna:
                i = eachline.strip().split('\t')
                info1 = re.split('=|;',i[8])
                parent = info1[3]
                mrna_n = info1[1]
                ref_set1 = []
                for a in range(0,8):
                        ref_set1.append(i[a])
                mrna_h = '\t'.join(ref_set1)
                if parent in ref_dict:
                        vs = ref_dict[parent]
                        head = '\t'.join(vs[0])
                        fout.write('%s\tID=%s;Name=%s\n'%(head,parent,vs[1]))
                        fout.write('%s\tID=%s;Parent=%s;Name=%s\n'%(mrna_h,mrna_n,parent,vs[1]))

fout.close()
