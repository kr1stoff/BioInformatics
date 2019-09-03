# -*- coding:utf-8 -*-

import numpy as np 
import pandas as pd


def dfep_filterbydiff(infile, outfile):
    dfep = pd.read_excel(infile, encoding='utf-8')
    dfep_filt = dfep[(dfep['HWE_p.value'] > 0.05) & (dfep['allelic_chisq_p.value'] < 0.05) & \
(dfep['het_case_count'] >= dfep['het_control_count']) & (dfep['althom_case_count'] >= dfep['althom_control_count'])]
    dfep_filt = dfep_filt.fillna('-')
    dfep_filt['het_diff'] = dfep_filt['het_case_count'] - dfep_filt['het_control_count']
    dfep_filt['hom_diff'] = dfep_filt['althom_case_count'] - dfep_filt['althom_control_count']

    dfep_filt.sort_values(by=['het_diff', 'hom_diff'], ascending=False, inplace=True)
    dfep_filt_cols = list(dfep_filt.columns)
    dfep_filt = pd.DataFrame(dfep_filt, columns=(dfep_filt_cols[:2]+dfep_filt_cols[-2:]+dfep_filt_cols[2:-2]))
    dfep_filt.to_excel(outfile, encoding='utf-8', index=False)


if __name__ == '__main__':
    infile = 'D:/xiangfu.meng/2019/癫痫科研/癫痫研发数据/output/cntr-ep差异分析结果_filtbygene.xlsx'
    outfile = 'D:/xiangfu.meng/2019/1909/0902/dfep_filtbydiff.xlsx'
    dfep_filterbydiff(infile, outfile)