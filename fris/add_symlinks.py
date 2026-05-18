#!/usr/bin/env python
import os
import subprocess


def main():
    start_year = 1
    end_year = 8
    #src_path = '/lustre/scratch4/turquoise/vankova/E3SM/scratch/chicoma-cpu/20230404.GMPAS-JRA1p5-DIB-ISMF-TMIX.TL319_ECwISC30to60E2r1.tide_tbl_jb230214.chicoma-cpu/run'
    #src_prefix = '20230404.GMPAS-JRA1p5-DIB-ISMF-TMIX.TL319_ECwISC30to60E2r1.tide_tbl_jb230214.chicoma-cpu'
    #dst_prefix = '20230413.GMPAS-JRA1p5-DIB-ISMF-TMIX.TL319_ECwISC30to60E2r1.tide_tbl_jb230214_branch01.chicoma-cpu'

    src_path = '/lustre/scratch4/turquoise/vankova/E3SM/scratch/chicoma-cpu/20230921.GMPAS-JRA1p5-DIB-PISMF.TL319_ECwISC30to60E2r1.ref.chicoma-cpu/run'
    src_prefix = '20230921.GMPAS-JRA1p5-DIB-PISMF.TL319_ECwISC30to60E2r1.ref.chicoma-cpu'
    dst_prefix = '20230922.GMPAS-JRA1p5-DIB-PISMF.TL319_ECwISC30to60E2r1.b9_sisc03.chicoma-cpu'

    streams = ['mpaso.hist.am.timeSeriesStatsMonthly',
               'mpaso.hist.am.timeSeriesStatsMonthlyMax',
               'mpaso.hist.am.timeSeriesStatsMonthlyMin',
               'mpassi.hist.am.timeSeriesStatsMonthly']

    for stream in streams:
        for year in range(start_year, end_year + 1):
            for month in range(1, 13):
                suffix = f'{stream}.{year:04d}-{month:02}-01.nc'
                src = f'{src_path}/{src_prefix}.{suffix}'
                dst = f'{dst_prefix}.{suffix}'

                args = ['ln', '-sfn', src, dst]
                print(' '.join(args))
                subprocess.run(args, check=True)



if __name__ == '__main__':
    main()
