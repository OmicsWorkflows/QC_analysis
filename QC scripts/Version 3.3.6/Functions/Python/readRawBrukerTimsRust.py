import timsrust_pyo3
import opentims_bruker_bridge
from opentimspy.opentims import OpenTIMS
import pathlib
import pandas as pd
from datetime import datetime
import sys, os, ast
import psutil
import time
import numpy as np

def writeSysInfo(file, process):
    cpu = psutil.cpu_percent()
    mem = psutil.virtual_memory()
    df = pd.DataFrame(
        {'process': process,
         'CPU.perc': cpu,
         'RAM.used': mem.used,
         'RAM.free': mem.free,
         'RAM.perc': mem.percent,
         'time': datetime.today().isoformat()
        }, index=[0])
    if os.path.isfile(file):
        df.to_csv(file, mode = 'a', header = False, index = False)
    else:
        df.to_csv(file, index = False)

def writeTicCsv(rt, intensity, time, filename):
    df = pd.DataFrame(
        {'RT': rt,
         'intensity': intensity,
         'time': time
        })
    df.to_csv(filename, index = False)

def getFileTime(f):
    f_path = pathlib.Path(f)
    D = OpenTIMS(f_path)

    try:
        import opentims_bruker_bridge
        
        metadict = D.table2dict('GlobalMetadata')
        k = metadict['Key'].tolist()
        v = metadict['Value'].tolist()
        date_raw = v[k.index('AcquisitionDateTime')]
        date_list = datetime.fromisoformat(date_raw)
        date_final = date_list.strftime('%Y-%m-%d %H:%M:%S')
        
        return(date_final)
    
    except ModuleNotFoundError:
        print('opentims_bruker_bridge not found')


def getQCpeptide(file, frames, m_list, mz_tol, info_path):
    rt = []
    int = [np.zeros(len(frames)) for _ in range(len(m_list))]

    metadata = timsrust_pyo3.Metadata(file + "/analysis.tdf")

    #writeSysInfo(info_path, process='started readQCpeptide fcn')

    for x, tmp in enumerate(frames):
        #writeSysInfo(info_path, process='calculating QC XIC')
        rt.append(tmp.rt / 60)
        mzs = np.array(metadata.resolve_mzs(tmp.tof_indices))

        for m_i, tmp_mz in enumerate(m_list):
            index = np.where((mzs > tmp_mz - mz_tol) & (mzs < tmp_mz + mz_tol))[0]
            
            if len(index) > 0:
                selected_ints = np.take(tmp.intensities, index)
                int[m_i][x] = np.sum(selected_ints)
            else:
                int[m_i][x] = 0
    return([rt, int])


def readQCpeptides(file, matrix, frames, mz_tol, output, info_path):
    #writeSysInfo(info_path, process='started readQCpeptides fcn')

    read = pd.read_csv(matrix)
    df = read[['id', 'precursor.mz']]
    df = df.drop_duplicates()

    mz_list = df['precursor.mz'].tolist()

    list_qc = getQCpeptide(file, frames, mz_list, mz_tol, info_path)

    for i in range(0, len(mz_list)):   
        output_df = pd.DataFrame(
            {'QC': df.at[i, 'id'],
            'rt': list_qc[0],
            'intensity': list_qc[1][i],
            'matrix_correlation': 'NA',
            'peak_start': 'NA',
            'peak_end': 'NA',
            'peak_x': 'NA',
            'baseline': 'NA'
            })

        if os.path.isfile(output):
            output_df.to_csv(output, mode = 'a', header = False, index = False)
        else:
            output_df.to_csv(output, index = False)


def readRawBrukerTims(input, file, msn,  qc_path, mz_tol, tic_output, bpc_output, qc_output, qc_only, info_path):
    #writeSysInfo(info_path, process='started readRawBrukerTims fcn')
    path = input+file
    time = getFileTime(path)
    #writeSysInfo(info_path, process='read file time')
    sample = file.replace('.d', '')

    reader = timsrust_pyo3.FrameReader(path)

    if msn == 'MS':
        all_frames = reader.read_ms1_frames()
    elif msn == 'MS2':
        all_frames = reader.read_ms2_frames()
    #writeSysInfo(info_path, process='read frames - ' + msn)

    if not qc_only:
        rt = []
        tic = []

        if bpc_output != '':
            bpc = []

        for i in range(0, len(all_frames)):
            tmp = all_frames[i]
            rt.append(tmp.rt/60)
            tic.append(sum(tmp.intensities))
            if bpc_output != '':
                bpc.append(max(tmp.intensities))
            #writeSysInfo(info_path, process='calculating TIC - ' + msn)

        writeTicCsv(rt, tic, time, tic_output+sample+'.csv')
        if bpc_output != '':
            writeTicCsv(rt, bpc, time, bpc_output+sample+'.csv')

    if qc_path != '':
        qc_fn = qc_output+sample+'.csv'
        readQCpeptides(path, qc_path, all_frames, mz_tol, qc_fn, info_path)


if __name__ == '__main__':
    input = sys.argv[1]
    file = sys.argv[2]
    msn = sys.argv[3]
    qc_path = sys.argv[4]
    mz_tol = float(sys.argv[5])
    tic_output = sys.argv[6]
    bpc_output = sys.argv[7]
    qc_output = sys.argv[8]
    qc_only = ast.literal_eval(sys.argv[9])
    info_path = tic_output + 'info.csv'

    #writeSysInfo(info_path, process='start')

    try:
        readRawBrukerTims(input, file, msn, qc_path, mz_tol, tic_output, bpc_output, qc_output, qc_only, info_path)
        print('ok')
    except:
        print('error')

input = "C:/Users/435328/Documents/QC workflow/data/"
file =  "6316_BRIMS1_DIA_108min_TFA_G18_2_1_GA7_1_24176.d"
msn = 'MS'
qc_path = "C:/Users/435328/Documents/QC workflow/ref_matrix.csv"
mz_tol = float("0.015")
tic_output = "C:/Users/435328/Documents/QC workflow/outputs/tic_data/MS/"
bpc_output =  "C:/Users/435328/Documents/QC workflow/outputs/bpc_data/"
qc_output = "C:/Users/435328/Documents/QC workflow/outputs/QC_data/"
qc_only = ast.literal_eval("False")
info_path = tic_output + 'info.csv'