import timsrust_pyo3
import opentims_bruker_bridge
from opentimspy.opentims import OpenTIMS
import pathlib
import pandas as pd
from datetime import datetime
import sys, os, ast
import psutil
import time

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


def getQCpeptide(file, frames, m_list, info_path):
    rt = []
    int = [[] for _ in range(0, len(m_list))]

    metadata = timsrust_pyo3.Metadata(file + "/analysis.tdf")

    writeSysInfo(info_path, process='started readQCpeptide fcn')

    for x in range(0, len(frames)):
        writeSysInfo(info_path, process='calculating QC XIC')
        tmp = frames[x]
        rt.append(tmp.rt/60)
        mzs = metadata.resolve_mzs(tmp.tof_indices)
        for m_i in range(0, len(m_list)):
            tmp_mz = m_list[m_i]
            try:
                index = [i for i, j in enumerate(mzs) if (j > tmp_mz-0.015) and (j < tmp_mz+0.015)]
                selected_ints = [tmp.intensities[k] for k in index]
                int[m_i].append((sum(selected_ints)))
            except ValueError:
                int[m_i].append(0)
    return([rt, int])

def readQCpeptides(file, matrix, frames, output, info_path):
    writeSysInfo(info_path, process='started readQCpeptides fcn')
    read = pd.read_csv(matrix)
    df = read[['id', 'precursor.mz']]
    df = df.drop_duplicates()

    mz_list = df['precursor.mz'].tolist()

    list_qc = getQCpeptide(file, frames, mz_list, info_path)

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


def readRawBrukerTims(input, file, msn, qc_path, output, qc_output, qc_only, info_path):
    writeSysInfo(info_path, process='started readRawBrukerTims fcn')
    path = input+file
    time = getFileTime(path)
    writeSysInfo(info_path, process='read file time')
    sample = file.replace('.d', '')
    fn = output+sample+'.csv'

    reader = timsrust_pyo3.FrameReader(path)

    if msn == 'MS':
        all_frames = reader.read_ms1_frames()
    elif msn == 'MS2':
        all_frames = reader.read_ms2_frames()
    writeSysInfo(info_path, process='read frames - ' + msn)

    if not qc_only:
        rt = []
        intensity = []

        for i in range(0, len(all_frames)):
            tmp = all_frames[i]
            rt.append(tmp.rt/60)
            intensity.append(sum(tmp.intensities))
            writeSysInfo(info_path, process='calculating TIC - ' + msn)

        writeTicCsv(rt, intensity, time, fn)

    if qc_path != '':
        qc_fn = qc_output+sample+'.csv'
        readQCpeptides(path, qc_path, all_frames, qc_fn, info_path)

if __name__ == '__main__':
    input = sys.argv[1]
    file = sys.argv[2]
    msn = sys.argv[3]
    qc_path = sys.argv[4]
    output = sys.argv[5]
    qc_output = sys.argv[6]
    qc_only = ast.literal_eval(sys.argv[7])
    info_path = output + 'info.csv'

    writeSysInfo(info_path, process='start')

    try:
        readRawBrukerTims(input, file, msn, qc_path, output, qc_output, qc_only, info_path)
        print('ok')
    except:
        print('error')

input = 'C:/Users/435328/Documents/QC workflow/data/'
file = '6087_BRIMS3_DIA2_RSLC_2col-108min_TFA_56_2_BE2_1_18895.d'
sys_info_output = 'C:/Users/435328/Documents/QC workflow/timsRustInfoMs2Spectra.csv'
path = 'C:/Users/435328/Documents/QC workflow/data/6087_BRIMS3_DIA2_RSLC_2col-108min_TFA_56_2_BE2_1_18895.d'
qc_path = 'C:/Users/435328/Documents/QC workflow/ref_matrix.csv'
msn = 'MS'
output = 'C:/Users/435328/Documents/QC workflow/'
qc_output = output
qc_only = False
info_path = output + 'info.csv'
reader = timsrust_pyo3.FrameReader(path)
all_frames = reader.read_ms1_frames()

read = pd.read_csv(qc_path)
df = read[['id', 'precursor.mz']]
df = df.drop_duplicates()

readRawBrukerTims(input, file, msn, qc_path, output, qc_output, qc_only, info_path)