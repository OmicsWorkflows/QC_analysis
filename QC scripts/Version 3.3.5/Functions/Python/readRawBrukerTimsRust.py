import timsrust_pyo3
import opentims_bruker_bridge
from opentimspy.opentims import OpenTIMS
import pathlib
import pandas as pd
import datetime
import sys, os, ast

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
        date_list = datetime.datetime.fromisoformat(date_raw)
        date_final = date_list.strftime('%Y-%m-%d %H:%M:%S')
        
        return(date_final)
    
    except ModuleNotFoundError:
        print('opentims_bruker_bridge not found')


def getQCpeptide(f, all_frames, m):
    reader = timsrust_pyo3.TimsReader(f)
    
    rt = []
    int = []

    for x in range(0, len(all_frames)):
        tmp = all_frames[x]
        if tmp.frame_type == 0:
            rt.append(tmp.rt/60)

            mzs = reader.resolve_mzs(tmp.tof_indices)

            try:
                index = [i for i, j in enumerate(mzs) if (j > m-0.015) and (j < m+0.015)]
                selected_ints = [tmp.intensities[k] for k in index]
                int.append(sum(selected_ints))
            except ValueError:
                int.append(0)
    
    return([rt, int])

def readQCpeptides(file, matrix, frames, output):
    read = pd.read_csv(matrix)
    df = read[['id', 'precursor.mz']]
    df = df.drop_duplicates()

    for index, row in df.iterrows():
        mz = row['precursor.mz']
        id = row['id']
        individual_qc = getQCpeptide(file, frames, mz)
        
        output_df = pd.DataFrame(
            {'QC': id,
            'rt': individual_qc[0],
            'intensity': individual_qc[1],
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

def readRawBrukerTims(input, file, msn, qc_path, output, qc_output, qc_only):
    path = input+file
    time = getFileTime(path)
    sample = file.replace('.d', '')
    fn = output+sample+'.csv'
    
all_frames = timsrust_pyo3.read_all_frames(path)

    if not qc_only:
        if msn == 'MS':
            ft = [0]
        elif msn == 'MS2':
            ft = [1, 2]

        rt = []
        intensity = []

        for i in range(0, len(all_frames)):
            tmp = all_frames[i]
            if tmp.frame_type in ft:
                rt.append(tmp.rt/60)
                intensity.append(sum(tmp.intensities))

        writeTicCsv(rt, intensity, time, fn)

    if qc_path != '':
        qc_fn = qc_output+sample+'.csv'
        readQCpeptides(path, qc_path, all_frames, qc_fn)

if __name__ == '__main__':
    input = sys.argv[1]
    file = sys.argv[2]
    msn = sys.argv[3]
    qc_path = sys.argv[4]
    output = sys.argv[5]
    qc_output = sys.argv[6]
    qc_only = ast.literal_eval(sys.argv[7])

    try:
        readRawBrukerTims(input, file, msn, qc_path, output, qc_output, qc_only)
        print('ok')
    except:
        print('error')

path = 'C:/Users/435328/Documents/QC workflow/data/6125_BRIMS3_C20_3_1_S6-C3_1_20993.d'
matrix = 'C:/Users/435328/Documents/QC workflow/ref_matrix.csv'
before = datetime.datetime.now()
readQCpeptides(path, matrix, all_frames, 'C:/Users/435328/Documents/QC workflow/test_timsrust.csv')
after = datetime.datetime.now()
