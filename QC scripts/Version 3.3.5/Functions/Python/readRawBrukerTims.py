import pathlib
import opentims_bruker_bridge
from opentimspy.opentims import OpenTIMS
import pandas as pd
import datetime
import sys

def writeTicCsv(timsobj, frames, intensityAll, time, filename):
    rt = timsobj.frame2retention_time(frames).tolist()
    rt_min = [r/60 for r in rt]
    index = [(x-1) for x in frames]
    intensity = [intensityAll[i] for i in index]
    df = pd.DataFrame(
        {'RT': rt_min,
         'intensity': intensity,
         'time': time
        })
    df.to_csv(filename, index = False)

def readRawBrukerTims(input, file, msn, output):
    path = pathlib.Path(input+file)
    D = OpenTIMS(path)

    sample = file.replace('.d', '')

    try:
        import opentims_bruker_bridge

        # get date and time of .d file creation
        metadict = D.table2dict('GlobalMetadata')
        k = metadict['Key'].tolist()
        v = metadict['Value'].tolist()
        date_raw = v[k.index('AcquisitionDateTime')]
        date_list = datetime.datetime.fromisoformat(date_raw)
        date_final = date_list.strftime('%Y-%m-%d %H:%M:%S')

        intensityAll = D.framesTIC().tolist()
        fn = output+sample+'.csv'

        if msn == 'MS':
            frames = [x for x in D.ms1_frames]
            writeTicCsv(D, frames, intensityAll, date_final, fn)
        
        if msn == 'MS2':
            frames = [x for x in D.ms2_frames]
            writeTicCsv(D, frames, intensityAll, date_final, fn)

    except ModuleNotFoundError:
        print('opentims_bruker_bridge not found')

if __name__ == '__main__':
    input = sys.argv[1]
    file = sys.argv[2]
    msn = sys.argv[3]
    output = sys.argv[4]
    
    try:
        readRawBrukerTims(input, file, msn, output)
        print('ok')
    except:
        print('error')