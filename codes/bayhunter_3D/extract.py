from glob import glob
import os
import numpy as np

for i in glob('./RESULTS/*'):

    lon = os.path.basename(i).split('_')[4]
    lat = os.path.basename(i).split('_')[5][:-4]
    print(i, lat, lon)
    try:
        depth = np.loadtxt(f'{i}/model.txt', usecols=0, skiprows=1, delimiter=',')
        vs = np.loadtxt(f'{i}/model.txt', usecols=1, skiprows=1, delimiter=',')
    
        np.savetxt(f"model_{lat}_{lon}_.txt", np.column_stack((depth, vs)), delimiter=',')
    except Exception as e:
        print(e)
