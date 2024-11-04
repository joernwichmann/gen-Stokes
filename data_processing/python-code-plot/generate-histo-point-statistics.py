import matplotlib.pyplot as plt
plt.rcParams['patch.edgecolor'] = 'none'
import seaborn as sns
import numpy as np
import pandas as pd
from tools_point_statistics import read_datafile, organize_output_single
from configs import lid_driven_point_statistics as cf

 
if __name__=="__main__":
    print(f"Start plot of histograms in dataformat '.{cf.HIST_FILEFORMAT}' with dpi '{cf.HIST_DPI}'")
    all_data_x = dict()
    all_data_y = dict()
    for expID in cf.EXPERIMENTS.keys():
        data_x = []
        data_y = []
        print(f'\nLoading data for {cf.EXPERIMENTS[expID]}') 
        for n in range(cf.NUMBER_SAMPLES):
            file_location = cf.ROOT_LOCATION + cf.EXPERIMENTS[expID] + cf.IND_LOCATION + f"_{n}" + cf.DATA_SOURCE
            complete_data = read_datafile(file_location)
            time, val_0, val_1 = organize_output_single(complete_data)
            for t, v0, v1 in zip(time,val_0,val_1):
                if t >= cf.STATIONARY_TIME[expID]:
                    data_x.append(v0)
                    data_y.append(v1)

        all_data_x[expID] = np.array(data_x)
        all_data_y[expID] = np.array(data_y)
        
        print(f'Plotting data for {cf.EXPERIMENTS[expID]}') 
        
        plt.figure()
        #sns.histplot(all_data[expID], bins="auto", stat="probability", kde=True, color=cf.COLOURS_MEAN[expID], log_scale=(True,False))
        sns.jointplot(x = all_data_x[expID], y=all_data_y[expID], color=cf.COLOURS_MEAN[expID], kind="kde")
        #plt.xlabel("Kinetic energy",fontsize=cf.LABEL_FONTSIZE)
        #plt.ylabel("Probability",fontsize=cf.LABEL_FONTSIZE)
        #plt.vlines(x=cf.STATIONARY_ENERGY,ymin=0,ymax=cf.YMAX[expID],colors="black",linestyles="solid")
        #plt.yticks(fontsize=cf.TICK_FONTSIZE)
        #plt.xticks(fontsize=cf.TICK_FONTSIZE)
        plt.plot(cf.STATIONARY_VAL_X[expID],cf.STATIONARY_VAL_Y[expID], marker = "o",color=cf.COLOURS_MEAN[expID])
        plt.tight_layout()
        plt.savefig(f"hist-point-{cf.EXPERIMENTS[expID]}.{cf.HIST_FILEFORMAT}",dpi=cf.HIST_DPI)
        plt.close()        
    

    #build pandas.dataframe
    prebuild = []
    for expID in cf.EXPERIMENTS.keys():
        for x,y in zip(all_data_x[expID],all_data_y[expID]):
            prebuild.append([f"p = {cf.P_VALUE[expID]}",x,y])
    build_data = pd.DataFrame(data=prebuild, columns=["experiment","x","y"])

    #build color palette
    colors = [cols for cols in cf.COLOURS_MEAN.values()]
    customPalette = sns.set_palette(sns.color_palette(colors))

    print(f'\nPlotting data for noise comparisons') 
    plt.figure()
    sns.jointplot(data = build_data, x = "x", y="y",  hue="experiment", kind="kde", palette=customPalette)
    for expID in cf.EXPERIMENTS.keys():
        plt.plot(cf.STATIONARY_VAL_X[expID],cf.STATIONARY_VAL_Y[expID], marker = "o", color=cf.COLOURS_MEAN[expID])
    plt.tight_layout()
    plt.savefig(f"hist-point-all.{cf.HIST_FILEFORMAT}",dpi=cf.HIST_DPI)
    plt.close()