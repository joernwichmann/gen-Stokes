import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams['patch.edgecolor'] = 'none'
import seaborn as sns
import os
from tools import read_datafile, organize_output


### select the experiments whose data will be visualised 
from configs import p_variation as cf
 
if __name__=="__main__":

    ### create output directory
    if not os.path.isdir(cf.OUTPUT_LOCATION):
            os.makedirs(cf.OUTPUT_LOCATION)
    print(f"Start plot of histograms in dataformat '.{cf.HIST_FILEFORMAT}' with dpi '{cf.HIST_DPI}'")
    all_data = {noise: [] for noise in cf.NOISE_TYPES}

    for noise in cf.NOISE_TYPES:
        print(f'\nLoading data for {noise}') 
        for n in range(cf.NUMBER_SAMPLES):
            file_location = cf.ROOT_LOCATION + noise + cf.IND_LOCATION + f"_{n}/" + cf.DATA_SOURCE
            complete_data = read_datafile(file_location)
            time, L1, L2, Linf, SD = organize_output(complete_data)
            for t, l2 in zip(time,L2):
                if t >= cf.STATIONARY_TIME[noise]:
                    all_data[noise].append(l2)
        
        print(f'Plotting data for {noise}')
        if cf.LINEAR_PLOT:
            print(f'\tGenerating linear plot of histogram') 
            plt.figure()
            sns.histplot(all_data[noise], bins="auto", stat="probability", kde=True, color=cf.COLOURS_MEAN[noise], log_scale=(cf.HIST_XAXIS_LOG,False))
            plt.xlabel("Kinetic energy",fontsize=cf.LABEL_FONTSIZE)
            plt.ylabel("Probability",fontsize=cf.LABEL_FONTSIZE)
            #plt.vlines(x=cf.STATIONARY_ENERGY[noise],ymin=0,ymax=cf.YMAX[noise],colors="black",linestyles="solid")
            plt.yticks(fontsize=cf.TICK_FONTSIZE)
            plt.xticks(fontsize=cf.TICK_FONTSIZE)
            plt.tight_layout()
            plt.savefig(cf.OUTPUT_LOCATION+f"hist-{noise}.{cf.HIST_FILEFORMAT}",dpi=cf.HIST_DPI)
            plt.close()
            print(f"Plot saved in '{cf.OUTPUT_LOCATION}hist-{noise}.{cf.HIST_FILEFORMAT}.{cf.TRAJ_FILEFORMAT}'")        
        
        if cf.LOG_PLOT:
            print(f'\tGenerating log plot of histogram') 
            plt.figure()
            sns.histplot(all_data[noise], bins="auto", stat="probability", kde=True, color=cf.COLOURS_MEAN[noise], log_scale=(True,True))
            plt.xlabel("Kinetic energy")
            #plt.vlines(x=cf.STATIONARY_ENERGY[noise],ymin=0,ymax=cf.YMAX[noise],colors="black",linestyles="solid")
            plt.tight_layout()
            plt.savefig(cf.OUTPUT_LOCATION+f"hist-{noise}-log.{cf.HIST_FILEFORMAT}",dpi=cf.HIST_DPI)
            plt.close()
            print(f"Plot saved in '{cf.OUTPUT_LOCATION}hist-{noise}-log.{cf.HIST_FILEFORMAT}.{cf.TRAJ_FILEFORMAT}'")
    
        

    print(f'\nPlotting data for noise comparisons') 
    ##legend
    # Create dummy Line2D objects for legend
    legendMarkers = [Line2D([0], [0.1], color=cf.COLOURS_MEAN[noise], linewidth = cf.LINEWIDTH_MEAN) for noise in cf.NOISE_TYPES]
    #legendMarkers.append(Line2D([0], [0.1], color=cf.BLACK, linewidth = cf.LINEWIDTH_MEAN))
    legendPvalues = [f"p = {cf.P_VALUE[noise]}" for noise in cf.NOISE_TYPES]
    #legendPvalues.append("det")
    if cf.LINEAR_PLOT:
        print(f'\tGenerating linear plot of histogram')
        plt.figure()
        for noise in cf.NOISE_TYPES:
            sns.histplot(all_data[noise], bins="auto", stat="probability", kde=True, color=cf.COLOURS_MEAN[noise], log_scale=(cf.HIST_XAXIS_LOG,False))
            #plt.vlines(x=cf.STATIONARY_ENERGY[noise],ymin=0,ymax=cf.YMAX[noise],colors="black",linestyles="solid")
        plt.xlabel("Kinetic energy")
        plt.legend(legendMarkers,legendPvalues)
        plt.tight_layout()
        plt.savefig(cf.OUTPUT_LOCATION+f"hist-{cf.EXPERIMENT_NAME}-all.{cf.HIST_FILEFORMAT}",dpi=cf.HIST_DPI)
        plt.close()
        print(f"Plot saved in '{cf.OUTPUT_LOCATION}hist-{cf.EXPERIMENT_NAME}-all.{cf.HIST_FILEFORMAT}'")

    if cf.LOG_PLOT:
        print(f'\tGenerating log plot of histogram')
        plt.figure()
        for noise in cf.NOISE_TYPES:
            sns.histplot(all_data[noise], bins="auto", stat="probability", kde=True, color=cf.COLOURS_MEAN[noise], log_scale=(True,True))
            #plt.vlines(x=cf.STATIONARY_ENERGY[noise],ymin=0,ymax=cf.YMAX[noise],colors="black",linestyles="solid")
        plt.xlabel("Kinetic energy")
        plt.legend(legendMarkers,legendPvalues)
        plt.tight_layout()
        plt.savefig(cf.OUTPUT_LOCATION+f"hist-{cf.EXPERIMENT_NAME}-all-log.{cf.HIST_FILEFORMAT}",dpi=cf.HIST_DPI)
        plt.close()
        print(f"Plot saved in '{cf.OUTPUT_LOCATION}hist-{cf.EXPERIMENT_NAME}-all-log.{cf.HIST_FILEFORMAT}'")