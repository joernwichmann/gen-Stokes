import matplotlib.pyplot as plt
from tools import read_datafile, organize_output
import configs as cf
from tqdm import tqdm

#TODO: TIKZ FONT size
if __name__=="__main__":
    print(f"Start plot of trajectories in dataformat '.{cf.TRAJ_FILEFORMAT}' with dpi '{cf.TRAJ_DPI}'")
    
    ################################################################### 
    ################## comparison between noise types
    ################################################################### 
    print(f'\nPlotting comparsion of different noise types')
    plt.figure()

    #plot individual energy trajectory
    print(f"\tPlotting individual trajectories")
    for n in tqdm(range(cf.NUMBER_SAMPLES)):
        for noise in cf.NOISE_TYPES:
            file_location = cf.ROOT_LOCATION + noise + cf.IND_LOCATION + f"_{n}/" + cf.DATA_SOURCE
            complete_data = read_datafile(file_location)
            time, L1, L2, Linf, SD = organize_output(complete_data)
            plt.plot(time,L2,color = cf.COLOURS_INDIVIDUAL[noise],linewidth=cf.LINEWIDTH_INDIVIDUAL,alpha=cf.LINEOPACITY_INDIVIDUAL)

    #plot mean and standard deviations
    print(f"\tPlotting mean and standard devations")
    for noise in cf.NOISE_TYPES:
        file_location = cf.ROOT_LOCATION + noise + cf.MEAN_LOCATION + cf.DATA_SOURCE
        complete_data = read_datafile(file_location)
        time, L1, L2, Linf, SD = organize_output(complete_data)
        conf_plus = []
        conf_minus = []
        for l2, sd in zip(L2, SD):
            conf_plus.append(l2 + sd)
            conf_minus.append(l2-sd) 
        #mean energy trajectory
        plt.plot(time,L2,color = cf.COLOURS_MEAN[noise],linewidth=cf.LINEWIDTH_MEAN,alpha=cf.LINEOPACITY_MEAN)
        #standard deviation
        plt.plot(time,conf_plus,color = cf.COLOURS_MEAN[noise],linewidth=cf.LINEWIDTH_SD,alpha=cf.LINEOPACITY_SD,linestyle=cf.LINESTYLE_SD)
        plt.plot(time,conf_minus,color = cf.COLOURS_MEAN[noise],linewidth=cf.LINEWIDTH_SD,alpha=cf.LINEOPACITY_SD,linestyle=cf.LINESTYLE_SD)

    #print(f"\tPlotting deterministic energy")
    #file_location = cf.ROOT_LOCATION + "intro-deterministic" + cf.MEAN_LOCATION + cf.DATA_SOURCE
    #complete_data = read_datafile(file_location)
    #time, L1, L2, Linf, SD = organize_output(complete_data)
    #plt.plot(time,L2,color = cf.BLACK,linewidth=1,alpha=cf.LINEOPACITY_MEAN)
    
    # styling the plot
    plt.ylabel("Kinetic energy")
    plt.xlabel("Time")
    plt.yscale("log")
    plt.xlim(0.0,1.0)
    plt.tight_layout()
    plt.savefig(f"intro-all-trajectories.{cf.TRAJ_FILEFORMAT}",dpi=cf.TRAJ_DPI)
    plt.close()
    print(f"Plot saved in 'intro-all-trajectories.{cf.TRAJ_FILEFORMAT}'")

    ################################################################### 
    ############# Plot of energies for fixed noise type
    ###################################################################     
    # #plot individual energy trajectory
    for noise in cf.NOISE_TYPES:
        print(f'\nGenerating trajectory plot for energies of {noise}')
        plt.figure()

        print(f"\tPlotting individual trajectories")
        for n in tqdm(range(cf.NUMBER_SAMPLES)):
            file_location = cf.ROOT_LOCATION + noise + cf.IND_LOCATION + f"_{n}/" + cf.DATA_SOURCE
            complete_data = read_datafile(file_location)
            time, L1, L2, Linf, SD = organize_output(complete_data)
            plt.plot(time,L2,color = cf.COLOURS_INDIVIDUAL[noise],linewidth=cf.LINEWIDTH_INDIVIDUAL,alpha=cf.LINEOPACITY_INDIVIDUAL)

        #plot mean and standard deviations
        print(f"\tPlotting mean and standard devations")
        file_location = cf.ROOT_LOCATION + noise + cf.MEAN_LOCATION + cf.DATA_SOURCE
        complete_data = read_datafile(file_location)
        time, L1, L2, Linf, SD = organize_output(complete_data)
        conf_plus = []
        conf_minus = []
        for l2, sd in zip(L2, SD):
            conf_plus.append(l2 + sd)
            conf_minus.append(l2-sd) 
        #mean energy trajectory
        plt.plot(time,L2,color = cf.COLOURS_MEAN[noise],linewidth=cf.LINEWIDTH_MEAN,alpha=cf.LINEOPACITY_MEAN)
        #standard deviation
        plt.plot(time,conf_plus,color = cf.COLOURS_MEAN[noise],linewidth=cf.LINEWIDTH_SD,alpha=cf.LINEOPACITY_SD,linestyle=cf.LINESTYLE_SD)
        plt.plot(time,conf_minus,color = cf.COLOURS_MEAN[noise],linewidth=cf.LINEWIDTH_SD,alpha=cf.LINEOPACITY_SD,linestyle=cf.LINESTYLE_SD)

        #print(f"\tPlotting deterministic energy")
        #file_location = cf.ROOT_LOCATION + "intro-deterministic" + cf.MEAN_LOCATION + cf.DATA_SOURCE
        #complete_data = read_datafile(file_location)
        #time, L1, L2, Linf, SD = organize_output(complete_data)
        #plt.plot(time,L2,color = cf.BLACK,linewidth=1,alpha=cf.LINEOPACITY_MEAN)
        

        ### styling the plot
        plt.ylabel("Kinetic energy", fontsize=cf.LABEL_FONTSIZE)
        plt.xlabel("Time", fontsize=cf.LABEL_FONTSIZE)
        plt.yscale("log")
        plt.yticks(fontsize=cf.TICK_FONTSIZE)
        plt.xticks(fontsize=cf.TICK_FONTSIZE)
        plt.xlim(0.0,1.0)
        plt.tight_layout()
        plt.savefig(f"{noise}-trajectories.{cf.TRAJ_FILEFORMAT}",dpi=cf.TRAJ_DPI)
        plt.close()
        print(f"Plot saved in '{noise}-trajectories.{cf.TRAJ_FILEFORMAT}'")





