import matplotlib.pyplot as plt
from tools_point_statistics import read_datafile, organize_output, organize_output_single
from configs import lid_driven_point_statistics as cf
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
        for expID in cf.EXPERIMENTS.keys():
            file_location = cf.ROOT_LOCATION + cf.EXPERIMENTS[expID] + cf.IND_LOCATION + f"_{n}" + cf.DATA_SOURCE
            complete_data = read_datafile(file_location)
            time, val_0, val_1 = organize_output_single(complete_data)
            plt.plot(val_0,val_1,color = cf.COLOURS_INDIVIDUAL[expID],linewidth=cf.LINEWIDTH_INDIVIDUAL,alpha=cf.LINEOPACITY_INDIVIDUAL)

    #plot mean and standard deviations
    print(f"\tPlotting mean, standard devations and deterministic")
    for expID in cf.EXPERIMENTS.keys():
        file_location = cf.ROOT_LOCATION + cf.EXPERIMENTS[expID] + cf.MEAN_LOCATION + cf.DATA_SOURCE
        complete_data = read_datafile(file_location)
        time, l1_0, l1_1, l2_0, l2_1, linf_0, linf_1, sd_0, sd_1 = organize_output(complete_data)
        l2_0_plus = []
        l2_0_minus = []
        l2_1_plus = []
        l2_1_minus = []
        for _l2, _sd in zip(l2_0, sd_0):
            l2_0_plus.append(_l2 + _sd)
            l2_0_minus.append(_l2-_sd)
        for _l2, _sd in zip(l2_1, sd_1):
            l2_1_plus.append(_l2 + _sd)
            l2_1_minus.append(_l2-_sd)

        #mean energy trajectory
        plt.plot(l1_0,l1_1,color = cf.COLOURS_MEAN[expID],linewidth=cf.LINEWIDTH_MEAN,alpha=cf.LINEOPACITY_MEAN)
        #standard deviation
        #plt.plot(l2_0_plus,l2_1_plus,color = cf.COLOURS_MEAN[expID],linewidth=cf.LINEWIDTH_SD,alpha=cf.LINEOPACITY_SD,linestyle=cf.LINESTYLE_SD)
        #plt.plot(l2_0_plus,l2_1_minus,color = cf.COLOURS_MEAN[expID],linewidth=cf.LINEWIDTH_SD,alpha=cf.LINEOPACITY_SD,linestyle=cf.LINESTYLE_SD)
        #plt.plot(l2_0_minus,l2_1_plus,color = cf.COLOURS_MEAN[expID],linewidth=cf.LINEWIDTH_SD,alpha=cf.LINEOPACITY_SD,linestyle=cf.LINESTYLE_SD)
        #plt.plot(l2_0_minus,l2_1_minus,color = cf.COLOURS_MEAN[expID],linewidth=cf.LINEWIDTH_SD,alpha=cf.LINEOPACITY_SD,linestyle=cf.LINESTYLE_SD)

        #plot det
        file_location = cf.ROOT_LOCATION + cf.EXPERIMENTS[expID] + cf.DET_LOCATION + cf.DATA_SOURCE
        complete_data = read_datafile(file_location)
        time, l1_0, l1_1, l2_0, l2_1, linf_0, linf_1, sd_0, sd_1 = organize_output(complete_data)
        plt.plot(l1_0,l1_1,color = cf.BLACK,linewidth=1,alpha=cf.LINEOPACITY_MEAN)
        plt.plot(cf.STATIONARY_VAL_X[expID],cf.STATIONARY_VAL_Y[expID], marker = "o",color=cf.COLOURS_MEAN[expID])
        
    
    # styling the plot
    plt.ylabel("y")
    plt.xlabel("x")
    plt.title('velocity @ [0.5,0.75]')
    #plt.xlim(0.0,1.0)
    plt.tight_layout()
    plt.savefig(f"point-{cf.EXPERIMENT_NAME}-all-trajectories.{cf.TRAJ_FILEFORMAT}",dpi=cf.TRAJ_DPI)
    plt.close()
    print(f"Plot saved in 'point-{cf.EXPERIMENT_NAME}-all-trajectories.{cf.TRAJ_FILEFORMAT}'")

    ################################################################### 
    ############# Plot of energies for fixed noise type
    ###################################################################     
    # #plot individual energy trajectory
    for expID in cf.EXPERIMENTS.keys():
        print(f'\nGenerating trajectory plot for Experiment {cf.EXPERIMENTS[expID]}')
        plt.figure()

        print(f"\tPlotting individual trajectories")
        for n in tqdm(range(cf.NUMBER_SAMPLES)):
            file_location = cf.ROOT_LOCATION + cf.EXPERIMENTS[expID] + cf.IND_LOCATION + f"_{n}" + cf.DATA_SOURCE
            complete_data = read_datafile(file_location)
            time, val_0, val_1 = organize_output_single(complete_data)
            plt.plot(val_0,val_1,color = cf.COLOURS_INDIVIDUAL[expID],linewidth=cf.LINEWIDTH_INDIVIDUAL,alpha=cf.LINEOPACITY_INDIVIDUAL)

        #plot mean and standard deviations
        print(f"\tPlotting mean and standard devations")
        file_location = cf.ROOT_LOCATION + cf.EXPERIMENTS[expID] + cf.MEAN_LOCATION + cf.DATA_SOURCE
        complete_data = read_datafile(file_location)
        time, l1_0, l1_1, l2_0, l2_1, linf_0, linf_1, sd_0, sd_1 = organize_output(complete_data)
        l2_0_plus = []
        l2_0_minus = []
        l2_1_plus = []
        l2_1_minus = []
        for _l2, _sd in zip(l2_0, sd_0):
            l2_0_plus.append(_l2 + _sd)
            l2_0_minus.append(_l2-_sd)
        for _l2, _sd in zip(l2_1, sd_1):
            l2_1_plus.append(_l2 + _sd)
            l2_1_minus.append(_l2-_sd)

        #mean energy trajectory
        plt.plot(l1_0,l1_1,color = cf.COLOURS_MEAN[expID],linewidth=cf.LINEWIDTH_MEAN,alpha=cf.LINEOPACITY_MEAN)
        #standard deviation
        #plt.plot(l2_0_plus,l2_1_plus,color = cf.COLOURS_MEAN[expID],linewidth=cf.LINEWIDTH_SD,alpha=cf.LINEOPACITY_SD,linestyle=cf.LINESTYLE_SD)
        #plt.plot(l2_0_plus,l2_1_minus,color = cf.COLOURS_MEAN[expID],linewidth=cf.LINEWIDTH_SD,alpha=cf.LINEOPACITY_SD,linestyle=cf.LINESTYLE_SD)
        #plt.plot(l2_0_minus,l2_1_plus,color = cf.COLOURS_MEAN[expID],linewidth=cf.LINEWIDTH_SD,alpha=cf.LINEOPACITY_SD,linestyle=cf.LINESTYLE_SD)
        #plt.plot(l2_0_minus,l2_1_minus,color = cf.COLOURS_MEAN[expID],linewidth=cf.LINEWIDTH_SD,alpha=cf.LINEOPACITY_SD,linestyle=cf.LINESTYLE_SD)

        print(f"\tPlotting deterministic energy")
        file_location = cf.ROOT_LOCATION + cf.EXPERIMENTS[expID] + cf.DET_LOCATION + cf.DATA_SOURCE
        complete_data = read_datafile(file_location)
        time, l1_0, l1_1, l2_0, l2_1, linf_0, linf_1, sd_0, sd_1 = organize_output(complete_data)
        plt.plot(l1_0,l1_1,color = cf.BLACK,linewidth=1,alpha=cf.LINEOPACITY_MEAN)
        plt.plot(cf.STATIONARY_VAL_X[expID],cf.STATIONARY_VAL_Y[expID], marker = "o",color=cf.COLOURS_MEAN[expID])
        

        ### styling the plot
        plt.ylabel("y", fontsize=cf.LABEL_FONTSIZE)
        plt.xlabel("x", fontsize=cf.LABEL_FONTSIZE)
        plt.title('velocity @ [0.5,0.75]', fontsize=cf.LABEL_FONTSIZE)
        plt.yticks(fontsize=cf.TICK_FONTSIZE)
        plt.xticks(fontsize=cf.TICK_FONTSIZE)
        #plt.xlim(0.0,1.0)
        plt.tight_layout()
        plt.savefig(f"point-{cf.EXPERIMENTS[expID]}-trajectories.{cf.TRAJ_FILEFORMAT}",dpi=cf.TRAJ_DPI)
        plt.close()
        print(f"Plot saved in 'point-{cf.EXPERIMENTS[expID]}-trajectories.{cf.TRAJ_FILEFORMAT}'")





