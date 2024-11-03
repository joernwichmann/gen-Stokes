import numpy as np
import os
import shutil
import csv



def _update_mean(samples: int, old_mean_vector: np.ndarray, update_vector: np.ndarray) -> np.ndarray:
    """Update the mean value."""
    if not old_mean_vector.shape == update_vector.shape:
        msg = "Shape of mean array und update array do not match."
        msg += f"\nShape mean array:\t {old_mean_vector.shape}"
        msg += f"\nShape update array:\t {update_vector.shape}"
        raise ValueError(msg)
    return samples/(samples + 1)*old_mean_vector + 1/(samples+1)*update_vector

def _update_square(samples: int, old_square_vector: np.ndarray, update_vector: np.ndarray) -> np.ndarray:
    """Update the normalised second moment."""
    if not old_square_vector.shape == update_vector.shape:
        msg = "Shape of mean array und update array do not match."
        msg += f"\nShape mean array:\t {old_square_vector.shape}"
        msg += f"\nShape update array:\t {update_vector.shape}"
        raise ValueError(msg)
    return samples/(samples + 1)*old_square_vector + 1/(samples+1)*np.power(update_vector,2)



class PointStatisticsObject:
    """Class that contains utilities for the computation of mean, second moment, and deviation of 'ref -> time -> functionAtPoint' dictionary."""
    def __init__(self, name: str, ref_to_time_grid: dict[int,list[float]], point: list[float]):
        self.name = name
        self.point = point
        self.point_dim = len(point)
        self.ref_to_time_grid = ref_to_time_grid
        self.ref_to_time_to_funcAtpoint_mean = {level: {time: np.ndarray((self.point_dim,)) for time in ref_to_time_grid[level]} for level in ref_to_time_grid.keys()}
        self.ref_to_time_to_funcAtpoint_square = {level: {time: np.ndarray((self.point_dim,)) for time in ref_to_time_grid[level]} for level in ref_to_time_grid.keys()}
        self.samples = 0
        

    def update(self,ref_to_time_to_function) -> None:
        """Add a sample to mean and second moment."""
        for level in self.ref_to_time_to_funcAtpoint_mean.keys():
            for time in self.ref_to_time_to_funcAtpoint_mean[level].keys():
                #update mean
                self.ref_to_time_to_funcAtpoint_mean[level][time][:] = _update_mean(self.samples,
                                                                                       self.ref_to_time_to_funcAtpoint_mean[level][time],
                                                                                       ref_to_time_to_function[level][time].at(self.point))
                #update square
                self.ref_to_time_to_funcAtpoint_square[level][time][:] = _update_square(self.samples,
                                                                                       self.ref_to_time_to_funcAtpoint_square[level][time],
                                                                                       ref_to_time_to_function[level][time].at(self.point))
        self.samples += 1
    
    @property
    def ref_to_time_to_funcAtpoint_deviation(self) -> dict[int,dict[float,np.ndarray]]:
        """Compute the standard deviation based on mean and second moment."""
        ref_to_time_to_funcAtpoint_dev = {level: {time: np.ndarray((self.point_dim,)) for time in self.ref_to_time_grid[level]} for level in self.ref_to_time_grid.keys()}
        for level in self.ref_to_time_to_funcAtpoint_mean.keys():
            for time in self.ref_to_time_to_funcAtpoint_mean[level].keys():
                #Analytically the deviation is always non-negative and there is no issue to take the square root. 
                #But numerically it might happen that the deviation is negative, leading to difficulties. We avoid this by first truncating negative values to 0.
                dev = self.ref_to_time_to_funcAtpoint_square[level][time] - np.power(self.ref_to_time_to_funcAtpoint_mean[level][time],2)
                zeros = np.zeros(dev.shape)
                dev_zeros = np.maximum(dev,zeros)
                ref_to_time_to_funcAtpoint_dev[level][time][:] = np.sqrt(dev_zeros)
        return ref_to_time_to_funcAtpoint_dev

    def _save(self, name_directory: str) -> None:
        """Save mean and standard deviation in csv format."""
        header = ["time"] + [f"Mean_{comp}" for comp in range(self.point_dim)] + [f"SD_{comp}" for comp in range(self.point_dim)]
        ref_to_data = {}
        for level in self.ref_to_time_to_funcAtpoint_mean.keys():
            ref_to_data[level] = [[time] + [self.ref_to_time_to_funcAtpoint_mean[level][time][comp] for comp in range(self.point_dim)] + [self.ref_to_time_to_funcAtpoint_deviation[level][time][comp] for comp in range(self.point_dim)]
                                   for time in self.ref_to_time_to_funcAtpoint_mean[level].keys()]
            
        if not os.path.isdir(name_directory):
            os.makedirs(name_directory)

        for level in self.ref_to_time_to_funcAtpoint_mean.keys():
            outfile = name_directory + "/refinement_" + str(level) + ".csv"
            with open(outfile,"w",newline="") as file:
                writer = csv.writer(file)
                writer.writerow(header)
                writer.writerows(ref_to_data[level])


    def save(self,name_directory: str) -> None:
        """Save mean and standard deviation in vtk format."""
        save_directory = name_directory + "/" + self.name
        #remove old data
        if os.path.isdir(save_directory):
                shutil.rmtree("./" + save_directory)
        self._save(save_directory)