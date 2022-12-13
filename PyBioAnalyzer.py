"""
Created on Wed Nov 20 20:20:00 2019

Originally from Dr.Georg Urtel's and Shunsuke Sumi's code.
Modified by Evgeniia Edeleva.
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
from glob import glob
from scipy import interpolate
#scipy 1.2.1
#python 3.7


class PyBioAnalyzer:
    def __init__(self, folder_name, assay_type):
        """
        folder_name: folder containing output csv files.
        assay_type : ['HS_DNA', 'pico_RNA', 'small_RNA']
        """
        self.folder_name = folder_name
        self.assay_type = assay_type
        if self.assay_type == 'HS_DNA':
            self.a = 2139
        elif self.assay_type == 'pico_RNA':
            self.a = 1438
        elif self.assay_type == 'small_RNA':
            self.a = 1218
        assert assay_type in ['HS_DNA', 'pico_RNA', 'small_RNA'], "Assay type is not supported. Choose 'HS_DNA', 'pico_RNA', or 'small_RNA'."

        BA_results = glob(os.path.join(f"{self.folder_name}", "*.csv"))
        assert BA_results != [], "Could'nt find the data folder! Check the folder location."
        self.BAfiles = BA_results
        
    def _load_bioanalyzer(self, bioanalyzer_file):
        """
        load [time, intensity value] from csv file
        input: path to the BA file 
        output: [[time, intensity],...]
        """
        with open(bioanalyzer_file, "rb") as f:
            l = f.readlines()[18:self.a]
            l = [list(map(float, line.decode("shift-jis").replace("\r\n", "").split(","))) for line in l]
        return pd.DataFrame(l)
    
    def _load_ladder_info(self):
        results = ""
        for file in self.BAfiles:
            if "Results" in file:
                results = file 
        assert results != "", "Couldn't find Results.csv! Check the folder name." 
        with open(results, "rb") as f:
            l = f.readlines()
            index = l.index(b'Sample Name,Ladder\r\n')
            ladder = []
            
            if self.assay_type == 'small_RNA':
                self.unit = "nt"
                size = [4, 20, 40, 60, 80, 100, 150] #from the small RNA Kit Guide
                for i in range(index+4, index+11):
                    line = l[i]
                    line = line.decode().split(",")
                    ladder.append([float(size[i-(index+4)]), float(line[-3])])  
        
            if self.assay_type == 'pico_RNA':
                self.unit = "nt"
                size = [25, 200, 500, 1000, 2000, 4000, 6000] #from the RNA Pico Kit Guide
                for i in range(index+4, index+11):
                    line = l[i]
                    line = line.decode().split(",")
                    ladder.append([float(size[i-(index+4)]), float(line[-3])])
                    
            if self.assay_type == 'HS_DNA':
                self.unit = "bp"
                size = [35, 50, 100, 150, 200, 300, 400, 500, 600, 700, 1000, 2000, 3000, 7000, 10380]
                for i in range(index+4, index+19):
                    line = l[i]
                    line = line.decode().split(",")
                    ladder.append([float(size[i-(index+4)]), float(line[-5])])
            
        ladder = pd.DataFrame(ladder, columns = ["Size", "Time_stamp"])    
        return ladder

    def _time2nucleotide(self):
        """
        linear interpolation.
        output: interpolated function.
        """
        ladder = self._load_ladder_info()
        f = interpolate.interp1d(ladder.Time_stamp, ladder.Size, fill_value = "extrapolate")
        return f

    def _load_all_bioanalyzer(self):
        """
        search bioanalyzer files and ladder file in folder_name and load data from the files.
        input: folder_name(where you saved files)
        output: bioanalyzer (n_sample signal traces), ladder (ladder trace) and time.
        """
        ladder_file = [file for file in self.BAfiles if "Ladder" in file][0]
        bioanalyzer = self._load_bioanalyzer(ladder_file)
        bioanalyzer.columns = ["Time_stamp", "Ladder"]
        for i in range(0, len(self.BAfiles)):
            for file in self.BAfiles:
                if "Sample" + str(i) +".csv" in file:
                    sample_data = self._load_bioanalyzer(file)
                    sample_data.columns = ["Time", file]
                    bioanalyzer = pd.concat([bioanalyzer, sample_data[file]], axis = 1)
        f = self._time2nucleotide()
        Size = [f(time) for time in bioanalyzer.Time_stamp]
        Size = pd.DataFrame(Size, columns = [f"Size[{self.unit}]"])
        bioanalyzer = pd.concat([bioanalyzer, Size], axis = 1)
        return bioanalyzer

    def linearity_check(self):
        ladder = self._load_ladder_info()
        f = self._time2nucleotide()
        plt.plot(ladder.Size, ladder.Time_stamp, marker = "o", color = "steelblue")
        plt.xlabel(f"Size[{self.unit}]", fontsize = 12)
        plt.ylabel("Migration time[s]", fontsize =12)
        plt.grid()
        return 
    
    def plot_samples(self, samples, plotrange, *labels, ladder = True):
        if type(samples) != list:
            samples = [samples]
        if labels == ():
            labels = samples
        elif labels != ():
            labels = labels[0]
        BATable = self._load_all_bioanalyzer()
        # print(BATable.head())
        # print(BATable.columns)
        plt.figure(figsize = (12,6))
        if ladder:
            plt.plot(BATable[f"Size[{self.unit}]"], BATable["Ladder"], color = "crimson", label = "Ladder")
        for f in samples:
            if (not "Ladder" in f) and (not "Results" in f):
                plt.plot(BATable[f"Size[{self.unit}]"], BATable[f], label = os.path.basename(f))
        plt.legend()
        plt.title(f"{self.folder_name}, {self.assay_type} kit", fontsize = 20)
        plt.xlabel(f"Size[{self.unit}]", fontsize = 16)
        plt.ylabel("Intensity", fontsize = 16)
        plt.xlim(plotrange)
        plt.grid(alpha = 0.3)
        return
        
if __name__ == '__main__':
    import argparse
    import matplotlib.pyplot as plt 

    parser = argparse.ArgumentParser()
    parser.add_argument("--folder_name", required = True)
    parser.add_argument("--assay_type", required = True, choices=['HS_DNA', 'pico_RNA', 'small_RNA'])
    parser.add_argument("--min_lim", default = 100, type = int)
    parser.add_argument("--max_lim", default = 500, type = int)
    parser.add_argument("--disable_ladder", action="store_false")
    args = parser.parse_args()

    pba = PyBioAnalyzer(args.folder_name, args.assay_type)
    print("Loaded following sample files:")
    for f in pba.BAfiles:
        print(f"\t{f}")
    pba.plot_samples(pba.BAfiles, [args.min_lim, args.max_lim], ladder = args.disable_ladder)
    plt.show()