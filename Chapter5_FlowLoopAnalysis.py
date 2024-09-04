# -*- coding: utf-8 -*-
"""
Script to monitor data during T2 reaction monitoring using Magritek Spinsolve spectrometer

If wanting to monitor data while expt is underway, do not use while T2Bulk expt is running
- will cause error

Script will not work if only first expt has completed, with no delays completed

User needs to input:
    exptfolderpath (Line 27) - directory with experiment data (can also be done with user input, see Lines 20 & 28)
    ExptOffset (Line 31) - time between start of "reaction" and start of NMR experiments (in minutes)
    saveFolder (Line 34) - folder where generated plots are saved

Created by Mark I. Grimes, 01 MAY 2024

"""
#%% Find folders for experimental data and initialise variables

#from tkinter import filedialog    # Uncomment if want user to select folder containing expt data
import os
import numpy as np
import matplotlib.pyplot as plt
import datetime

## Choose folder containing expt data - must be clean, only expt folders of interest
exptfolderpath = r"C:\Users\MRRC User\Desktop\MarkGrimes\Expts\MG-E030\Bioreactor run\Expt"
#exptfolderpath = filedialog.askdirectory()   # Opens window for user to select folder

## Input offset between start of "reaction" and start of experiments in minutes
ExptOffset = 60

## Select save folder for plot images when generated
saveFolder = r"C:\Users\MRRC User\Desktop\MarkGrimes\Expts\MG-E030\Bioreactor run\Image folder"

## Pulls out all folder names within expt data folder
AllFolders = [folder for folder in os.listdir(exptfolderpath) if os.path.isdir(os.path.join(exptfolderpath, folder))]

## Keep only folders of interest
X = [folder for folder in AllFolders if "Batch" not in folder]
subFolders = X

## Find position of last T2Bulk expt
LastInst = max([i for i, s in enumerate(subFolders) if 'T2Bulk' in s])
n = int(LastInst/3) # Number of expts is (n + 1)

## Remove any folders after last T2Bulk expt
subFolders[LastInst+1:] = []

## Create array of data folders
X2 = [folder for folder in subFolders if "T2Bulk" in folder]
DataFolders = X2

## Create array of wait time folders
X3 = [folder for folder in subFolders if "WaitTime" in folder]
WaitFolders = X3
n2 = len(WaitFolders)

del X, X2, X3   # Remove unneeded variables

## Create general file path for functions
genfilepath = exptfolderpath + "\\"

## Initialise variables
datapathname = []
T2Vals = np.zeros(len(DataFolders))
CIVals = np.zeros([len(DataFolders),2])
ExptDurVals = np.zeros(len(DataFolders))
WaitDurVals = np.zeros(len(WaitFolders))

print(' ')
print('Data folder names loaded into script.')
print(' ')

#%% Functions to pull out data

## Function to pull out and fit T2 data and acquisition parameters
def T2BulkAnalysis_funcData(pathname):
    
    from scipy.optimize import curve_fit
    
    ## Pull out acquisition parameters
    filename = os.path.join(pathname, 'acqu.par')
    with open(filename, 'r') as file:
        lines = file.readlines()
        
    ## Put acquisition parameters into dictionary
    acqu = {}
    for line in lines:
            key, value = line.strip().split('=')
            acqu[key.strip()] = value.strip()
                     
    ## Pull out data from pt1 file
    filepath = os.path.join(pathname, 'T2CPMG.pt1')
    with open(filepath, 'rb') as file:
        FirstData = np.fromfile(file, dtype=np.float32)
    
    ## Trim excess points
    nEchoes = int(acqu['nrEchoes'])
    data = FirstData[(nEchoes + 419)::2] 
    
    ## Create x-axis values
    x = np.arange(1, len(data) + 1)
    tEcho = float(acqu['echoTime'])
    t = x * tEcho * 1e-6
    
    ## Define the fitting function
    def F(t, A, T2):
        return A * 1e2 * np.exp(-t / T2)    # 1e2 added as scaling factor, may need to be changed as req'd

    ## Initial guesses for fitting parameters - second variable is guess for T2 in s
    x0 = [1, 2.5]

    ## Perform the curve fitting
    popt, pcov = curve_fit(F, t, data, p0=x0)
    
    ## Generate 95% confidence intervals of T2 value
    perr = np.sqrt(np.diag(pcov))
    CI = popt[1] - 1.96 * perr[1], popt[1] + 1.96 * perr[1]

    ## Extract T2 value
    T2 = popt[1]

    return T2, CI


#%% Collate experimental data

## Collect T2 data, CIs and expt length in array
for k in range(n+1):
    datapathname_k = os.path.join(genfilepath, DataFolders[k])
    T2, CI = T2BulkAnalysis_funcData(datapathname_k)
    T2Vals[k] = T2
    CIVals[k,:] = CI
    print('Data from expt ' + str(k+1) + ' loaded.')
    
del T2, CI    # Remove unneeded variables

## Generate R2 values, and associated CIs
R2 = 1./T2Vals
CIVals_R2 = 1./CIVals
CI_R2 = (CIVals_R2[:,0] - CIVals_R2[:,1])/2

print(' ')
print('Data analysis completed.')

#%% Determine experimental time points

FolderTS = []
TSFolders = AllFolders
TSFolders[LastInst+2:] = [] # Remove any folders after last T2Bulk expt

## Generate Python-compatible timestamps for each folder using timestamp in folder name
for k in range(len(AllFolders)):
    tempTS = TSFolders[k].split(" ")
    tempTS = tempTS[0]
    tempTS = '20' + tempTS[0:2] + '-' + tempTS[2:4] + '-' + tempTS[4:6] + ' ' + tempTS[7:9] + ':' + tempTS[9:11] + ':' + tempTS[11:13]
    tempTS = datetime.datetime.strptime(tempTS,"%Y-%m-%d %H:%M:%S")
    FolderTS.append(tempTS)
    
TimeVals = np.zeros(len(AllFolders))
TimeVals[0] = ExptOffset*60 # Set first time value as expt offset

## Determine other time values from timestamps
for k in range(1, len(FolderTS)):
    tempTimeVal = FolderTS[k] - FolderTS[k-1]
    tempTimeVal = tempTimeVal.total_seconds()
    TimeVals[k] = tempTimeVal
    
## Set first time point value as expt offset + time between batch start and first expt start
timepoints = np.zeros(len(DataFolders))
timepoints[0] = TimeVals[0] + TimeVals[1]
DetTimeVals = TimeVals[2:]

## Combine remaining time values to create time points up to last T2Bulk expt
for k in range(1,len(DataFolders)):
    timepoints[k] = DetTimeVals[(k*3)-1] + DetTimeVals[(k*3)-2] + DetTimeVals[(k*3)-3] + timepoints[k-1]
    
del tempTS, tempTimeVal, DetTimeVals    # Remove unneeded variables

print(' ')
print('Expt time points calculated.')

    
#%% Display data and save image of generated plot

## Generate time for plot title and image filename
ct = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
plotname = (ct + '.png')
plotname = plotname.replace(":", "-")

## Plot data, with x-axis units changed based on maximum time point
## Saves image of plot with time stamp as filename in chosen folder
if max(timepoints) < 10800:  # t < 3 h, plot in minutes
    t_plot = timepoints/60
    LimBuffer = (max(t_plot)*1.02) - max(t_plot)
    plt.figure(dpi = 800)
    plt.plot(t_plot, R2, 'o:k', mec = '#00B1C1', mfc = '#00B1C1')   # Dotted line added to guide the eye
    plt.errorbar(t_plot, R2, yerr = CI_R2, ecolor = 'k', linestyle = 'none')
    plt.xlabel('Time [min]')
    plt.ylabel('$R_2$($^1$H$_2$O) [s$^{-1}$]')
    plt.xlim((0 - LimBuffer), (max(t_plot) + LimBuffer))
    plt.title(ct)
    plt.savefig(os.path.join(saveFolder, plotname))
        
elif max(timepoints) >= 10800 and max(timepoints) < 86400:  # 3 h <= t < 24 h, plot in hours
    t_plot = timepoints/3600
    LimBuffer = (max(t_plot)*1.02) - max(t_plot)
    plt.figure(dpi = 800)
    plt.plot(t_plot, R2, 'o:k', mec = '#00B1C1', mfc = '#00B1C1')   # Dotted line added to guide the eye
    plt.errorbar(t_plot, R2, yerr = CI_R2, ecolor = 'k', linestyle = 'none')
    plt.xlabel('Time [h]')
    plt.ylabel('$R_2$($^1$H$_2$O) [s$^{-1}$]')
    plt.xlim((0 - LimBuffer), (max(t_plot) + LimBuffer))
    plt.title(ct)
    plt.savefig(os.path.join(saveFolder, plotname))
            
elif max(timepoints) >= 86400:  # t >= 24 h, plot in days
    t_plot = timepoints/86400
    LimBuffer = (max(t_plot)*1.02) - max(t_plot)
    plt.figure(dpi = 800)
    plt.plot(t_plot, R2, '.:k', mec = '#00B1C1', mfc = '#00B1C1')   # Dotted line added to guide the eye
    plt.errorbar(t_plot, R2, yerr = CI_R2, ecolor = 'k', linestyle = 'none')
    plt.xlabel('Time [days]')
    plt.ylabel('$R_2$($^1$H$_2$O) [s$^{-1}$]')
    plt.xlim((0 - LimBuffer), (max(t_plot) + LimBuffer))
    plt.title(ct)
    plt.savefig(os.path.join(saveFolder, plotname))

print(' ')
print('Plot generated at ' + ct +'.')
print('Save location: ' + saveFolder)    