# Lampropholis temperature and yolk experiment

Understanding the impact of maternal resource investment and temperature on lizard behaviour and performance 

# Navigation
All analyses can be reproduced using the `lizard_behav_analysis.R` file in the `R/` folder. Functions that are used for analyses are sourced within this file but can be found in the `func.R` script in the `R` folder.

# Data files
Raw data used for analyses is in the `data` folder named `dataset_reclean_git.csv`. These data contain all the behavioural and performance measurements for lizards. The `data` folder also contains the `morphol.csv` file which is used for alll morphology analyses. 

The `dataset_reclean_git.csv` data file is processed using code in the `lizard_behav_analysis.R` script. Two processed data files, cleaned and ready for analysis, can us used to reproduce analyses. These are found in the `output/data/` folder and are named as follows:

- `dat2.csv`: This file contains the performance, behaviour, treatment information, lizard id and body size (SVL) measurements for *Lampropholis delicata*.
- `dat3.csv`: This file contains the performance, behaviour, treatment information, lizard id and body size (SVL) measurements for *Lampropholis guichenoti*.

# Column names
There are a series of column names within the data files that are relevant. These are detailed below:
filename: `morphol.csv`
- "id" : Unique lizard identity       
- "clutch" : Unique clutch identity
- "temp" : Temperature treatment of the lizard ("cold" or "hot")     
- "egg_treat" : Egg treatment of the lizard ("control" or "ablated")
- "sp" : Species, Deli = *Lampropholis delicata* and Guich = *Lampropholis guichenoti*
- "Weigth"  : Weight in grams of each lizard
- "SVL"  : Snout-vent length (mm) of each lizard     
- "Tail" : Tail length (mm) of each lizard       
- "Total" : Total length (mm) of each lizard. SVL + Tail     
- "age" : Age, in days, of each lizard

The behaviour and performance data come from merged data files and have a bit more information:

filename: `dat2.csv` and `dat3.csv`
- "id"	: Unique lizard identity
- "clutch": Unique clutch identity
- "day_hatch"
- "day_photo"
- "temp" : Temperature treatment of the lizard ("cold" or "hot")  
- "egg_treat" : Egg treatment of the lizard ("control" or "ablated")
- "sp" : Species, Deli = *Lampropholis delicata* and Guich = *Lampropholis guichenoti*
- "mass.hatch" : Weight in grams of each lizard at hatching
- "SVL_hatch" : Snout-vent length (mm) of each lizard at hatching
- "Tail_length_hatch" Tail length (mm) of each lizard at hatching
- "Weigth": Weight in grams of each lizard
- "Distance.moved": Total distance moved during assay in cm
- "Time_hide_sec": Time taken for lizard to enter the hide (seconds)
- "Time_snout_sec" : Time taken for lizard to emerge snout from the hide (seconds)
- "Time_emerge_sec" : Time taken for lizard to completely emerge from the hide (seconds)
- "SVL": SVL of lizard at assay time (mm)
- "Tail": Tail length of lizard at assay time (mm)
- "Total": Total length of lizard at assay time (mm)
- "Tail_intact": Binary, did the lizard have an intact tail or did it lose it
- "day": day of assay
- "speed_1m_s": 1m sprint speed (m/s)
- "burst_25cm": fastest 25cm speed (cm/s)
- "X25fast"
- "age"	: Age of lizard in days
- "Notes": any relevant notes
- "time_tohide":  Time taken for lizard to enter the hide (seconds)
- "time_hiding":Time taken for lizard to emerge snout from the hide (seconds)
- "time_active": Time taken for lizard to completely emerge from the hide (seconds)
- "maternal": Egg treatment of the lizard ("control" or "ablated")
- "scaleage": Z-transformed age
- "svldeli": SVL of L delicata
- "logTimeSnout": Time taken for lizard to emerge snout from the hide (seconds) on a log scale
- "logspeed_1m":1m sprint speed (m/s) on log scale
- "logspeed_burst": fastest 25cm speed (cm/s) on log scale
- "logTime_emerge_sec":  Time taken for lizard to completely emerge from the hide (seconds) on a log scale
- "z_svl": Z-transformed SVL