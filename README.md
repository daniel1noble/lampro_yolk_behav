# lampro_yolk_behav

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