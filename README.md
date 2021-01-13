## Hurricane Maria Impacts on Forests
This script is used for producing the forest disturbance intensity map after Hurricane Maria made landfall on Puerto Rico island on Sep. 20th, 2017.
The script produced first part of the results in this peer-reviewed [paper](https://www.sciencedirect.com/science/article/abs/pii/S0034425720303102) and this [preprint](https://peerj.com/preprints/26597/).

Preprocessing steps of Landsat images include:
1. cloud masking
2. topographic illumination correction
3. radiometric calibration using invariant targets

#### Landsat 8 images before and after Hurricane Maria
Pre-hurricane image
![pre-hurricane](https://github.com/ylfeng93/Hurricane-Maria-Impacts-on-Forests/blob/main/docs:imgs/Screen%20Shot%202021-01-13%20at%2010.19.18%20AM.png)
Post-hurricane image
![post-hurricane](https://github.com/ylfeng93/Hurricane-Maria-Impacts-on-Forests/blob/main/docs:imgs/Screen%20Shot%202021-01-13%20at%2010.19.31%20AM.png)

Specture mixture analysis is a great measurement that split using remote sensing images to the fractions of endmembers. Landsat-derived DNPV is a good representation of forest disturbance intensity.

Hurricane Disturbance Intensity Metric - DNPV
![metrics](https://github.com/ylfeng93/Hurricane-Maria-Impacts-on-Forests/blob/main/docs:imgs/Screen%20Shot%202021-01-13%20at%2010.20.26%20AM.png)
