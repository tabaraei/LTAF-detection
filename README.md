# Long-term Atrial Fibrillation Detection

## Installation
- This project was developed and tested on a __Windows__ operating system. When trying to install the WFDB package on __MacOS__, some unexpected errors encountered and could not replicate the same results on MacOS. So, <ins>it is recommended to run the project on Windows</ins>.
- Clone the project using `git clone https://github.com/tabaraei/LTAF-detection.git`, or directly download the zip file, and store the `LTAF-detection` codebase on desired directory in your local machine.
> [!IMPORTANT]
> In order to be able to run the code, make sure to:
> 1. In MATLAB console, go to the current directory using `cd <PATH-TO-LTAF-detection>` where the `main.m` file exists.
> 2. Download the WFDB package from [here](https://physionet.org/physiotools/matlab/wfdb-app-matlab/wfdb-app-toolbox-0-10-0.zip), unzip the folder, and put the contents of `mcode` folder into the project directory as `<PATH-TO-LTAF-detection>/packages/WFDB`.
> 3. Download the LTAF dataset from [here](https://physionet.org/static/published-projects/ltafdb/long-term-af-database-1.0.0.zip) and put all the contents in `<PATH-TO-LTAF-detection>/data/LTAF`.

## Paper
This project implements an atrial fibrillation (AF) detector based on the paper titled "Low-complexity detection of atrial fibrillation in continuous long-term monitoring" [^1]. There are 5 main stages to detect the AF episodes as propsed below:
__1. Preprocessing (forward_backward_averager):__ This function computes the forward-backward filtering to achieve a linear phase on a given series X, and it aims to compute the exponential average to better track the trend using 0<alpha<1 as the degree of smoothing. This procedure results in the estimation of the mean RR interval, which can be used as a feature in the AF detector.
__2. Preprocessing (median_filter):__ In order to reduce the effect of ectopic beats in the RR series, we may use a simple 3-point median filter as `rm(n)=median{r(n-1),r(n),r(n+1)}`. This filter is also useful to reject the outlier RR intervals.
__3. RR Irregularity Detection (irregularity_detector):__ In order to distinguish RR irregularities, we should use a sliding detection window of length N, located at time n, and compute the number of all pairwise RR interval combinations differing more than gamma seconds, and normalize them with their maximum value. To achieve this, first, we will introduce the heaviside function "H", where it simply denotes whether the pairwise difference of given RR intervals are below a certain threshold "gamma" or not. Then, for each RR interval "n", we will compute the count of pairwise RR interval combinations differing more than "gamma", and we will normalize this value w.r.t its maximum value. Lastly, We compute the ratio "I" between the smoothed version of M utilizing exponential averaging, and the RR interval trend stored in "rt".
__4. Bigeminy Supression (bigeminy_supressor):__ Bigeminy is a cardiac arrhythmia characterized by a pattern in which  every normal heartbeat is followed by a premature beat, creating a regular pattern of two beats, which can be incorrectly interpreted as AF when the detection is RR-based.
__5. Signal Fusion and Detection (signal_fusion):__ The decision function O(n) is produced through a simple signal fusion, which is identical to "Bt" unless it exceeds a fixed threshold delta when, instead it becomes identical to "It".

## Dataset
The algorithm was developed and tested on the Long Term Atrial Fibrillation (LTAF) database [^2][^3], consisting of 84 ECG recordings from patients with paroxysmal or persistent AF, mostly with a 24-hour duration. The entire database comprises nearly 9 million beats, with 59% occurring during AF.



## Algorithm
![ECG 1](https://github.com/tabaraei/LTAF-detection/blob/master/plots/ecg1.png)
![ECG 2](https://github.com/tabaraei/LTAF-detection/blob/master/plots/ecg2.png)
![RR](https://github.com/tabaraei/LTAF-detection/blob/master/plots/rr.png)
![Median Filter](https://github.com/tabaraei/LTAF-detection/blob/master/plots/medianfilter.png)
![Averager](https://github.com/tabaraei/LTAF-detection/blob/master/plots/averager.png)
![M and I](https://github.com/tabaraei/LTAF-detection/blob/master/plots/M(n)%20%26%20I(n).png)
![B](https://github.com/tabaraei/LTAF-detection/blob/master/plots/B(n).png)
![O](https://github.com/tabaraei/LTAF-detection/blob/master/plots/O(n).png)
![O zoom](https://github.com/tabaraei/LTAF-detection/blob/master/plots/O(n)%20zoom.png)
![Results0](https://github.com/tabaraei/LTAF-detection/blob/master/plots/results0.png)
![Results1](https://github.com/tabaraei/LTAF-detection/blob/master/plots/results1.png)
![Results2](https://github.com/tabaraei/LTAF-detection/blob/master/plots/results2.png)
![Results3](https://github.com/tabaraei/LTAF-detection/blob/master/plots/results3.png)
![Results4](https://github.com/tabaraei/LTAF-detection/blob/master/plots/results4.png)

> [!NOTE]
> Useful information that users should know, even when skimming content.

> [!TIP]
> Helpful advice for doing things better or more easily.

> [!WARNING]
> Urgent info that needs immediate user attention to avoid problems.

> [!CAUTION]
> Advises about risks or negative outcomes of certain actions.

Here is a simple footnote[^1].

A footnote can also have multiple lines[^2].

[^1]: A. Petrėnas, V. Marozas, L. Sörnmo, Low-complexity detection of atrial fibrillation in continuous long-term monitoring, Computers in Biology and Medicine, 65, 184-191, 2015
[^2]: S. Petrutiu, A.V. Sahakian, S. Swiryn, Abrupt changes in fibrillatory wave characteristics at the termination of paroxysmal atrial fibrillation in humans, Europace 9 (2007) 466–470
[^3]: A.L. Goldberger, L.A. Amaral, L. Glass, J.M. Hausdorff, P.C. Ivanov, R.G. Mark, J.E. Mietus, G.B. Moody, C.K. Peng, H.E. Stanley, PhysioBank, PhysioToolkit, and PhysioNet: components of a new research resource for complex physiologic signals, Circulation 101 (2000) E215–E220
