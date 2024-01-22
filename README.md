# Atrial Fibrillation Detector
## Paper
This project implements an atrial fibrillation (AF) detector based on the paper titled "Low-complexity detection of atrial fibrillation in continuous long-term monitoring" [^1].

## Dataset
The algorithm was developed and tested on the Long Term Atrial Fibrillation (LTAF) database [^2][^3], consisting of 84 ECG recordings from patients with paroxysmal or persistent AF, most with a 24-hour duration. The entire database comprises nearly 9 million beats, with 59% occurring during AF.

> [!IMPORTANT]
> In order to be able to run the code, make sure to:
> 1) Download the WFDB package from [here](https://physionet.org/physiotools/matlab/wfdb-app-matlab/wfdb-app-toolbox-0-10-0.zip), unzip the folder, and put the contents of `mcode` folder into the project directory as `<LTAF-detection>/packages/wfdb`.
> 2) Download the LTAF dataset from [here](https://physionet.org/static/published-projects/ltafdb/long-term-af-database-1.0.0.zip) and put all the contents in `<LTAF-detection>/data/LTAF`.

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
