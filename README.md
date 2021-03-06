# IRIS
#### A method for predicting ensembles of **I**n *vivo*  **R**NA secondary structures using PAR**IS** data ([GitHub Link](https://github.com/JY-Zhou/IRIS))

RNA secondary structures play a pivotal role in posttranscriptional regulation and the functions of non-coding RNAs, yet *in vivo* RNA secondary structures remain enigmatic. PARIS (Psoralen Analysis of RNA Interactions and Structures) is a recently developed high-throughput sequencing-based approach that enables direct capture of RNA duplex structures *in vivo*. However, the existence of incompatible, fuzzy pairing information obstructs the integration of PARIS data with the existing tools in the construction of RNA secondary structure models at the single base resolution. Here, we introduce IRIS, a method for predicting RNA secondary structure ensembles based on PARIS data. IRIS generates a large set of candidate RNA secondary structure models under the guidance of redistributed PARIS reads and then uses a Bayesian model to identify the optimal ensemble, according to both thermodynamic principles and PARIS data. We verified the predicted ensembles based on the evidence from evolutionary conservation and consistency with other experimental RNA structural data. 

## Usage

#### Import the IRIS API into Python scripts

We strongly recommend that users use the IRIS API by writing some simple Python scripts. Here, we prepare three Jupyter notebooks that utilize the IRIS API as examples. If [JupyterLab](https://jupyter.org/) is installed on the system (generally, JupyterLab will be automatically installed when setting up a Python environment using [Anaconda](https://www.anaconda.com/products/individual)), users can make full use of the following notebooks interactively and can freely and flexibly change the codes. 

- [Experiments.ipynb](Experiments.ipynb) 

    This notebook uses the IRIS API to evaluate the performance of IRIS. It includes the main results and figures in the IRIS article, and it is easy to reproduce these results using JupyterLab.

- [Case_study_U2_snRNA.ipynb](Case_study_U2_snRNA.ipynb) 

    This notebook takes the U2 snRNA as an example to perform IRIS, which is the main case discussed in the article for explain how IRIS works. This notebook also helps to reproduce the results and figures in the article concerning the U2 snRNA.

- [Case_study_RMRP.ipynb](Case_study_RMRP.ipynb)
    
    As a supplementary case, the execution of IRIS on the RMRP is included in this notebook. Here, we make the description and code as concise as possible to set up this notebook as an example of how to use the IRIS API in Python scripts.

#### Run IRIS in the command-line interface ([IRIS_CLI.py](IRIS_CLI.py))

We also wrap IRIS as a traditional command-line tool. The usage is shown below.

`python IRIS_CLI.py [-h] -i [FASTA file] -r [BAM file] -o [output file] -K [int] [-C [int]] [-F [float]] [--krange  [...]]`

Optional arguments
 - `-h, --help`       &nbsp;&nbsp; show help messages and exit
 - `-i [FASTA file]`  &nbsp;&nbsp; the input RNA sequence in FASTA format
 - `-r [BAM file]`    &nbsp;&nbsp; the input mapped PARIS reads in BAM format
 - `-o [output file]` &nbsp;&nbsp; the output file of the predicted ensemble
 - `-K [int]`         &nbsp;&nbsp; the expected number of representative structures in the ensemble
 - `-C [int]`         &nbsp;&nbsp; the number of candidate structures generated &nbsp; (*default*: 100)
 - `-F [float]`       &nbsp;&nbsp; the fraction threshold for filtering stems &nbsp; (*default*: 0.75)
 - `--krange  [ ...]` &nbsp;&nbsp; the range of length *k* of stems  &nbsp; (*default*: 4 5 6 7)
 
An example of executing IRIS CLI using default parameters.

`python IRIS_CLI.py -i sample/sample.fa -r sample/sample.bam -o sample/sample.out -K 2`

The parameters can be specified by the following command.

`python IRIS_CLI.py -i sample/sample.fa -r sample/sample.bam -o sample/sample.out -K 2 -C 100 -F 0.75 --krange 4 5 6 7`


## Requirements

IRIS is implemented in Python 3. We suggest users to install the Python 3.x environment using [Anaconda](https://www.anaconda.com/products/individual). Then, the following required Python libraries can be easily installed using `conda`.

| Library | Version | Installation |
|---|---|---|
| viennarna | >= 2.4.14 | `conda install -c bioconda viennarna` |
| pysam | >= 0.15.3 | `conda install -c bioconda pysam` |
| biopython | >= 1.74 | `conda install biopython` |
| numpy | >= 1.17.2 | `conda install numpy` |
| scipy | >= 1.3.1 | `conda install scipy` |
| scikit-learn | >= 0.21.3 | `conda install scikit-learn` |
| matplotlib | >= 3.1.1 | `conda install matplotlib` |
| pillow | >= 6.1.0 | `conda install pillow` |
| tabulate | >= 0.8.5 | `conda install tabulate` |