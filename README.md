# H2_ALSH: Homocentric Hypersphere partition Asymmetric LSH

Version: 1.0.0

Release date: 27-05-2018


Introduction
--------

This package is written in the C++ programming language. It provides an 
internal Asymmetric LSH scheme based on Homocentric Hypersphere partition
(H2-ALSH) for c-Approximate Maximum Inner Product (or simply c-AMIP) 
search.


How to use this Package?
--------

We provide a Makefile and a script (i.e., run_mnist.sh) as a running example 
for comipling and running this package. Before start running this package, 
please ensure the input format of the dataset and query set is correct. We 
provide a sample dataset and query set (i.e., Mnist) for your reference.

We also provide the scripts (i.e., run_sift.sh, run_gist.sh, run_netflix.sh, 
run_yahoo.sh, and run_para.sh) for the users who would like to reproduce our 
results presented in SIGKDD 2018. The datasets Sift, Gist, Netflix, and Yahoo 
we used can be downloaded from the following links:

* Sift: https://drive.google.com/open?id=1dAFbjQWoBIAW30lGzTXaf_w7DxBSUzoN

* Gist: https://drive.google.com/open?id=1r1rsSm6-IdWX2-8eFJkChFfP0Ej7TYM4

* Netflix: https://drive.google.com/open?id=1bJQftqxlC8u4ijDf5gpEnw1tJE2nLfoG

* Yahoo: https://drive.google.com/open?id=18k0ISgjtQhHHqoGi8A96Fm-Q2gX_jxgw


Authors
--------

* **Qiang Huang**

  Smart Systems Institute, National University of Singapore (NUS),
  
  Singapore, 119613 
  
  huangq2011@gmail.com, huangq25@mail2.sysu.edu.cn
  
  https://sites.google.com/site/qianghuang2017/
  

* **Guihong Ma**

  School of Data and Computer Science, Sun Yat-Sen University (SYSU),
  
  Guangzhou, China, 510006
  
  maguihong@vip.qq.com


Relevant Papers
--------

The paper for the package of H2-ALSH has been published in SIGKDD 2018 which is 
displayed as follows:

* **Qiang Huang, Guihong Ma, Jianlin Feng, Qiong Fang, and Anthony K. H. Tung. 
Accurate and Fast Asymmetric Locality-Sensitive Hashing Scheme for Maximum 
Inner Product Search. SIGKDD, pages 1561-1570, 2018.**

If you use the package for publications, please cite the paper above.
