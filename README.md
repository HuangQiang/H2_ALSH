# H2_ALSH: Homocentric Hypersphere partition Asymmetric LSH

## Introduction

This package provides several internal asymmetric LSH schemes for ```k-Maximum Inner Product Search (k-MIPS)```. It is also an implementation of ```H2_ALSH``` from the paper as follows:
```
Qiang Huang, Guihong Ma, Jianlin Feng, Qiong Fang, and Anthony K. H. Tung. 
Accurate and Fast Asymmetric Locality-Sensitive Hashing Scheme for Maximum 
Inner Product Search. SIGKDD, pages 1561-1570, 2018.
```

## Compilation

The package requires ```g++``` with ```c++11``` support. To download and compile the code, type:
```
$ git clone https://github.com/HuangQiang/H2_ALSH.git
$ cd H2_ALSH
$ make
```

## Datasets

We four real-life datasets [```Sift```](https://drive.google.com/open?id=1dAFbjQWoBIAW30lGzTXaf_w7DxBSUzoN), [```Gist```](https://drive.google.com/open?id=1r1rsSm6-IdWX2-8eFJkChFfP0Ej7TYM4), [```Netflix```](https://drive.google.com/open?id=1bJQftqxlC8u4ijDf5gpEnw1tJE2nLfoG), and [```Yahoo```](https://drive.google.com/open?id=18k0ISgjtQhHHqoGi8A96Fm-Q2gX_jxgw) for comparison. The statistics of the datasets are summarized in the following table:

| Datasets | #Objects | #Queries | Dimensionality | Data Size |
| -------- | -------- | -------- | --------       | --------  |
| Sift     | ```1,000,000``` | ```1,000```  | 128 | 337.8 MB  |
| Netflix  | ```17,770```    | ```1,000```  | 300 | 50.8 MB   |
| Yahoo    | ```624,961```   | ```1,000```  | 300 | 2.3 GB    |
| Gist     | ```1,000,000``` | ```1,000```  | 960 | 4.0 GB    |

## Run Experiments

```
Usage: alsh [OPTIONS]

This package supports 12 options to evaluate the performance of H2_ALSH, L2_ALSH, L2_ALSH2, XBOX, Sign_ALSH, Simple_LSH and Linear_Scan for k-MIPS. The parameters are introduced as follows.

  -alg    integer    options of algorithms (0 - 11)
  -n      integer    cardinality of dataset
  -d      integer    dimensionality of dataset and query set
  -qn     integer    number of queries
  -K      integer    number of hash tables for Sign_ALSH and Simple_LSH
  -m      integer    extra dimension for L2_ALSH, L2_ALSH2, and Sign_ALSH
  -U      float      a value in (0,1] for L2_ALSH, L2_ALSH2, and Sign_ALSH
  -c0     float      approximation ratio for NN Search (c0 > 1)
  -c      float      approximation ratio for MIP Search (0 < c < 1)
  -ds     string     address of data  set
  -qs     string     address of query set
  -ts     string     address of truth set
  -op     string     output path
```

We provide all scripts to repeat all experiments. A quick example is shown as follows (run ```H2_ALSH``` on ```Mnist```):
```
./alsh -alg 1 -n 60000 -qn 1000 -d 50 -c0 2.0 -c 0.5 -ds data/Mnist/Mnist.ds -qs data/Mnist/Mnist.q -ts data/Mnist/Mnist.mip -op results/Mnist/
```

If you would like to get more information to run other algorithms, please check the scripts in the package. When you run the package, please ensure that the path for the dataset, query set, and truth set is correct. Since the package will automatically create folder for the output path, please keep the path as short as possible. 


## Related publication

If you use this package for publications, please cite the paper as follows.
```
@inproceedings{huang2018accurate,
    title={Accurate and Fast Asymmetric Locality-Sensitive Hashing Scheme for Maximum Inner Product Search}
    author={Huang, Qiang and Ma, Guihong and Feng, Jianlin and Fang, Qiong and Tung, Anthony KH},
    booktitle={Proceedings of the 24th ACM SIGKDD International Conference on Knowledge Discovery \& Data Mining},
    pages={1561--1570},
    year={2018},
    organization={ACM}
}
```
