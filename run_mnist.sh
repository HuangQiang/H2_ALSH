#!/bin/bash
make clean
make -j

# ------------------------------------------------------------------------------
#  Parameters
# ------------------------------------------------------------------------------
dname=Mnist
n=60000
d=50

qn=1000
K=512
m=3
U1=0.83
U2=0.85
c0=2.0
c=0.5
dPath=./data/${dname}/${dname}
oPath=./results/${dname}/

# # ------------------------------------------------------------------------------
# #  Ground-Truth
# # ------------------------------------------------------------------------------
# ./alsh -alg 0 -n ${n} -qn ${qn} -d ${d} -ds ${dPath}.ds -qs ${dPath}.q \
#     -ts ${dPath}.mip

# ------------------------------------------------------------------------------
#  Algorithms for c-k-AMIP search
# ------------------------------------------------------------------------------
./alsh -alg 1 -n ${n} -qn ${qn} -d ${d} -c0 ${c0} -c ${c} -ds ${dPath}.ds \
    -qs ${dPath}.q -ts ${dPath}.mip -op ${oPath}

# ./alsh -alg 2 -n ${n} -qn ${qn} -d ${d} -m ${m} -U ${U1} -c0 ${c0} \
#     -ds ${dPath}.ds -qs ${dPath}.q -ts ${dPath}.mip -op ${oPath}

# ./alsh -alg 3 -n ${n} -qn ${qn} -d ${d} -m ${m} -U ${U1} -c0 ${c0} \
#     -ds ${dPath}.ds -qs ${dPath}.q -ts ${dPath}.mip -op ${oPath}

# ./alsh -alg 4 -n ${n} -qn ${qn} -d ${d} -c0 ${c0} -ds ${dPath}.ds \
#     -qs ${dPath}.q -ts ${dPath}.mip -op ${oPath}

# ./alsh -alg 5 -n ${n} -qn ${qn} -d ${d} -K ${K} -m ${m} -U ${U2} \
#     -ds ${dPath}.ds -qs ${dPath}.q -ts ${dPath}.mip -op ${oPath}

# ./alsh -alg 6 -n ${n} -qn ${qn} -d ${d} -K ${K} -ds ${dPath}.ds \
#     -qs ${dPath}.q -ts ${dPath}.mip -op ${oPath}

# ./alsh -alg 7 -n ${n} -qn ${qn} -d ${d} -ds ${dPath}.ds -qs ${dPath}.q \
#     -ts ${dPath}.mip -op ${oPath}

# # ------------------------------------------------------------------------------
# #  Precision-Recall Curves of Algorithms for c-k-AMIP search
# # ------------------------------------------------------------------------------
# ./alsh -alg 8 -n ${n} -qn ${qn} -d ${d} -c0 ${c0} -c ${c} -ds ${dPath}.ds \
#     -qs ${dPath}.q -ts ${dPath}.mip -op ${oPath}

# ./alsh -alg 9 -n ${n} -qn ${qn} -d ${d} -K ${K} -m ${m} -U ${U2} \
#     -ds ${dPath}.ds -qs ${dPath}.q -ts ${dPath}.mip -op ${oPath}

# ./alsh -alg 10 -n ${n} -qn ${qn} -d ${d} -K ${K} -ds ${dPath}.ds \
#     -qs ${dPath}.q -ts ${dPath}.mip -op ${oPath}

# # ------------------------------------------------------------------------------
# #  Norm Distribution
# # ------------------------------------------------------------------------------
# ./alsh -alg 11 -n ${n} -d ${d} -ds ${dPath}.ds -op ${oPath}
