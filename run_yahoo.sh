#!/bin/bash
make
make clean

# ------------------------------------------------------------------------------
#  Parameters
# ------------------------------------------------------------------------------
dname=Yahoo
n=624961
d=300

qn=1000
K=512
m=3
U1=0.83
U2=0.85
c0=2.0
c=0.5
dPath=./data/${dname}/${dname}
oFolder=./results/${dname}/

# ------------------------------------------------------------------------------
#  Ground-Truth
# ------------------------------------------------------------------------------
./alsh -alg 0 -n ${n} -qn ${qn} -d ${d} -ds ${dPath}.ds -qs ${dPath}.q \
    -ts ${dPath}.mip

# ------------------------------------------------------------------------------
#  Algorithms for c-k-AMIP search
# ------------------------------------------------------------------------------
./alsh -alg 1 -n ${n} -qn ${qn} -d ${d} -c0 ${c0} -c ${c} -ds ${dPath}.ds \
    -qs ${dPath}.q -ts ${dPath}.mip -of ${oFolder}

./alsh -alg 2 -n ${n} -qn ${qn} -d ${d} -m ${m} -U ${U1} -c0 ${c0} \
    -ds ${dPath}.ds -qs ${dPath}.q -ts ${dPath}.mip -of ${oFolder}

./alsh -alg 3 -n ${n} -qn ${qn} -d ${d} -m ${m} -U ${U1} -c0 ${c0} \
    -ds ${dPath}.ds -qs ${dPath}.q -ts ${dPath}.mip -of ${oFolder}

./alsh -alg 4 -n ${n} -qn ${qn} -d ${d} -c0 ${c0} -ds ${dPath}.ds \
    -qs ${dPath}.q -ts ${dPath}.mip -of ${oFolder}

./alsh -alg 5 -n ${n} -qn ${qn} -d ${d} -K ${K} -m ${m} -U ${U2} -c0 ${c0} \
    -ds ${dPath}.ds -qs ${dPath}.q -ts ${dPath}.mip -of ${oFolder}

./alsh -alg 6 -n ${n} -qn ${qn} -d ${d} -K ${K} -c0 ${c0} -ds ${dPath}.ds \
    -qs ${dPath}.q -ts ${dPath}.mip -of ${oFolder}

./alsh -alg 7 -n ${n} -qn ${qn} -d ${d} -ds ${dPath}.ds -qs ${dPath}.q \
    -ts ${dPath}.mip -of ${oFolder}

# ------------------------------------------------------------------------------
#  Precision-Recall Curves of Algorithms for c-k-AMIP search
# ------------------------------------------------------------------------------
./alsh -alg 8 -n ${n} -qn ${qn} -d ${d} -c0 ${c0} -c ${c} -ds ${dPath}.ds \
    -qs ${dPath}.q -ts ${dPath}.mip -of ${oFolder}

./alsh -alg 9 -n ${n} -qn ${qn} -d ${d} -K ${K} -m ${m} -U ${U2} -c0 ${c0} \
    -ds ${dPath}.ds -qs ${dPath}.q -ts ${dPath}.mip -of ${oFolder}

./alsh -alg 10 -n ${n} -qn ${qn} -d ${d} -K ${K} -c0 ${c0} -ds ${dPath}.ds \
    -qs ${dPath}.q -ts ${dPath}.mip -of ${oFolder}

# ------------------------------------------------------------------------------
#  Norm Distribution
# ------------------------------------------------------------------------------
./alsh -alg 11 -n ${n} -d ${d} -ds ${dPath}.ds -of ${oFolder}
