#!/bin/bash
make
make clean

# ------------------------------------------------------------------------------
#  Tuning approximation ratio c for c-AMIP search
# ------------------------------------------------------------------------------
qn=1000
c0=2.0
for c in 0.1 0.25 0.5 0.8
do
    # --------------------------------------------------------------------------
    #  Sift
    # --------------------------------------------------------------------------
    dname=Sift
    n=1000000
    d=128
    dPath=./data/${dname}/${dname}
    oFolder=./results/${dname}/para/c=${c}/

    ./alsh -alg 1 -n ${n} -qn ${qn} -d ${d} -c0 ${c0} -c ${c} -ds ${dPath}.ds \
        -qs ${dPath}.q -ts ${dPath}.mip -of ${oFolder}

    # --------------------------------------------------------------------------
    #  Gist
    # --------------------------------------------------------------------------
    dname=Gist
    n=1000000
    d=960
    dPath=./data/${dname}/${dname}
    oFolder=./results/${dname}/para/c=${c}/

    ./alsh -alg 1 -n ${n} -qn ${qn} -d ${d} -c0 ${c0} -c ${c} -ds ${dPath}.ds \
        -qs ${dPath}.q -ts ${dPath}.mip -of ${oFolder}

    # --------------------------------------------------------------------------
    #  Netflix
    # --------------------------------------------------------------------------
    dname=Netflix
    n=17770
    d=300
    dPath=./data/${dname}/${dname}
    oFolder=./results/${dname}/para/c=${c}/

    ./alsh -alg 1 -n ${n} -qn ${qn} -d ${d} -c0 ${c0} -c ${c} -ds ${dPath}.ds \
        -qs ${dPath}.q -ts ${dPath}.mip -of ${oFolder}

    # --------------------------------------------------------------------------
    #  Yahoo
    # --------------------------------------------------------------------------
    dname=Yahoo
    n=624961
    d=300
    dPath=./data/${dname}/${dname}
    oFolder=./results/${dname}/para/c=${c}/

    ./alsh -alg 1 -n ${n} -qn ${qn} -d ${d} -c0 ${c0} -c ${c} -ds ${dPath}.ds \
        -qs ${dPath}.q -ts ${dPath}.mip -of ${oFolder}
done

# ------------------------------------------------------------------------------
#  Tuning approximation ratio c0 for c0-ANN search
# ------------------------------------------------------------------------------
qn=1000
c=0.5
for c0 in 1.5 2.0 2.5 3.0
do
    # --------------------------------------------------------------------------
    #  Sift
    # --------------------------------------------------------------------------
    dname=Sift
    n=1000000
    d=128
    dPath=./data/${dname}/${dname}
    oFolder=./results/${dname}/para/c0=${c0}/

    ./alsh -alg 1 -n ${n} -qn ${qn} -d ${d} -c0 ${c0} -c ${c} -ds ${dPath}.ds \
        -qs ${dPath}.q -ts ${dPath}.mip -of ${oFolder}

    # --------------------------------------------------------------------------
    #  Gist
    # --------------------------------------------------------------------------
    dname=Gist
    n=1000000
    d=960
    dPath=./data/${dname}/${dname}
    oFolder=./results/${dname}/para/c0=${c0}/

    ./alsh -alg 1 -n ${n} -qn ${qn} -d ${d} -c0 ${c0} -c ${c} -ds ${dPath}.ds \
        -qs ${dPath}.q -ts ${dPath}.mip -of ${oFolder}

    # --------------------------------------------------------------------------
    #  Netflix
    # --------------------------------------------------------------------------
    dname=Netflix
    n=17770
    d=300
    dPath=./data/${dname}/${dname}
    oFolder=./results/${dname}/para/c0=${c0}/

    ./alsh -alg 1 -n ${n} -qn ${qn} -d ${d} -c0 ${c0} -c ${c} -ds ${dPath}.ds \
        -qs ${dPath}.q -ts ${dPath}.mip -of ${oFolder}

    # --------------------------------------------------------------------------
    #  Yahoo
    # --------------------------------------------------------------------------
    dname=Yahoo
    n=624961
    d=300
    dPath=./data/${dname}/${dname}
    oFolder=./results/${dname}/para/c0=${c0}/

    ./alsh -alg 1 -n ${n} -qn ${qn} -d ${d} -c0 ${c0} -c ${c} -ds ${dPath}.ds \
        -qs ${dPath}.q -ts ${dPath}.mip -of ${oFolder}
done
