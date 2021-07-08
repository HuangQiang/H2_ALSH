#!/bin/bash
make clean
make -j

qn=1000
dataPath=../data
resultPath=../results

# ------------------------------------------------------------------------------
#  Tuning approximation ratio c for c-AMIP search
# ------------------------------------------------------------------------------
c0=2.0
for c in 0.1 0.25 0.5 0.8
do
    # --------------------------------------------------------------------------
    #  Sift
    # --------------------------------------------------------------------------
    dname=Sift
    n=1000000
    d=128
    dPath=${dataPath}/${dname}/${dname}
    oPath=${resultPath}/${dname}/para/c=${c}/

    ./alsh -alg 1 -n ${n} -qn ${qn} -d ${d} -c0 ${c0} -c ${c} -ds ${dPath}.ds \
        -qs ${dPath}.q -ts ${dPath}.mip -op ${oPath}

    # --------------------------------------------------------------------------
    #  Gist
    # --------------------------------------------------------------------------
    dname=Gist
    n=1000000
    d=960
    dPath=${dataPath}/${dname}/${dname}
    oPath=${resultPath}/${dname}/para/c=${c}/

    ./alsh -alg 1 -n ${n} -qn ${qn} -d ${d} -c0 ${c0} -c ${c} -ds ${dPath}.ds \
        -qs ${dPath}.q -ts ${dPath}.mip -op ${oPath}

    # --------------------------------------------------------------------------
    #  Netflix
    # --------------------------------------------------------------------------
    dname=Netflix
    n=17770
    d=300
    dPath=${dataPath}/${dname}/${dname}
    oPath=${resultPath}/${dname}/para/c=${c}/

    ./alsh -alg 1 -n ${n} -qn ${qn} -d ${d} -c0 ${c0} -c ${c} -ds ${dPath}.ds \
        -qs ${dPath}.q -ts ${dPath}.mip -op ${oPath}

    # --------------------------------------------------------------------------
    #  Yahoo
    # --------------------------------------------------------------------------
    dname=Yahoo
    n=624961
    d=300
    dPath=${dataPath}/${dname}/${dname}
    oPath=${resultPath}/${dname}/para/c=${c}/

    ./alsh -alg 1 -n ${n} -qn ${qn} -d ${d} -c0 ${c0} -c ${c} -ds ${dPath}.ds \
        -qs ${dPath}.q -ts ${dPath}.mip -op ${oPath}
done

# ------------------------------------------------------------------------------
#  Tuning approximation ratio c0 for c0-ANN search
# ------------------------------------------------------------------------------
c=0.5
for c0 in 1.5 2.0 2.5 3.0
do
    # --------------------------------------------------------------------------
    #  Sift
    # --------------------------------------------------------------------------
    dname=Sift
    n=1000000
    d=128
    dPath=${dataPath}/${dname}/${dname}
    oPath=${resultPath}/${dname}/para/c0=${c0}/

    ./alsh -alg 1 -n ${n} -qn ${qn} -d ${d} -c0 ${c0} -c ${c} -ds ${dPath}.ds \
        -qs ${dPath}.q -ts ${dPath}.mip -op ${oPath}

    # --------------------------------------------------------------------------
    #  Gist
    # --------------------------------------------------------------------------
    dname=Gist
    n=1000000
    d=960
    dPath=${dataPath}/${dname}/${dname}
    oPath=${resultPath}/${dname}/para/c0=${c0}/

    ./alsh -alg 1 -n ${n} -qn ${qn} -d ${d} -c0 ${c0} -c ${c} -ds ${dPath}.ds \
        -qs ${dPath}.q -ts ${dPath}.mip -op ${oPath}

    # --------------------------------------------------------------------------
    #  Netflix
    # --------------------------------------------------------------------------
    dname=Netflix
    n=17770
    d=300
    dPath=${dataPath}/${dname}/${dname}
    oPath=${resultPath}/${dname}/para/c0=${c0}/

    ./alsh -alg 1 -n ${n} -qn ${qn} -d ${d} -c0 ${c0} -c ${c} -ds ${dPath}.ds \
        -qs ${dPath}.q -ts ${dPath}.mip -op ${oPath}

    # --------------------------------------------------------------------------
    #  Yahoo
    # --------------------------------------------------------------------------
    dname=Yahoo
    n=624961
    d=300
    dPath=${dataPath}/${dname}/${dname}
    oPath=${resultPath}/${dname}/para/c0=${c0}/

    ./alsh -alg 1 -n ${n} -qn ${qn} -d ${d} -c0 ${c0} -c ${c} -ds ${dPath}.ds \
        -qs ${dPath}.q -ts ${dPath}.mip -op ${oPath}
done
