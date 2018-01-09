make
rm *.o

#-------------------------------------------------------------------------------
#  Usage of the package of internal approximate MIP search
#-------------------------------------------------------------------------------
#    -alg  {integer}  options of algorithms (0 - 13)
#    -n    {integer}  cardinality of the dataset
#    -d    {integer}  dimensionality of the dataset
#    -qn   {integer}  number of queries
#    -K    {integer}  number of hash tables for Sign_ALSH and Simple_LSH
#    -m    {integer}  additional dimension for L2_ALSH, L2_ALSH2, and Sign_ALSH
#    -U    {real}     a real value (0,1] for L2_ALSH, L2_ALSH2, and Sign_ALSH
#    -c    {real}     approximation ratio for nn search (c > 1)
#    -C    {real}     approximation ratio for mip search (0 < C < 1) for H2_ALSH 
#    -ds   {string}   address of the dataset
#    -qs   {string}   address of the query set
#    -ts   {string}   address of the ground truth set
#    -of   {string}   output folder to store output results
#
#-------------------------------------------------------------------------------
#  The options of algorithms are:
#-------------------------------------------------------------------------------
#    0 - Generating ground truth results
#        Parameters: -alg 0 -n -qn -d -ds -qs -ts
#
#    1 - Running AMIP search by L2_ALSH
#        Parameters: -alg 1 -n -qn -d -m -U -c -ds -qs -ts -of
#
#    2 - Running AMIP search by L2_ALSH2
#        Parameters: -alg 2 -n -qn -d -m -U -c -ds -qs -ts -of
#
#    3 - Running AMIP search by XBox and XBox2
#        Parameters: -alg 3 -n -qn -d -c -ds -qs -ts -of
#
#    4 - Running AMIP search by H2_ALSH
#        Parameters: -alg 4 -n -qn -d -c -C -ds -qs -ts -of
#
#    5 - Running AMIP search by Sign_ALSH
#        Parameters: -alg 5 -n -qn -d -K -m -U -c -ds -qs -ts -of
#
#    6 - Running AMIP search by Simple_LSH
#        Parameters: -alg 6 -n -qn -d -K -c -ds -qs -ts -of
#
#    7 - Running MIP search by Linear_Scan
#        Parameters: -alg 7 -n -qn -d -ds -qs -ts -of
#
#    8 - Running Precision-Recall Curve of AMIP search by L2_ALSH
#        Parameters: -alg 1 -n -qn -d -m -U -c -ds -qs -ts -of
#
#    9 - Running Precision-Recall Curve of AMIP search by L2_ALSH2
#        Parameters: -alg 2 -n -qn -d -m -U -c -ds -qs -ts -of
#
#    10 - Running Precision-Recall Curve of AMIP search by XBox and XBox2
#        Parameters: -alg 3 -n -qn -d -c -ds -qs -ts -of
#
#    11 - Running Precision-Recall Curve of AMIP search by H2_ALSH
#        Parameters: -alg 4 -n -qn -d -c -C -ds -qs -ts -of
#
#    12 - Running Precision-Recall Curve of AMIP search by Sign_ALSH
#        Parameters: -alg 5 -n -qn -d -K -m -U -c -ds -qs -ts -of
#
#    13 - Running Precision-Recall Curve of AMIP search by Simple_LSH
#        Parameters: -alg 6 -n -qn -d -K -c -ds -qs -ts -of
#
#-------------------------------------------------------------------------------

dname=Sift
n=1000000
d=128

qn=1000
K=512
m=3
U1=0.83
U2=0.85
c=2.0
C=0.25
dPath=../../data/${dname}/${dname}
oFolder=../../results/${dname}/XBox3_Mem/

./alsh -alg 0 -n ${n} -qn ${qn} -d ${d} -ds ${dPath}.ds -qs ${dPath}.q -ts ${dPath}ID.mip

./alsh -alg 1 -n ${n} -qn ${qn} -d ${d} -m ${m} -U ${U1} -c ${c} -ds ${dPath}.ds -qs ${dPath}.q -ts ${dPath}ID.mip -of ${oFolder}

./alsh -alg 2 -n ${n} -qn ${qn} -d ${d} -m ${m} -U ${U1} -c ${c} -ds ${dPath}.ds -qs ${dPath}.q -ts ${dPath}ID.mip -of ${oFolder}

./alsh -alg 3 -n ${n} -qn ${qn} -d ${d} -c ${c} -ds ${dPath}.ds -qs ${dPath}.q -ts ${dPath}ID.mip -of ${oFolder}

./alsh -alg 4 -n ${n} -qn ${qn} -d ${d} -c ${c} -C ${C} -ds ${dPath}.ds -qs ${dPath}.q -ts ${dPath}ID.mip -of ${oFolder}

./alsh -alg 5 -n ${n} -qn ${qn} -d ${d} -K ${K} -m ${m} -U ${U2} -c ${c} -ds ${dPath}.ds -qs ${dPath}.q -ts ${dPath}ID.mip -of ${oFolder}

./alsh -alg 6 -n ${n} -qn ${qn} -d ${d} -K ${K} -c ${c} -ds ${dPath}.ds -qs ${dPath}.q -ts ${dPath}ID.mip -of ${oFolder}

./alsh -alg 7 -n ${n} -qn ${qn} -d ${d} -ds ${dPath}.ds -qs ${dPath}.q -ts ${dPath}ID.mip -of ${oFolder}

./alsh -alg 8 -n ${n} -qn ${qn} -d ${d} -m ${m} -U ${U1} -c ${c} -ds ${dPath}.ds -qs ${dPath}.q -ts ${dPath}ID.mip -of ${oFolder}

./alsh -alg 9 -n ${n} -qn ${qn} -d ${d} -m ${m} -U ${U1} -c ${c} -ds ${dPath}.ds -qs ${dPath}.q -ts ${dPath}ID.mip -of ${oFolder}

./alsh -alg 10 -n ${n} -qn ${qn} -d ${d} -c ${c} -ds ${dPath}.ds -qs ${dPath}.q -ts ${dPath}ID.mip -of ${oFolder}

./alsh -alg 11 -n ${n} -qn ${qn} -d ${d} -c ${c} -C ${C} -ds ${dPath}.ds -qs ${dPath}.q -ts ${dPath}ID.mip -of ${oFolder}

./alsh -alg 12 -n ${n} -qn ${qn} -d ${d} -K ${K} -m ${m} -U ${U2} -c ${c} -ds ${dPath}.ds -qs ${dPath}.q -ts ${dPath}ID.mip -of ${oFolder}

./alsh -alg 13 -n ${n} -qn ${qn} -d ${d} -K ${K} -c ${c} -ds ${dPath}.ds -qs ${dPath}.q -ts ${dPath}ID.mip -of ${oFolder}
