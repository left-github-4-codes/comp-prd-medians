A.(general instructions)
--------------------------------------------------------------------
codes for all tables, figures, examples, and simulations in 

"Exact computation of projection regression depth and fast computation of depth induced median and estimators".
----------------------------------------------------------------------

are downloadable. To download a specific file menthioned below with name (xxx.R or xxx.m)
just type the link: https://www.stt.msu.edu/users/zuo/Codes/2020/file-name 
 

B. (details on codes)
------------------------------
Table1 (used codes):  Ex_UF_2plus_Z19.R, Ex_UF_2plus_ZZ_v1.R, AA_UF_3.m

Figure 1: codes-4-figure1-1.R

Table2: test_cpp-1.cpp and cpp_test-1.R 

#(note when p=2, use all B to form the convex hull to generate beta while p>2, just use the (p+1)
#deepest beta to form the convex hull to generate beta from it) 
(old ones, codes-4-table2-v7.R; codes-4-table2-v6.R)

tuning paramaters: 
RepN =1000 (the total number of replications, it is fixed for tables 2)

get_N_beta function which gets beta's determined by p points
that produce a hyperplane by y=x'\beta, 
it is fixed for table 2 and could be increased, which might be better for T^*_{PRD})

In p=2 case:
ND=300 (the total number of random directions used for UF, it should increase with p)

Nbet=500 (the total number of beta selected in the convex hull formed by (p+1) minimum
UF points, should increase with p, in p=2 case, the convex formed by all RN beta's)

RN=800, (the number of beta in function get_N_beta), it should be at least (p+1) greater than Nbet.

N1=900 (used for determine N below, the minimum value of N)

N=min((n+300)(p-1), choose(n, p), N1) (the total number of beta you try to select from total (n choose p))

The following are the values of all parameters used in Table 2 (of course, they might not be optimal)

p=2, RepN =1000, ND=300, Nbet=500, RN=800, N1=900
p=3, RepN =1000, ND=500, Nbet=700, RN=800, N1=900
p=4, RepN =1000, ND=500, Nbet=700, RN=1000, N1=1000
p=6, RepN =1000, ND=600, Nbet=800, RN=1200, N1=1200


Table 3: cpp_test-1-4-table3.R and test_cpp-1.cpp

Using the same values for the tuning parameters as those in table 2

Figures 2 and 3: cpp_test-1-4-table3-4-boxplot.R and test_cpp-1.cpp
(tuning parameters include p, ND, Nbet, RN, N1, update them accordingly)

Table 4 cpp_test-1-4-table4-v2.R and test_cpp-1.cpp

Using the same values for the tuning parameters as those in table 2.