A.(general instructions)
--------------------------------------------------------------------
codes of all tables, figures, examples, and simulations in 

"Computation of Projection Regression depth and its induced median".
----------------------------------------------------------------------

are downloadable. To download a specific file menthioned below with name (xxx.R or xxx.m)
just type the link: https://www.stt.msu.edu/users/zuo/Codes/file-name 
 

B. (details on codes)
---------------------------------
For Table 3 in Section 3.2.
---------------------------------

The main codes for it is codes-4-table-3.7.R (or codes-4-table-3.6.R)

To run the code, first run all function part, then the main body part at the top.
The code includes

tuning paramaters: 
RepN =1000 (the total number of replications, it is fixed for table 3)

Nbet1 =300 (the total number of extra beta's used after 
get_N_beta function which gets beta's determined by p points
that produce a hyperplane by y=x'\beta, it is fixed for table 3 and
could be increased, which might be better for T^*_{PRD})

ND=300 (the total number of random directions used for UF, it should increase with p)

Nbet=400 (the total number of beta selected in the convex hull formed by (p+1) minimum
UF points, should increase with p)

RN=800, (the number of beta in function get_N_beta), it should be at least (p+1) greater than Nbet.

N1=900 (used for determine N below, the minimum value of N)

N=min((n+300)(p-1), choose(n, p), N1) (the total number of beta you try to select from total (n choose p))

The following are the values of all parameters used in Table 3 (of course, they might not be optimal)

p=2, RepN =1000, Nbet1 =300, ND=300, Nbet=500, RN=800, N1=900
p=3, RepN =1000, Nbet1 =300, ND=500, Nbet=700, RN=800, N1=900
p=4, RepN =1000, Nbet1 =300, ND=500, Nbet=700, RN=1000, N1=1000
p=6, RepN =1000, Nbet1 =300, ND=600, Nbet=800, RN=1200, N1=1200

C. (details on codes)
---------------------------------
For Tables 1 and 4
---------------------------------

The main codes for these two tables are matlab codes.

They are:

Ex_UF_2plus_2.m 
(for exact computation of UF in 2D or higher)

Two base building functions are

ufv.m
update_UF.m 


main_4_AAUF_v1.m (for the Table 2, simulation of three AA's)

AA_UF_1.m
AA_UF_2.m
AA_UF_3.m

---------------------------------
For Tables 2 and 3
---------------------------------

table 2: 
Ex_UF_HD_no_UN.m 
and 
AA_UF_1.m
AA_UF_2.m
AA_UF_3.m

table 3:
main_4_AA_VS_Ex.m

--------------------------
For tables 5 and 6
--------------------------
table 5:

test_cpp-1.cpp and cpp_test-1.R 

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


Table 6: cpp_test-1-4-table3.R and test_cpp-1.cpp

Using the same values for the tuning parameters as those in table 5

Table 7:  cpp_test-1-4-table4-v2.R and test_cpp-1.cpp

Using the same values for the tuning parameters as those in table 5.


D.(details on code)
----------------------------------
For Figure 1
----------------------------------
codes-4-figure1.R or codes-4-fig1.txt

----------------------
For Figure 2
---------------------------
circular-sequence-changes.R

angular-region-1.R

----------------------------------
For Figure 4 
-----------------------------------
the matlab codes for the Figure is:
code_4_fig_4.m
-----------------------------
or txt file
-----------------------------
codes-4-fig-4-9-22-19.txt


-----------------------------
For Figure 5, example 3-1
-----------------------------
codes-4-example-3-1-added.R

-------------------------
For Figure 6, example 3-2
-------------------------
codes-4-example-3-2.R

------------------------------
For Figures 7, 8, example 3.3.3
--------------------------------
cpp_test-1-4-table3-4-boxplot.R and test_cpp-1.cpp
(tuning parameters include p, ND, Nbet, RN, N1, update them accordingly)


E.(details on codes)
----------------------------
For Table 4 and 5
----------------------------
efficiency-simulation.v3.R (or ....v4.R)
or
efficiency-optim-optimx.v1.speed-up.R
##################################################
many functions in above codes are limited to p=2,
general versions are given in the codes-4-table-3.7.R
#####################################################