!!! Gillespie f90 code to simulate 'miRNA-mediated negative feedback loop' model 
!!! Developer : Raunak Adhikary 
!!! Manuscript name : Effects of microRNA-mediated negative feedback on gene expression noise 
!!! Authors : Raunak Adhikary1+, Arnab Roy1+, Mohit Kumar Jolly2*, Dipjyoti Das1++
!!! 1 Department of Biological Sciences, Indian Institute of Science Education And Research Kolkata Mohanpur, Nadia 741246, West Bengal, India.
!!! 2 Centre for BioSystems Science and Engineering, Indian Institute of Science, Bengaluru 560012, India.
!!! ++ Corresponding authors
!!! Email : ra19rs093@iiserkol.ac.in, dipjyoti.das.@iiserkol.ac.in, mkjolly@iisec.ac.in
!!! Date : June 2020
!!! this f90 code generates mRNA, protein, miRNA and mRNA-miRNA complex numbers corresponding to ON and OFF state at any time t


program gillespie
implicit none 


!!! --------------  First we declare the variables. The description of the variables are shown on right hand side. ----------------------- 


real, dimension(13) :: c 				!!! c's are the reaction parameters
real, dimension(11) :: a				!!! a's are the propensity funtions

integer :: mu,nu,seed 					!!! mu, nu are required to choose a reaction 
							!!! seed is required to generate random number from uniform distribution
							   
real :: m,p,s,com,g_c                                   !!! variables m, p, s, com compute number of mRNA, protein, miRNA and 
							!!! mRNA-miRNA complex molecules at any time t. 
							!!! g_c is a binary variable that takes only two values 0 and 1,
							!!! corresponding to miRNA gene OFF and ON states.

real :: mon,moff,pon,poff,son,soff,comon,comoff		!!! variables mon and moff calculates number of mRNA in ON and OFF states, 
							!!! respectively, at any time t. m = moff + moff.
							!!! Similar for other variables. 
							   
real :: T,T2,TINT,TPRINT,tao,steady_t 			!!! T = system evolution time, T2 = stopping time,
							!!! TINT = interval time, TPRINT = time required for printing outputs in datafile,
							!!! steady_t = time considered as starting time of steady state 

real:: summ,R1,R2,R2A0,A0 				!!! standard parameters required for gillespie algorithm

real :: gamma0,alpha					!!! system parameter gamma and alpha
       
real :: run_t1,run_t2					!!! variables required to count total time required for simulation

call cpu_time(run_t1)           


!!!! ------------ assigning initial values and parameter values to the variables -------------------------------------------

T2=250000; TINT=0.5;					!!! stopping time is 250000 second, TINT is 0.5 second

T=0; steady_t = 50000;					!!! starting time of steady state is 50000 second

mon=0; moff=0; m=0; pon=0; poff=0; p=0			!!! initial counts of mRNA, protein, miRNA and mRNA-miRNA complex 
son=0; soff=0; s=0; comon=0; comoff=0; com=0; 

g_c=0;							!!! initial state of miRNA coding gene

gamma0 = 0.0004; alpha = 0.95                           !!! alpha = cataliticity parameter
							!!! gamma0 = mRNA-miRNA complex overall degradation rate


 c(1)= 0.0007                               		!!! c(1) = mRNA transcription rate (k_r)
 c(2)= 0.0004                               		!!! c(2) = mRNA degradation rate (g_r)
 c(3)= 0.4                                  		!!! c(3) = protein translation rate (k_p)
 c(4)= 0.0005                               		!!! c(4) = protein degradation rate (g_p)
 c(5)= 0                                    		!!! c(5) = miRNA gene activation rate (k_act)
 c(6)= 0.01                                 		!!! c(6) = miRNA gene deactivation rate (k_deact)
 c(7)= 0.0005                               		!!! c(7) = miRNA gene basal transcription rate (k_s_0)
 c(8)= 0.5                                  		!!! c(8) = miRNA gene enhanced transcription rate (k_s)
 c(9)= 0.0003                               		!!! c(9) = miRNA degradation rate (g_s)
 c(10)= 0.1                                 		!!! c(10) = mRNA-miRNA association rate (k_+)
 c(11)= 0.1                                 		!!! c(11) = mRNA-miRNA dissociation rate (k_-)
 c(12)= alpha*gamma0                        		!!! c(12) =alpha*gamma
 c(13)= (1-alpha)*gamma0                    		!!! c(13) =(1-alpha)*gamma

!!!! assign value of T to TPRINT --------------------------------------------------------------

TPRINT=T

!!! open a data file and write the values------------------------------------------------------
 
open(2, File= 'parameter.txt',status='unknown', position='append')		!!! open a datafile to print the system parameters 

open(3, File= 'output.dat',status='unknown', position='append')			!!! open a datafile to print molecule numbers at any time

write(2,5) c(1),c(2),c(3),c(4),c(5),c(6),c(7),c(8),c(9),c(10),c(11),c(12),c(13),gamma0,alpha,T2,TINT
5 format("##","k_r=",1X,F12.7,1X,"; g_r=",1X,F12.7,1X,"; k_p=",1X,F12.7,1X,"; g_p=",1X,F12.7,1X, &
"k_act=",1X,F15.9,1X,"; k_deact=",1X,F15.9,1X,"; k_s_0=",1X,F12.7,1X,"; k_s=",1X,F12.7,1X, &
"; gamma_s=",1X,F12.7,1X,"k_+ =",1X,F12.7,1X,"; k_- =",1X,F12.7,1X,"; alpha*gamma0 =",1X,F15.9,1X, &
"; (1-alpha)*gamma0 =",1X,F15.9,"; gamma0 =",1X,F12.7,1X,"; alpha =",1X,F12.7,1X,'## T2 = ',1X,F12.4,1X, &
'TINT = ',1X,F12.7)
flush(2)

write(2,*)'## (1)TPRINT   (2)mon   (3)moff   (4)m   (5)pon   (6)poff   (7)p   (8)son   (9)soff   (10)s   &
		(11)comon   (12)comoff   (13)com    (14)gon'
flush(2)

!!!! call random number generator -------------------------------------------------------------------

   call SYSTEM_CLOCK(COUNT=seed)          	!!! SYSTEM_CLOCK is a inbuild fortran subroutine to call a random number generator 
                                          	!!!    and      'count=seed'
                                          	!!! gives the seed for the random number generator srand
   call srand(seed)


!!!! a(i)'s are the propensity functions. Here we calculate the values of a(i) ------------------------------------------

10 a(1) = c(1)                              	!!! a(1) = mRNA transcription rate (k_r)
   a(2) = c(2)*m                            	!!! a(2) = mRNA degradation rate * mRNA number (g_r * m)
   a(3) = c(3)*m                            	!!! a(3) = protein transcription rate * mRNA number (k_p * m)
   a(4) = c(4)*p                            	!!! a(4) = protein degradation rate * protein number (g_p * p)


   if (g_c.EQ.0) then 				!!! if miRNA gene is OFF	
           a(5) = c(5)*p                    	!!! a(5) = miRNA gene activation rate * protein number (k_act * p)
   else                                     	!!! else if miRNA gene is ON
           a(5) = c(6)                      	!!! a(6) = miRNA gene deactivation rate (k_deact)
           
   end if 

   if (g_c.EQ.0) then				!!! if miRNA gene is OFF
           a(6) = c(7)                      	!!! a(6) = miRNA gene basal transcription rate (k_s_0)
   else                                     	!!! else if miRNA gene is ON
           a(6) = c(8)                      	!!! a(6) = miRNA gene enhanced transcription rate (k_s)
           
   end if 
   
   a(7) = c(9)*s                            	!!! a(7) = miRNA degradation rate * miRNA number (g_s * s)
   a(8) = c(10)*m*s                         	!!! a(8) = mRNA-miRNA association rate * mRNA number * miRNA number (k_+ * m * s)
   a(9) = c(11)*com                         	!!! a(9) = mRNA-miRNA dissociation rate * complex number (k_- * c)
   a(10) = c(12)*com                        	!!! a(10) = alpha * gamma * complex number (alpha * gamma * c)
   a(11) = c(13)*com                        	!!! a(11) = gamma * (1 - alpha) * complex number (gamma * (1 - alpha) * c)


   a0 = a(1)+a(2)+a(3)+a(4)+a(5)+a(6)+a(7)+a(8)+a(9)+a(10)+a(11)    		!!! total propensity 
 
!!!! calculate R1 & R2 by random number generator ---------------------------------------------------

20  R1 = rand(0)
    R2 = rand(0)
    
    if (R1.LT.1E-30.OR.R2.LT.1E-30) then	!!! if one of R1 or R2 is less than 10^(-30) then calculate R1 and R2 again
       go to 20
    end if
    
!!!! increment T by tao ------------------------------------------------------------------------
   
   tao = alog(1/R1)/a0

21 T = T + alog(1/R1)/a0                     	!!! increase T by tao
   
22 if (T.LT.TPRINT) then
   go to 25
   end if
   
!!!! printing the molecule numbers -------------------------------------------------------------
!!!! Note that, we only print the steady state molecule numbers --------------------------------

23 if (TPRINT.GE.steady_t) then

	write(3,9) TPRINT,mon,moff,m,pon,poff,p,son,soff,s,comon,comoff,com,g_c
9 	format(F20.10,1X,F30.10,1X,F30.10,1X,F30.10,1X,F30.10,1X,F30.10,1X,F30.10,1X,F30.10,1X,F30.10,1X,F30.10,1X, &
	F30.10,1X,F30.10,1X,F30.10,1X,F30.10)
	flush(3)

   end if

!!! increase TPRINT ----------------------------------------------------------------------

 if (TPRINT.GE.T2) go to 46
  
  TPRINT=TPRINT + TINT
  
  go to 22

!!! calculate mu ------------------------------------------------------------------ 

25 R2A0=R2*A0
   summ=0

26 do nu=1,11
      mu=nu
      summ=summ+ a(nu)
      if(summ.GE.R2A0) go to 30
   end do

!!! find the reaction to happen and update molecule numbers according to reaction ---------------------

30 go to(31,32,33,34,35,36,37,38,39,40,41), mu

31 m=m+1                                          !!! reaction: unconstituted gene -> free mRNA (k_r) ------------
   
   if (g_c.EQ.0) then
        moff = m; poff = p; soff = s; comoff = com
        mon = 0; pon = 0; son = 0; comon = 0
   end if
   if (g_c.EQ.1) then
        mon = m; pon = p; son = s; comon = com
        moff = 0; poff = 0; soff = 0; comoff = 0
   end if
   
   go to 45

32 m=m-1                                          	!!!reaction: free mRNA -> degraded mRNA (g_r)  --------------
   
   if (g_c.EQ.0) then
        moff = m; poff = p; soff = s; comoff = com
        mon = 0; pon = 0; son = 0; comon = 0
   end if
   if (g_c.EQ.1) then
        mon = m; pon = p; son = s; comon = com
        moff = 0; poff = 0; soff = 0; comoff = 0
   end if
   
   go to 45

33 p=p+1                                          	!!!reaction: free mRNA -> free protein (k_p)  -----------------
   
   if (g_c.EQ.0) then
        moff = m; poff = p; soff = s; comoff = com
        mon = 0; pon = 0; son = 0; comon = 0
   end if
   if (g_c.EQ.1) then
        mon = m; pon = p; son = s; comon = com
        moff = 0; poff = 0; soff = 0; comoff = 0
   end if
  
   go to 45

34 p=p-1                                          	!!!reaction: free protein -> degraded protein (g_p) -------------
   
   if (g_c.EQ.0) then
        moff = m; poff = p; soff = s; comoff = com
        mon = 0; pon = 0; son = 0; comon = 0
   end if
   if (g_c.EQ.1) then
        mon = m; pon = p; son = s; comon = com
        moff = 0; poff = 0; soff = 0; comoff = 0
   end if
   
   go to 45

35 if (g_c.EQ.0) then
   
   p=p-1
   g_c=g_c+1				!!! reaction: deactivated miRNA-gene + free protein -> activated miRNA-gene (k_act) ---------
   
   mon = m; pon = p; son = s; comon = com;            
   moff = 0; poff = 0; soff = 0; comoff = 0 
   
   go to 45

   else
 
   p=p+1                             	!!! reaction: activated miRNA-gene -> deactivated miRNA-gene + free protein (k_deact) -------
   g_c=g_c-1

   moff = m; poff = p; soff = s; comoff = com
   mon = 0; pon = 0; son = 0; comon = 0
   
   go to 45

   end if
 
36 if (g_c.EQ.0) then

   s=s+1                                    		!!! reaction: deactivate miRNA-gene -> free miRNA (k_s_0) ---------
   
   moff = m; poff = p; soff = s; comoff = com
   mon = 0; pon = 0; son = 0; comon = 0
   
   go to 45

   else
 
   s=s+1                                     		!!! reaction: activated miRNA-gene -> free miRNA (k_s) -----------
 
   mon = m; pon = p; son = s; comon = com
   moff = 0; poff = 0; soff = 0; comoff = 0
   
   go to 45 
 
   end if

37 s=s-1                                              	!!! reaction: free miRNA -> degraded miRNA (g_s)  --------------
   
   if (g_c.EQ.0) then
        moff = m; poff = p; soff = s; comoff = com
        mon = 0; pon = 0; son = 0; comon = 0
   end if
   if (g_c.EQ.1) then
        mon = m; pon = p; son = s; comon = com
        moff = 0; poff = 0; soff = 0; comoff = 0
   end if
   
   go to 45

38 m=m-1
   s=s-1                                              	!!! reaction: free mRNA + free miRNA -> mRNA-miRNA complex (k_+)  --------
   com=com+1                                              
   
   if (g_c.EQ.0) then
        moff = m; poff = p; soff = s; comoff = com
        mon = 0; pon = 0; son = 0; comon = 0
   end if
   if (g_c.EQ.1) then
        mon = m; pon = p; son = s; comon = com
        moff = 0; poff = 0; soff = 0; comoff = 0
   end if
   
   go to 45

39 m=m+1                                              	!!! reaction: mRNA-miRNA complex -> free mRNA + free miRNA (k_-) --------
   s=s+1
   com=com-1
   
   if (g_c.EQ.0) then
        moff = m; poff = p; soff = s; comoff = com
        mon = 0; pon = 0; son = 0; comon = 0
   end if
   if (g_c.EQ.1) then
        mon = m; pon = p; son = s; comon = com
        moff = 0; poff = 0; soff = 0; comoff = 0
   end if
   
   go to 45

40 com=com-1                                         	!!! reaction: mRNA-miRNA complex -> degraded mRNA-miRNA complex (alpha * gamma)
  
   if (g_c.EQ.0) then
        moff = m; poff = p; soff = s; comoff = com
        mon = 0; pon = 0; son = 0; comon = 0
   end if
   if (g_c.EQ.1) then
        mon = m; pon = p; son = s; comon = com
        moff = 0; poff = 0; soff = 0; comoff = 0
   end if
   
   go to 45

41 com=com-1                                          	!!! reaction: mRNA-miRNA complex -> free miRNA (gamma * (1-alpha)) -------
   s=s+1
   
   if (g_c.EQ.0) then
        moff = m; poff = p; soff = s; comoff = com
        mon = 0; pon = 0; son = 0; comon = 0
   end if
   if (g_c.EQ.1) then
        mon = m; pon = p; son = s; comon = com
        moff = 0; poff = 0; soff = 0; comoff = 0
   end if
   
   go to 45

45 if (T.LT.T2) then 					!!! check if system time exceeds stopping time 
   go to 10                           			!!! if system time < stopping time then, loop continues 
   end if

46 call cpu_time(run_t2)
write(2,*) '# total time taken =',(run_t2 - run_t1)	!!! print total time 
write(2,*)

end program gillespie					!!! end of the program

