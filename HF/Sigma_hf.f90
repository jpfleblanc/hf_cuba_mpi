      PROGRAM DISP_HF

! Hamiltonian parameters
      DOUBLE PRECISION :: temp, beta, beta2      ! temperature, inv. temperature
      DOUBLE PRECISION :: amu                    ! chemical potential 
      DOUBLE PRECISION :: mass                   ! mass of particle 
      DOUBLE PRECISION :: k_F, E_F               ! Fermi momentum and Energy
      DOUBLE PRECISION :: R_s, E_Ha              ! r_s and Hartree energy = e^2/a_B = m*e^4
      DOUBLE PRECISION :: density                ! total density n_{\up}+n_{down}=(k_F)^3/(3*pi^2) 
      DOUBLE PRECISION :: e2                     ! e^2      
      DOUBLE PRECISION :: Co_Co                  ! 4*pi*e^2  (V(q)=4*pi*e^2/(q^2+kappa^2))  
      DOUBLE PRECISION :: plasma_F               ! plasma frequency
      DOUBLE PRECISION :: kappa                  ! Screening parameter: V(q)=4*pi*e^2/(q^2+kappa^2)
      DOUBLE PRECISION :: kc, kqc                ! cut-offs for G0 and V0/Pi 
      
! Momentum grig
      INTEGER          :: powerp                 ! number of p-points 
      INTEGER          :: Np, Kp1, Kp2, Kp3, Np1 ! Np is number of p-points in G
      INTEGER          :: Nq, Kq1, Kq2, Kq3, Nq1 ! Nq is number of q-points in Pi      
      INTEGER          :: ipow, ipow1, ipow2, ipow3      
      INTEGER          :: ik_F                   ! index corresponding to p=k_F
      DOUBLE PRECISION :: dimp                   ! step of original equidist momentum grid 
      INTEGER, ALLOCATABLE          :: ptoi(:)   ! see explanations below
      DOUBLE PRECISION, ALLOCATABLE :: pgrid(:)  ! momentum grid for G     
      INTEGER, ALLOCATABLE          :: qtoi(:)   ! see explanations below
      DOUBLE PRECISION, ALLOCATABLE :: qgrid(:)  ! momentum grid for Pi (auxiliary here)        

! Constants
      DOUBLE PRECISION, PARAMETER :: pi=3.141592653589793d0 
      DOUBLE PRECISION, PARAMETER :: eilera=0.5772156649d0 
      DOUBLE PRECISION, PARAMETER :: katalana=0.9159655942d0 
      DOUBLE PRECISION, PARAMETER :: e_number=2.718281828459d0 
      DOUBLE PRECISION, PARAMETER :: pi_d2=pi/2.d0 
      DOUBLE PRECISION, PARAMETER :: pi_m2=pi*2.d0 
      DOUBLE PRECISION, PARAMETER :: pi_m4=pi*4.d0

! Adjustment of chemical potential 
      LOGICAL          :: adj_mu                          ! adjust mu if .TRUE.   
      INTEGER          :: flagmu1, flagmu2, flagmu3   
      DOUBLE PRECISION :: ntarget, dmu, acc_mu            ! desired density, initial step, accuracy
      DOUBLE PRECISION :: mun1, mun2, nmu1, nmu2, nmu3    ! service variables for adjusting mu
      DOUBLE PRECISION :: dnmu                            ! service variable 

! Various
      INTEGER :: ipint, cou_obs  
      INTEGER :: NdevidK, NdevidNq
      DOUBLE PRECISION :: amomc, pp, dp, dps, delint, E_k
      DOUBLE PRECISION :: kk, qq, argu, argd, funq, step  
      DOUBLE PRECISION :: bolv1, bolv2, upli, dnli, sinte
      DOUBLE PRECISION :: sglazh, coudump 
      DOUBLE PRECISION :: fund(32769)                     ! function for Simpson integral
      DOUBLE PRECISION :: mumc(0:100000)                  ! chem pot vs MC steps
      DOUBLE PRECISION :: denc(0:100000)                  ! density vs MC
      LOGICAL :: flag_adj

! Sigma_HF
      DOUBLE PRECISION, ALLOCATABLE :: Sigma_HF0(:)
      DOUBLE PRECISION, ALLOCATABLE :: Sigma_HF(:)        ! Hartree-Fock Sigma(p) (Hartree=0 here)
      DOUBLE PRECISION, ALLOCATABLE :: SigHF_IN(:)        ! same as above, but for interpolations
      
! Indices
      INTEGER :: i, j, ik, ip, iprn

! Prints
      CHARACTER :: cxr
      
      cxr=','

!     ======================================================
!     Input
!     ======================================================

      open(1,file='input_hf.dat')
      READ(1,*) temp          ; print*,' temp    =', temp 
      READ(1,*) R_s           ; print*,' R_s     =', R_s
      READ(1,*) kappa         ; print*,' kappa   =', kappa     
      READ(1,*) powerp, Kp1, Kq1, Kp3                       ! momentum-grid parameters 
                print*, ' p-grid #=', powerp, Kp1, Kq1, Kp3    
      close(1)

      beta=1.d0/temp; beta2=beta/2.d0

!     ======================================================
!     Units/Parameters (\hbar=1)
!     ======================================================

      k_F=1.d0                                              ! Fermi momentum
      mass=0.5d0                                            ! mass
      E_F=0.5d0*k_F*k_F/mass                                ! Fermi energy = 1 now
      density=(k_F**3)/(3.d0*pi*pi)                         ! total density for both spins
      amu=E_F                                               ! Initial chemical potential

      e2=(R_s*k_F/mass)*((4.d0/(9.d0*pi))**(1.d0/3.d0))     ! e^2, depends on R_s
      E_Ha=mass*e2*e2                                       ! Hartree energy = e^2/a_B; a_B=1/m*e^2; \hbar=1 
      Co_Co=4.d0*pi*e2                                      ! V(q)=4*pi*e^2/q^2
      plasma_F=dsqrt(4.d0*pi*e2*density/mass)               ! plasma frequency
      q_vol=pi_m2*pi_m2*pi_m2                               ! (2\pi)^3

      print*,'          '
      print*,' k_F     =', k_F      
      print*,' E_F     =', E_F
      print*,' mass    =', mass      
      print*,' e^2     =', e2
      print*,' E_Ha    =', E_Ha
      print*,' plasma_F=', plasma_F
      print*,' q-volume=', q_vol
      print*,' mu      =', amu 
      
!     ======================================================
!     momentum grid
!     ======================================================

      ipow3=2**powerp                                       ! powerp is # of point in the orig q,p-grid till k_F
      dimp=k_F/dfloat(ipow3)                                ! step of original momentum grid 

      IF(dimp*k_F/mass > 0.3333d0*temp) stop 'dp*v_F > T/3'

      Np1=2**Kp1                                            ! # of small steps in each p interval 
      Kp2=powerp-Kp1                                        ! # of p-intervals till k_F
      Nq1=2**Kq1                                            ! # of small steps in each q interval 
      Kq2=powerp-Kq1                                        ! # of q-intervals till k_F
      
      Kq3=Kq2+Kp3+1; Kp3=Kp2+Kp3                          ! # of q/p-intervals between k_F and end of file

      Np=Np1*(Kp2+Kp3)+1; Nq=Nq1*(Kq2+Kq2+Kq3)+1            ! # of p-points (ip=0,Np) and # of q-points (iq=0,Nq)

      print*, ' '
      print*, ' mom  grid ',' dimp=',dimp,' Kq2=',Kq2,' Kq3=',Kq3
      print*, ' Np =',Np,' q-grid is irrelevant here' 
      print*, ' '

      ipow=powerp+(Kp3-Kp2) 
      ipow=(2**powerp)+(2**ipow)-Np1                        ! dimp*ipow is the definition of kc (see in setgrids)
      print*, ' ptoi', ipow, dimp*ipow 

!     Array ptoi() gives the closest to the given momentum p (in the equidistant grid)  
!     momentum p* in the pgrid() array (k_p=p/dimp; j=ptoi(k_p); p*=pgrid(j)). It is 
!     defined from 0 to max index ipow, corresponding to kc/dimp, where kc is the U-V 
!     cut-off for momentum p (for Greens function). See FUNCTION SIGMAF(amom) below 
!     (interpolates Sigma for arbitrary momentum from 0 to max_p=pgrid(Np)) on how to 
!     practically use ptoi().

      allocate(ptoi(0:ipow)) 
      
      ipow1=powerp+(Kq3-Kp2) 
      ipow1=(2**powerp)+(2**powerp)+(2**ipow1)-Np1          ! dimp*ipow1 is the definition of kqc (see in setgrids)
      print*, ' qtoi', ipow1, dimp*ipow1

!     Array qtoi() has the same meaning as ptoi(). It is for Pi/W and auxiliary here.

      allocate(qtoi(0:ipow1))      

      CALL SETGRIDS                                         ! inside pgrid() and qgrid() are allocated

      allocate(Sigma_HF0(0:Np))                             ! for averaging
      allocate(Sigma_HF(0:Np))                              ! answer, Sigma_HF(p) on pgrid()
      allocate(SigHF_IN(0:Np))                              ! Sigma array to use in interpolations

      open(1,file='pgrid.dat')
      write(1,*) Np
      do ip=0,Np
       write(1,*) ip, pgrid(ip)
      enddo
      close(1)
      
      open(1,file='ptoi.dat')
      write(1,*) ipow
      write(1,*) dimp 
      do ip=0,ipow
       write(1,*) ip, ptoi(ip)
      enddo     
      close(1)
      
      print*,' ik_F=', ik_F
      
!     ======================================================
!     Density adjustment parameters
!     ======================================================

      mun1=-100.d0; nmu1=0.d0; mun2=-100.d0; nmu2=0.d0; dnmu=-100.d0

      ntarget=density                                       ! target density = (k_F**3)/(3.d0*pi*pi)

      print*,'         '
      print*,' ntarget        =', ntarget

      adj_mu=.TRUE.                                         ! do mu adjustment if .true.
      dmu=0.1d0                                             ! initial step for adjustment
      acc_mu=1.d-14                                         ! accuracy of adjustment

      mumc=0.d0                                             ! mu(MC step) array
      denc=0.d0                                             ! density(MC step) array

!     ======================================================
!     Initial Sigma for iterations
!     ======================================================

      Sigma_HF=0.d0; SigHF_IN=0.d0

!     ======================================================
!     Parameters for integrations
!     ======================================================

      NdevidK=8192                                          ! to use in calculations of Sigma
      NdevidNq=8192                                        ! to use in calculations of density

!     ======================================================
!     Initial density (should be the same as "density" above)
!     ======================================================

      print*,' Initial density=', N_R0(NdevidNq)            ! NdevidNq points for each p-step in pgrid
      print*,'                 '

!     ======================================================
!     Parameters for self-consisency loop
!     ======================================================

      cou_obs=0

      flag_adj=.FALSE.
      
      coudump=0.d0 
      
      step=0.d0

      mumc(cou_obs)=amu; denc(cou_obs)=N_R0(NdevidNq)
      
      print*,' step=',cou_obs,'mu=',amu,'dens=',denc(cou_obs)

!     ======================================================
      DO                                                    ! Start self-consistency loop
!     ======================================================

      step=step+1.d0
      
      Sigma_HF0=0.d0                                         ! Sigma Hartree Fock (ik)

      DO ik=0,Np
       kk=pgrid(ik)                                         ! External momentum
                                     
       bolint=0.d0
       DO ip=0,Np-1                                         ! Integral from 0 to "infinity"
        pp=pgrid(ip)                                        ! internal momentum from pgrid for integration
        dp=pgrid(ip+1)-pp                                   ! current interval in a given scale of p-grid 
        dps=dp/NdevidK                                      ! (will be devided into smaller steps dps)

        fund=0.d0
        DO ipint=1,NdevidK+1                                ! take SIMPS integral over p within each small step
         qq=pp+(ipint-1)*dps
!        =========================================================
!        Open the line below to compute everything numerically
!        =========================================================
!!!         funq=N_q(qq)                                    
!        =========================================================
!        Close the block below to compute everything numerically
!        =========================================================
         if(qq .le. k_F) then
          funq=N_q(qq)-1.d0; if((pp+dps)>k_F) funq=N_q(qq)
         else 
          funq=N_q(qq)
         endif      
!        =========================================================         
         argu=(kk+qq)*(kk+qq)+kappa*kappa
         argd=(kk-qq)*(kk-qq)+kappa*kappa            
         fund(ipint)=qq*funq*DLOG(argu/argd)  
        ENDDO
        upli=pp; dnli=pp+dp
        CALL SIMPINT(upli,dnli,NdevidK,fund,sinte)          ! Int dqq F(qq,kk)

        bolint=bolint+sinte
       ENDDO                                                ! internal momentum cycle (q)

       bolv1=-bolint*e2/(2.d0*pi*kk)     

!      =========================================================
!      Close the line below to compute everything numerically
!      =========================================================  

       bolv2=MAIN_INT(kk)                                   ! Computes the integral at T=0 analytically

!      =========================================================
!      Open the line below to compute everything numerically
!      =========================================================

!!!       bolv2=0.d0

!      =========================================================

       Sigma_HF0(ik)=bolv1+bolv2
      ENDDO                                                 ! external momentum cycle (q)

      coudump=coudump+1.d0
      DO ik=0,Np
       sglazh=Sigma_HF(ik) 
       sglazh=(sglazh*(coudump-1.d0)+Sigma_HF0(ik))/coudump
       Sigma_HF(ik)=sglazh
      ENDDO
      if(.NOT.flag_adj) coudump=0 
      
      SigHF_IN=Sigma_HF                                     ! array for interpolating in SIGMAF(k)

!     ======================================================
!     Print Sigma_HF and full Hartre-Fock Energy
!     ======================================================

      cou_obs=cou_obs+1
      mumc(cou_obs)=amu; denc(cou_obs)=N_R0(NdevidNq)

      open(1,file='All_mu.dat')
      do iprn=0,cou_obs
       write(1,*) iprn, mumc(iprn), denc(iprn)
      enddo
      close(1)
      
      print*,' step=',cou_obs,'mu=',mumc(cou_obs),'dens=',denc(cou_obs)

      open(1,file='Sigma_HF.dat')
      write(1,*) amu, N_R0(NdevidNq), density 
      do ik=0,Np; kk=pgrid(ik) 
       write(1,*) kk, Sigma_HF(ik)
      enddo
      close(1)

!!!      open(1,file='Sigma_HF-Kun.dat')
!!!      do ik=0,Np; kk=pgrid(ik) 
!!!       write(1,*) kk, (Sigma_HF(ik)-Sigma_HF(ik_F))
!!!      enddo
!!!      close(1)

      open(1,file='Energy_grid.dat')      
      do ik=0,Np; kk=pgrid(ik) 
       E_k=0.5d0*(kk*kk/mass)+Sigma_HF(ik)-amu              ! Full Hartree-Fock dispersion
       write(1,*) kk, E_k
      enddo
      close(1)

!     =======================================================
!     mu adjustment for density - it has to be ntarget = 
!     (k_F^3)/(3.d0*pi^2) with nonzero Sigma_HF as well 
!     (Lattinger Th).
!     =======================================================
      
      IF(adj_mu) THEN

       mun1=amu; nmu1=N_R0(NdevidNq)
       IF(DABS(nmu1-ntarget)<acc_mu) go to 55
       flagmu1=-1; IF(nmu1>ntarget) flagmu1=1

       do
        mun2=mun1+dmu; IF(flagmu1==1) mun2=mun1-dmu  
        amu=mun2; nmu2=N_R0(NdevidNq)       
        flagmu2=-1; IF(nmu2>ntarget) flagmu2=1 
        IF(flagmu1*flagmu2==-1) exit                        ! target crossed
        dmu=dmu*2.d0                                        ! not crosssed yet
       enddo
         
       do
        amu=(mun1+mun2)/2.d0; nmu3=N_R0(NdevidNq)        
        flagmu3=-1; IF(nmu3>ntarget) flagmu3=1 
        IF(DABS(nmu3-ntarget)<acc_mu) exit                  ! on target
        IF(flagmu1*flagmu3==1) then
         mun1=amu
        else
         mun2=amu
        endif
       enddo
        
       dmu=DABS(mun2-mun1)

 55    flag_adj=.TRUE.  

      ENDIF

!     ======================================================
!     We adjusted chemical pottential. Lets repeat calcula-  
!     tions of Sigma_HF, and so on - self-consistency loop.
!     ======================================================
      
      if(step>1.d10) EXIT

!     ======================================================
      ENDDO                                                 ! Self-consistency loop
!     ======================================================

 10   format(e19.12,a1,10(e19.12,a1),e19.12)      



      CONTAINS



!     ======================================================
!     Functions and Subroutines used in the main code:
!     ======================================================


!
!     ======================================================
!     Occupation number n_q
!     ======================================================
!
      DOUBLE PRECISION FUNCTION N_q(qq)
      IMPLICIT NONE  
      DOUBLE PRECISION, INTENT (IN) :: qq 
      DOUBLE PRECISION :: E_q, occ

      N_q=0.d0

      E_q=0.5d0*(qq*qq/mass)+SIGMAF(qq)-amu
      IF(E_q<0.d0) THEN
       occ=1.d0/(1.d0+dexp(E_q*beta))                       ! occ number            
      ELSE		        
       occ=dexp(-E_q*beta)/(1.d0+dexp(-E_q*beta))                         
      ENDIF 

      N_q=occ

      END FUNCTION N_q

!
!     ======================================================
!     Occupation number n(r=0) obtained by integrating n_q 
!     over {\bf q}. This is to compare with "density" when 
!     performing the self-consistency loop for Sigma_HF.
!     ======================================================
!
      DOUBLE PRECISION FUNCTION N_R0(NpoNq)
      IMPLICIT NONE 
      INTEGER, INTENT (IN) :: NpoNq     
      INTEGER :: ik, ipin 
      DOUBLE PRECISION :: kk, dk, dks, qcur, bolin 
      DOUBLE PRECISION :: upli, dnli, sitec  
      DOUBLE PRECISION :: fun_nq(32769)   

      N_R0=0.d0
      
      bolin=0.d0
      DO ik=0,Np-1 
       kk=pgrid(ik)    
       dk=pgrid(ik+1)-kk 
       dks=dk/NpoNq       

       fun_nq=0.d0
       DO ipin=1,NpoNq+1
        qcur=kk+(ipin-1)*dks 
        fun_nq(ipin)=qcur*qcur*N_q(qcur)
       ENDDO

       upli=kk; dnli=kk+dk
       CALL SIMPINT(upli,dnli,NpoNq,fun_nq,sitec)

       bolin=bolin+sitec
      enddo

      N_R0=bolin*8.d0*pi/q_vol                              ! 2 from spin and 4*pi from angles      
      
      END FUNCTION N_R0  

!
!     ======================================================
!     Interpolations for Sigma_HF(momentum)
!     ======================================================
!
      DOUBLE PRECISION FUNCTION SIGMAF(amom)
      IMPLICIT NONE  
      DOUBLE PRECISION, INTENT (IN) :: amom
      INTEGER :: ip_ind, ipc, ipc1, ipc2, ipc3
      DOUBLE PRECISION :: pc1, pc2, pc3
      DOUBLE PRECISION :: Ap1, Ap2, Ap3       
      DOUBLE PRECISION :: Fp1, Fp2, Fp3 
      
      SIGMAF=0.d0

!     Select 3 ponts from pgrid() for quadratic interpolation

      ip_ind=amom/dimp; ipc=ptoi(ip_ind)                    ! index of closest to amom momentum in pgrid

      IF(ipc<0 .OR. ipc>Np) THEN
       print*,' ipc<0 .OR. ipc>Np - SIGMAF!!!'
       print*,' ip_ind=',ip_ind,' ipc=',ipc
       print*,' mom=', amom
       STOP
      ENDIF

      if(ipc==0) then
       ipc1=0; ipc2=1; ipc3=2
      else if(ipc==Np) then
       ipc1=Np-2; ipc2=Np-1; ipc3=Np
      else
       ipc1=ipc-1; ipc2=ipc; ipc3=ipc+1
      endif
      pc1=pgrid(ipc1); pc2=pgrid(ipc2); pc3=pgrid(ipc3)

      Fp1=SigHF_IN(ipc1); Fp2=SigHF_IN(ipc2); Fp3=SigHF_IN(ipc3)

      Ap1=(amom-pc2)*(amom-pc3)/((pc1-pc2)*(pc1-pc3)) 
      Ap2=(amom-pc1)*(amom-pc3)/((pc2-pc1)*(pc2-pc3)) 
      Ap3=(amom-pc1)*(amom-pc2)/((pc3-pc1)*(pc3-pc2)) 
    
      SIGMAF=Fp1*Ap1+Fp2*Ap2+Fp3*Ap3 
      
      END FUNCTION SIGMAF  

!
!     ======================================================
!     Simpson integral
!     ======================================================
!
      SUBROUTINE SIMPINT(ac,bc,nc,fc,si) 
      IMPLICIT NONE  
      INTEGER, INTENT (IN) :: nc
      DOUBLE PRECISION, INTENT (IN) :: ac, bc 
      DOUBLE PRECISION, INTENT (IN) :: fc(32769) 
      DOUBLE PRECISION :: hc, s1, s2, si 
      INTEGER :: ico, ns

      si=0.d0 
      
      hc=dabs((bc-ac)/nc)
      ns=nc/2

      s1=0.d0
      do ico=1,ns
       s1=s1+4.d0*fc(2*ico)
      enddo

      s2=0.d0
      do ico=2,ns
       s2=s2+2.d0*fc(2*ico-1)
      enddo
  
      si=hc*(fc(1)+s1+s2+fc(nc+1))/3.d0 

      END SUBROUTINE SIMPINT

!
!     =======================================================
!     Analytical integral without n_q (a part of Sigma_HF 
!     integral that is taken analytically - it gives the answer 
!     at T=0 and main part of answer at small T). 
!     =======================================================
!
      DOUBLE PRECISION FUNCTION MAIN_INT(kmom)
      IMPLICIT NONE  
      DOUBLE PRECISION, INTENT (IN) :: kmom 
      DOUBLE PRECISION              :: argu, argd 
      DOUBLE PRECISION              :: bolv1, bolv2, bolv3 

      MAIN_INT=0.d0
      
      bolv1=(k_F*k_F+kappa*kappa-kmom*kmom)/(4.d0*kmom)
      
      argu=(kmom+k_F)*(kmom+k_F)+kappa*kappa
      argd=(kmom-k_F)*(kmom-k_F)+kappa*kappa
      bolv2=DLOG(argu/argd)
      
      argu=(k_F+kmom)/kappa
      argd=(k_F-kmom)/kappa   
      bolv3=DATAN(argu)+DATAN(argd)
      
      MAIN_INT=-e2*(k_F+bolv1*bolv2-kappa*bolv3)/pi
      
      END FUNCTION MAIN_INT

!
!     ======================================================
!     Prepare p- and q-grids for G and Pi. In our calculs we 
!     only need pgrid() and ptoi() (but qgrid is auxiliary 
!     here).
!     ======================================================
!
      SUBROUTINE SETGRIDS
      IMPLICIT NONE
      INTEGER :: i, j, k       
      INTEGER :: Iqa1, Iqa2, Ntau2, Nppp  
      DOUBLE PRECISION :: delta, x
      INTEGER :: imy

!     P-grid
    
      allocate(pgrid(0:Np))
      allocate(qgrid(0:Nq))
      
      Ntau2=Np1*Kp2
      pgrid=-100.d0
      kc=k_F+dimp*Np1*(2**Kp3-1.d0)                         ! kc adjusted to p step 
      print*,' '
      print*, ' kc  =', kc, Np

      delta=dimp
      k=0
      qgrid(k)=0.d0; x=qgrid(k)      
      do j=1,Kp2
       Nppp=Np1
       do i=1,Nppp
        qgrid(k+i)=x+delta*i
       enddo
       k=k+Nppp
       x=qgrid(k); delta=delta*2.d0    
      enddo

      i=Ntau2+1; qgrid(i)=k_F; ik_F=i      

      IF(qgrid(Ntau2) < 0.d0) stop 'check text above; Pgrid'
      IF(DABS(qgrid(Ntau2)-k_F+Np1*dimp)>1.d-14) stop 'grid P'      

      do j=0,i
       pgrid(i-j)=k_F-qgrid(j)
      enddo
          
      delta=dimp
      k=i
      x=pgrid(i)   
      do j=1,Kp3
       Nppp=Np1 
       do i=1,Nppp
        pgrid(k+i)=x+delta*i
       enddo 
       k=k+Nppp
       x=pgrid(k); delta=delta*2.d0    
      enddo                                      

      IF(DABS(pgrid(Np)-kc)>1.d-14) stop 'pgrid Pc'       
   
      i=powerp+(Kp3-Kp2); i=2**powerp+2**i-Np1 
      IF(DABS(dimp*i-kc)>1.d-12) stop 'check ptoi limits' 

      pgrid(0)=pgrid(1)/1000.d0 

      k=0
      do j=0,i-1
       x=j*dimp+1.d-12
       IF(pgrid(k+1) < x) k=k+1
       ptoi(j)=k
      enddo
      ptoi(i)=Np-1  

!     ======================================================
!     Going back from momentum p to grid index: (1) k_p=p/dimp; 
!     then j=ptoi(k_p) gives the index of closest to p momentum 
!     p* in pgrid: p*=pgrid(j). 
!     ======================================================
     
!     Q-grid

      Ntau2=Nq1*Kq2
      qgrid=-100.d0
      kqc=2.d0*k_F+dimp*Nq1*(2**Kq3-1.d0)                   ! kqc adjusted to q step 
      print*, ' kqc =', kqc, Nq
           
      delta=dimp
      k=0
      qgrid(k)=0.d0; x=qgrid(k)
      do j=1,Kq2
       do i=1,Nq1
        qgrid(k+i)=x+delta*i 
       enddo      
       k=k+Nq1
       x=qgrid(k); delta=delta*2.d0    
      enddo                                                       

      IF(qgrid(Ntau2) < 0.d0) stop 'check text above; qgrid'
      IF(DABS(qgrid(Ntau2)-k_F+Nq1*dimp)>1.d-14) stop 'grid Q'          

      i=2*Ntau2+1
      do j=0,Ntau2
       qgrid(i-j)=2.d0*k_F-qgrid(j)
      enddo
      
!     ======================================================      
      qgrid(0)=min(2.d-06, dimp/2.d0)                       ! first q point 
!     ======================================================      
      
      delta=dimp
      k=i
      x=qgrid(i)
      do j=1,Kq3
       do i=1,Nq1
        qgrid(k+i)=x+delta*i
       enddo      
       k=k+Nq1
       x=qgrid(k); delta=delta*2.d0    
      enddo                                      
      IF(DABS(qgrid(Nq)-kqc)>1.d-14) stop 'tgrid Qc' 

      i=powerp+(Kq3-Kq2); i=2**powerp+2**powerp+2**i-Nq1
      k=0
      do j=0,i-1
       x=j*dimp+1.d-12
       IF(qgrid(k+1) < x) k=k+1
       qtoi(j)=k
      enddo 
      qtoi(i)=Nq-1  

!     ======================================================
!     Going back from momentum q to grid index: (1) q_p=q/dimp; 
!     then jj=qtoi(q_p) gives the index of closest to q momentum  
!     q* in qgrid: q*=qgrid(jj).         
!     ======================================================  

      END SUBROUTINE SETGRIDS

      
      
      END PROGRAM DISP_HF


