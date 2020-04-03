C FILE: obs_py.f

        program test_obs
        !subroutine test_obs(OBSR)
        !comment(and uncomment) line above(bellow) when compiling for python(f2py)
	  !subroutine test_obs(isn,wt,q2t,thi,phii,eps,lmx,t_obs,nob,
c     >  parf,np,ragami,iagami,rTGa,iTGa,cmet,p_count )
	  
	  implicit none
	  
	  real*8 OBSR(68),parf(500),t_obs, t_amp
	  real*8 wt,q2t,thi,phii,eps,cme(260),cmet(260)
	  real*8 ragami(3,5,260,4,13,2),iagami(3,5,260,4,13,2)
	  real*8 rTGa(13,13,5,4,260,3,2),iTGa(13,13,5,4,260,3,2)
	  real*8 rmat,imat
	  complex*8 a_gami(3,5,260,4,13,2),TGa(13,13,5,4,260,3,2)
	  complex*8 m_gam(4,5,13,3,2),mat
	  integer ic1,ic2,kj,iv,iem,ien,j,nob,iq2,i
	  integer par_cnt(13,5,4,3),ipar,np,lmx,isn,p_count(13,5,4,3)
	  character par_file*12
	  common /mtrx/ a_gami,TGa,cme,par_cnt
	  
Cf2py intent(in) isn,wt,q2t,thi,phii,eps,lmx,nob
Cf2py intent(in) parf
Cf2py intent(in) ragami,iagami,rTGa,iTGa,cmet,p_count 
Cf2py intent(out) t_obs
	!---------------------------------------------------------------
	  parf=0.	
	  call creater(1,ragami,iagami,rTGa,iTGa,cmet)
	  call par_count(par_cnt)
	  par_file='par_strt.txt'
	  call get_par(par_file,parf,np)
	!---------------------------------------------------------------
	
	  do ic1=1,13
	  do kj=1,5
	  do iv=1,4
	  do iem=1,2
	  do j=1,3
	  do ien=1,260
	  rmat=ragami(j,kj,ien,iv,ic1,iem)
	  imat=iagami(j,kj,ien,iv,ic1,iem)
	  mat=complex(rmat,imat)
	  a_gami(j,kj,ien,iv,ic1,iem)=mat
	  do ic2=1,13
	  rmat=rTGa(ic1,ic2,kj,iv,ien,j,iem)
	  imat=iTGa(ic1,ic2,kj,iv,ien,j,iem)
	  mat=complex(rmat,imat)
	  TGa(ic1,ic2,kj,iv,ien,j,iem)=mat
	  enddo
	  enddo
	  enddo
	  enddo
	  enddo
	  enddo
	  enddo
	  
	  cme=cmet
	  !par_cnt=p_count
	  
	  m_gam=0.
	  OBSR=0.
	  
	!---------------------------------------------------------------
	  isn=1
	  wt=1233.
	  q2t=2
	  thi=90.
	  phii=90.
	  eps=0.95
	  lmx=4
	  nob=42
	  t_obs=0.
	 
	!---------------------------------------------------------------
       do i=1,100
       q2t=float(i-1)*0.04
	  call get_obs(wt,q2t,thi,phii,eps,isn,parf(1:np),OBSR,lmx)
	  t_obs= OBSR(nob)
	  print*,t_obs, q2t
	 enddo 
	  
	!---------------------------------------------------------------
	  
	  !print*,t_obs
	  
	  !call get_M(q2t,wt,parf(1:np),m_gam,lmx)  
	  !t_amp=aimag(m_gam(3,2,1,2,1))
	  
	  
	  
	  !print*,m_gam(3,2,1,2,1)
	  
	  
	  end
      

      
      !===============================================================
       
      subroutine get_par(par_file,parf,np)
      implicit none
      integer np
      real*8 parf(500)
      character par_file*12
Cf2py intent(in) par_file      
Cf2py intent(out) parf,np      
      
	  open(unit=44,file=par_file)
	  
	  np=188
		   
	  parf=0.5
		   
	  read(44,*)parf(1:np)
	  
	  np=500
		  
	  close(unit=44)
	  end
      
       
      !=============================================================== 
      subroutine get_obs(wt,q2t,thi,phii,eps,isn,parf,OBSR,lmx)
          
      !---------------------------------------------------------------   
	!   Calculate observables for a given set of multipole.
	!   Helicity amplitudes and observable are constructed
	!   from expresions used in PRC 96, 025210 (2017).
	!--------------------------------------------------------------- 
       
	implicit none
	integer imp,kj,iv,iem,isn,k_j,i_v,il,lmx,np
	real*8 wt,q2t,th,phi,rpi,thi,phii,eps
	complex*8 m_gam(4,5,13,3,2),m_gam_pi0p(2,5,3)
	complex*8 m_gam_pipn(2,5,3)
	complex*8 H1,H2,H3,H4,H5,H6
	complex*8 Ep(0:5),Em(0:5),Mp(0:5),Mm(0:5),Lp(0:5),Lm(0:5)
	real*8 OBSR(68),parf(*)
      logical there
      
      
Cf2py intent(in) isn,wt,q2t,thi,phii,eps,lmx,nob
Cf2py intent(in) parf
Cf2py intent(out) OBSR
      
      !print*,parf
       
	rpi=dacos(-1.d0)
	th=thi/180.*rpi
	phi=phii/180.*rpi
	
	m_gam=0.
	m_gam_pi0p=0.
	m_gam_pipn=0.
	
	Ep=0.
	Em=0.
	Mp=0.
	Mm=0.
	Lp=0.
	Lm=0. 
	
	call get_M(q2t,Wt,parf,m_gam,lmx)   
	!print*,q2t,Wt, m_gam(3,2,1,2,1)	
	do kj=1,5
	do iv=1,2
	il=iv+1+kj-3
	if(il.le.lmx)then
	do iem=1,3 !multipole index; ie_m=1 for electric, =2 for magnetic  and =3 for longitudinal
	m_gam_pi0p(iv,kj,iem)=m_gam(iv,kj,1,iem,1)+
     >2./3.*m_gam(iv+2,kj,1,iem,1)! pi_0 proton final state
	m_gam_pipn(iv,kj,iem)=sqrt(2.)*( m_gam(iv,kj,1,iem,1)-
     >1./3.*m_gam(iv+2,kj,1,iem,1))! pi_+ neutron final state
	enddo  !i_em
	endif
	enddo  !iv
	enddo  !kj
	
	! Generate  set of multiploles E,M, and L expanded to l=5 for a given center of mass energy  W
	
	do k_j=1,5
	do i_v=1,2
	il=i_v+1+k_j-3    !pion orbital angular momentum
	if(isn.eq.1)then !pi_0 proton multipoles:
	if(i_v.eq.1)then
	Ep(il)=m_gam_pi0p(i_v,k_j,1)
	Mp(il)=m_gam_pi0p(i_v,k_j,2)
	Lp(il)=m_gam_pi0p(i_v,k_j,3)
	elseif(i_v.eq.2)then
	Em(il)=m_gam_pi0p(i_v,k_j,1)
	Mm(il)=m_gam_pi0p(i_v,k_j,2)
	Lm(il)=m_gam_pi0p(i_v,k_j,3)
	endif
	else  !pi_+ neutron multipoles:
	if(i_v.eq.1)then
	Ep(il)=m_gam_pipn(i_v,k_j,1)
	Mp(il)=m_gam_pipn(i_v,k_j,2)
	Lp(il)=m_gam_pipn(i_v,k_j,3)
	elseif(i_v.eq.2)then
	Em(il)=m_gam_pipn(i_v,k_j,1)
	Mm(il)=m_gam_pipn(i_v,k_j,2)
	Lm(il)=m_gam_pipn(i_v,k_j,3)
	endif
	endif
	enddo
	enddo
       
      rpi=dacos(-1.d0)
      
      OBSR=0.
       
	call h_to_obs(Wt,q2t,Ep,Em,Mp,Mm,Lp,Lm,th,phi,eps,
     > H1,H2,H3,H4,H5,H6,OBSR,lmx)
	
66    end 

	!==============================================================
      
      Subroutine h_to_obs(wt,q2,Ep,Em,Mp,Mm,Lp,Lm,th,phi,eps,
     >H1,H2,H3,H4,H5,H6,OBSR,lmx)
       
      implicit none
      integer ibeg,lmx
      real*8 wt, q2,WCM
	real*8 M_P, M_PI0, kcm,qcm, RHO, CONV, CONV2
	real*8 th,SQ2,rpi,phi,eps,epsi,epfc1,epfc2,epfc3
	complex*8 Ep(0:5),Em(0:5),Mp(0:5),Mm(0:5),Lp(0:5),Lm(0:5)
	complex*8 H1,H2,H3,H4,H5,H6
	complex*8 CH1, CH2, CH3, CH4, CH5, CH6
	real*8 RT00,RT0Y,RTYP0,RTXPX,RTXPZ,RTZPX,RTZPZ,
     > RL00,RL0Y,RLXPX,RLZPX,RLT00,RLT0X,RLT0Y,RLT0Z,RLTXP0,RLTZP0, 
     > RLTXPX,RLTZPX,RTT00,RTT0X,RTT0Z,RTTXP0,RTTZP0,RLTP00,RLTP0X, 
     > RLTP0Y,RLTP0Z,RLTPXP0,RLTPZP0,RLTPXPX,RLTPZPX,RTTP0X,RTTP0Z,
     > RTTPXP0,RTTPZP0,dsT,dsL,dsTT,dsTL,dsTLP,DSG,T,P,TXP,LXP,TZP,
     > LZP,S,H,G,OX,OZ,F,E,CXP,CZP,A0Y,A0YP,dTLs,PXh,PY,PZh,RTLs,
     > RTLsn,AE,AT,AET,OBSR(68)
     
     
     
	OBSR=0.d0
	 
	M_P=0.93827		!Proton mass in GeV		
	M_PI0=0.13498	!Pion mass in GeV
	WCM=wt/1000.		!Center of mass energy in GeV
	rpi=dacos(-1.d0)
       
	SQ2 =dsqrt(2.d0)
       
      !---------------------------------------------------------------
      ! Select method for generating helicity amplitudes
      ! 
	! ibeg = 1  use Helicity amplitudes generated in Mathematica (see mltp_tohel_obs4.nb).
      ! ibeg = 2  Amplitudes calculated here from Legendre polynomials
      !
      !---------------------------------------------------------------
	ibeg=1 
       
	if(ibeg.eq.1)then
	call gen_h5(Ep,Em,Mp,Mm,Lp,Lm,th,H1,H2,H3,H4,H5,H6)
	else
	call gen_h(Ep,Em,Mp,Mm,Lp,Lm,th,H1,H2,H3,H4,H5,H6,lmx)
	endif
	 
	 	!    Here, complex conjugate for H amplitudes
      
      CH1=CONJG(H1)
      CH2=CONJG(H2)
      CH3=CONJG(H3)
      CH4=CONJG(H4)
      CH5=CONJG(H5)
      CH6=CONJG(H6)
      
      !   Observables. Using  notation introduced by Ron in test.f
      
      RT00 = 0.5d0*(H1*CH1+H2*CH2+H3*CH3+H4*CH4)
      RT0Y = -Aimag(H2*CH1+H4*CH3)
      RTYP0 = Aimag(H3*CH1+H4*CH2)
      RTXPX = Real(H4*CH1+H3*CH2)
      RTXPZ = Real(H3*CH1-H4*CH2)
      RTZPX = Real(H2*CH1-H4*CH3)
      RTZPZ = 0.5d0*(H1*CH1-H2*CH2-H3*CH3+H4*CH4)
      
      RL00 = Real(H5*CH5+H6*CH6)
      RL0Y = -2.d0*Aimag(H6*CH5)
      RLXPX = Real(-H5*CH5+H6*CH6)
      RLZPX = 2.d0*Real(H6*CH5)
      
      RLT00 = Real( (H1-H4)*CH5+(H2+H3)*CH6 )/SQ2
      RLT0X = Aimag( (H3-H2)*CH5-(H1+H4)*CH6 )/SQ2
      RLT0Y = -Aimag( (H2+H3)*CH5-(H1-H4)*CH6 )/SQ2
      RLT0Z = -Aimag( (H1+H4)*CH5-(H2-H3)*CH6 )/SQ2
      RLTXP0 = Aimag( (H2-H3)*CH5-(H1+H4)*CH6 )/SQ2
      RLTZP0 = -Aimag( (H1+H4)*CH5+(H2-H3)*CH6 )/SQ2
      RLTXPX = -Real( (H1-H4)*CH5-(H2+H3)*CH6 )/SQ2
      RLTZPX = Real( (H2+H3)*CH5+(H1-H4)*CH6 )/SQ2
      
      RTT00 = Real(H3*CH2-H4*CH1)
      RTT0X = Aimag(H3*CH1-H4*CH2)
      RTT0Z = -Aimag(H4*CH1+H3*CH2)
      RTTXP0 = Aimag(H2*CH1-H4*CH3)
      RTTZP0 = Aimag(H3*CH2-H4*CH1)
      
      RLTP00 = -Aimag( (H1-H4)*CH5+(H2+H3)*CH6 )/SQ2
      RLTP0X = Real( (H2-H3)*CH5+(H1+H4)*CH6 )/SQ2
      RLTP0Y = -Real( (H2+H3)*CH5+(H4-H1)*CH6 )/SQ2
      RLTP0Z = Real( (H1+H4)*CH5+(H3-H2)*CH6 )/SQ2
      RLTPXP0 = Real( (H3-H2)*CH5+(H1+H4)*CH6 )/SQ2
      RLTPZP0 = Real( (H1+H4)*CH5+(H2-H3)*CH6 )/SQ2
      RLTPXPX = Aimag( (H1-H4)*CH5-(H2+H3)*CH6 )/SQ2
      RLTPZPX = -Aimag( (H2+H3)*CH5+(H1-H4)*CH6 )/SQ2
      
      RTTP0X =Real(H2*CH1+H4*CH3)
      RTTP0Z =0.5d0*(H1*CH1-H2*CH2+H3*CH3-H4*CH4)
      RTTPXP0 =Real(H3*CH1+H4*CH2)
      RTTPZP0 =0.5d0*(H1*CH1-H2*CH2-H3*CH3-H4*CH4)
      
      !     Taken from test.f code from Ron
      
      KCM = (WCM**2 - M_P**2)/(2.*M_P)*(M_P/WCM)
      RHO = (WCM**2 - M_P**2 + M_PI0**2)/(2.*WCM)
      RHO = sqrt(RHO**2 - M_PI0**2)/KCM
            
      CONV2 = sqrt(q2)*2.*WCM!/(WCM**2-q2-M_P**2) !sqrt(EPSL/EPS)
      epfc1 = sqrt(2.*eps*(1.+eps))
      epfc2 = sqrt(2.*eps*(1.-eps))
      epfc3 = sqrt(1.-eps**2)
      
        
      dsT = RHO*RT00
      dsL = RHO*RL00*CONV2**2
      dsTT = RHO*RTT00
      dsTL = RHO*RLT00*CONV2
      dsTLP = RHO*RLTP00*CONV2      
      
      DSG=dsT+eps*dsL+epfc1*dsTL*cos(phi)
     >  +eps*dsTT*cos(2.*phi) 
      
      !Photoproduction Observables
      
      T=RT0Y/RT00
      P=RTYP0/RT00
      TXP=RTXPX/RT00
      LXP=-RTXPZ/RT00
      TZP=RTZPX/RT00
      LZP= RTZPZ/RT00
      S=-RTT00/RT00
      H=RTT0X/RT00
      G=-RTT0Z/RT00
      OX=RTTXP0/RT00
      OZ=RTTZP0/RT00
      F=RTTP0X/RT00
      E=-RTTP0Z/RT00
      CXP=-RTTPXP0/RT00
      CZP=-RTTPZP0/RT00
      
      A0Y=-epfc1*RLT00/(RT00+eps*(RL00+RTT00))
      A0YP=(dsTLP/DSG)*epfc2
      dTLs=dsT+eps*dsL
      PXh=-epfc2*RLTP0X/(RT00+eps*RL00)
      PY=-epfc2*RLT0Y/(RT00+eps*RL00)
      PZh=-epfc3*RTTP0Z/(RT00+eps*RL00)
      RTLs=RT00+eps*RL00
      RTLsn=RTYP0-eps*RL0Y
      
      AE=epfc2*CONV2*(RLTP00/DSG)*sin(phi)
      
	AT=(sin(th)*(cos(phi)*(epfc1*RLT0X*sin(phi)+
     >  eps*RTT0X*sin(2.*phi))+sin(phi)*(RLT0Y+eps*(CONV2*RL0Y
     >  + epfc1*RLT0Y*cos(phi)-RTYP0*cos(2.*phi))))
     >  +cos(th)*eps*(epfc1*RLT0Z*sin(phi)+RTT0Z*sin(2.*phi)))/DSG
     
	AET=0.!-eps*(sin(th)*(epfc2*RLTP0X*cos(phi)**2
      !>   +epfc3/eps*RTTP0X*cos(phi)+epfc2*RLTP0Y*sin(th)**2)
      !>   +cos(th)*(epfc2*RLTP0Z*cos(phi)+epfc3/eps* RTTP0Z))/DSG
	


       OBSR(1:68)=(/RT00,RT0Y,RTYP0,RTXPX,RTXPZ,RTZPX,RTZPZ,RL00,  
     > RL0Y,RLXPX,RLZPX,RLT00,RLT0X,RLT0Y,RLT0Z,RLTXP0,RLTZP0,RLTXPX, 
     > RLTZPX,RTT00,RTT0X,RTT0Z,RTTXP0,RTTZP0,RLTP00,RLTP0X,RLTP0Y,  
     > RLTP0Z,RLTPXP0,RLTPZP0,RLTPXPX,RLTPZPX,RTTP0X,RTTP0Z,RTTPXP0,
     > RTTPZP0,dsT,dsL,dsTT,dsTL,dsTLP,DSG,T,P,TXP,LXP,TZP,LZP,
     > S,H,G,OX,OZ,F,E,CXP,CZP,A0Y,A0YP,dTLs,PXh,PY,PZh,RTLs,RTLsn,
     > AE,AT,AET/)*(0.1d0)**2
      
      End Subroutine h_to_obs
         
	!================================================================
      
      subroutine gen_h5(Ep,Em,Mp,Mm,Lp,Lm,th,H1,H2,H3,H4,H5,H6)
	!----------------------------------------------------------------
	! Helicity amplitudes written explicitly in terms of multipoles
	! after procedure outline in PRC 96, 025210 (2017)
      ! valid up to l = 5
      !---------------------------------------------------------------- 
      
      implicit none
      
      complex*8 Ep(0:5),Em(0:5),Mp(0:5),Mm(0:5),Lp(0:5),Lm(0:5),
     >  H1,H2,H3,H4,H5,H6
      real*8 th
     
       H1=0.+(1.0606601717798214*Em(2)-1.060660171779821*Em(3)+  
     > 1.988737822087165*Em(4)-1.9887378220871632*Em(5)-         
     > 1.0606601717798214*Ep(1)+1.0606601717798212*Ep(2)-  
     > 1.988737822087165*Ep(3)+1.9887378220871632*Ep(4)-  
     > 2.9002426572104483*Ep(5)+1.0606601717798214*Mm(2)-  
     > 1.060660171779821*Mm(3)+1.988737822087165*Mm(4)-  
     > 1.9887378220871632*Mm(5)+1.0606601717798214*Mp(1)-  
     > 1.060660171779821*Mp(2)+1.988737822087165*Mp(3)-  
     > 1.9887378220871632*Mp(4)+  
     > 2.9002426572104483*Mp(5))*dSin(0.5*th)+  
     > (1.0606601717798214*Em(2)+1.590990257669732*Em(3)-  
     > 0.6629126073623883*Em(4)+2.6516504294495533*Em(5)-  
     > 1.0606601717798214*Ep(1)-1.5909902576697321*Ep(2)+  
     > 0.6629126073623883*Ep(3)-2.651650429449555*Ep(4)+  
     > 1.7401455943262665*Ep(5)+1.0606601717798214*Mm(2)+  
     > 1.590990257669732*Mm(3)-0.6629126073623883*Mm(4)+  
     > 2.6516504294495533*Mm(5)+1.0606601717798214*Mp(1)+  
     > 1.590990257669732*Mp(2)-0.6629126073623883*Mp(3)+  
     > 2.6516504294495533*Mp(4)-  
     > 1.7401455943262647*Mp(5))*dSin(1.5*th)+  
     > (2.6516504294495533*Em(3)+1.988737822087165*Em(4)-  
     > 2.6516504294495533*Ep(2)-1.988737822087165*Ep(3)+  
     > 8.881784197001252e-16*Ep(4)-3.190266922931496*Ep(5)+  
     > 2.6516504294495533*Mm(3)+1.988737822087165*Mm(4)+  
     > 2.6516504294495533*Mp(2)+1.988737822087165*Mp(3)+  
     > 3.190266922931496*Mp(5))*dSin(2.5*th)+  
     > (4.640388251536718*Em(4)+2.32019412576836*Em(5)-  
     > 4.640388251536718*Ep(3)-2.32019412576836*Ep(4)-  
     > 0.8700727971631341*Ep(5)+4.640388251536718*Mm(4)+  
     > 2.32019412576836*Mm(5)+4.640388251536718*Mp(3)+  
     > 2.32019412576836*Mp(4)+0.8700727971631359*Mp(5))*dSin(3.5*th)  
     > +(6.960582377305078*Em(5)-6.960582377305078*Ep(4)-  
     > 2.610218391489404*Ep(5)+6.960582377305078*Mm(5)+  
     > 6.960582377305078*Mp(4)+2.610218391489404*Mp(5))*dSin(4.5*th)  
     > +(-9.570800768794482*Ep(5)+  
     > 9.570800768794482*Mp(5))*dSin(5.5*th)
     
      H2=0.+dCos(5.5*th)*(-13.399121076312273*Ep(5)-  
     > 9.570800768794484*Mp(5))+  
     > dCos(3.5*th)*(-4.640388251536718*Em(4)-  
     > 0.7733980419227855*Em(5)-7.733980419227864*Ep(3)-  
     > 1.1600970628841782*Ep(4)-6.090509580141944*Ep(5)+  
     > 7.733980419227864*Mm(4)+1.16009706288418*Mm(5)-  
     > 4.640388251536718*Mp(3)-0.7733980419227873*Mp(4)-  
     > 4.350363985815697*Mp(5))+  
     > dCos(1.5*th)*(-1.0606601717798214*Em(2)-  
     > 0.5303300858899109*Em(3)-1.988737822087165*Em(4)-  
     > 1.3258252147247767*Em(5)-3.181980515339464*Ep(1)-  
     > 1.060660171779821*Ep(2)-3.3145630368119416*Ep(3)-  
     > 1.988737822087165*Ep(4)-4.060339720094632*Ep(5)+  
     > 3.1819805153394642*Mm(2)+1.060660171779821*Mm(3)+  
     > 3.3145630368119416*Mm(4)+1.9887378220871668*Mm(5)-  
     > 1.0606601717798214*Mp(1)-0.5303300858899105*Mp(2)-  
     > 1.9887378220871632*Mp(3)-1.3258252147247767*Mp(4)-  
     > 2.9002426572104554*Mp(5))+  
     > dCos(0.5*th)*(-0.35355339059327373*Em(2)-  
     > 1.0606601717798214*Em(3)-1.1932426932522997*Em(4)-  
     > 1.9887378220871668*Em(5)-1.4142135623730951*Ep(0)-  
     > 1.0606601717798214*Ep(1)-2.121320343559642*Ep(2)-  
     > 1.988737822087165*Ep(3)-2.98310673313075*Ep(4)-  
     > 2.9002426572104483*Ep(5)+1.4142135623730951*Mm(1)+  
     > 1.0606601717798214*Mm(2)+2.121320343559642*Mm(3)+  
     > 1.9887378220871632*Mm(4)+2.983106733130743*Mm(5)-  
     > 0.35355339059327395*Mp(1)-1.0606601717798219*Mp(2)-  
     > 1.1932426932522953*Mp(3)-1.9887378220871614*Mp(4)-  
     > 2.071601898007451*Mp(5))+  
     > dCos(2.5*th)*(-2.6516504294495533*Em(3)-  
     > 0.6629126073623883*Em(4)-3.0935921676911464*Em(5)-  
     > 5.303300858899107*Ep(2)-1.1048543456039805*Ep(3)-  
     > 4.640388251536717*Ep(4)-2.030169860047309*Ep(5)+  
     > 5.303300858899107*Mm(3)+1.1048543456039805*Mm(4)+  
     > 4.6403882515367165*Mm(5)-2.6516504294495533*Mp(2)-  
     > 0.6629126073623883*Mp(3)-3.0935921676911384*Mp(4)-  
     > 1.4501213286052383*Mp(5))+  
     > dCos(4.5*th)*(-6.960582377305078*Em(5)-  
     > 10.440873565957618*Ep(4)-1.21810191602839*Ep(5)+  
     > 10.440873565957615*Mm(5)-6.9605823773050775*Mp(4)-  
     > 0.8700727971631359*Mp(5))
      
      H3=0.+dCos(2.5*th)*(-2.6516504294495533*Em(3)+  
     > 1.988737822087165*Em(4)+8.881784197001252e-16*Em(5)-  
     > 2.6516504294495533*Ep(2)+1.988737822087165*Ep(3)+  
     > 8.881784197001252e-16*Ep(4)+3.190266922931496*Ep(5)-  
     > 2.6516504294495533*Mm(3)+1.988737822087165*Mm(4)+  
     > 8.881784197001252e-16*Mm(5)+2.6516504294495533*Mp(2)-  
     > 1.988737822087165*Mp(3)-3.190266922931496*Mp(5))+  
     > dCos(0.5*th)*(1.0606601717798214*Em(2)+  
     > 1.0606601717798212*Em(3)+1.988737822087165*Em(4)+  
     > 1.9887378220871632*Em(5)+1.0606601717798214*Ep(1)+  
     > 1.0606601717798212*Ep(2)+1.988737822087165*Ep(3)+  
     > 1.9887378220871632*Ep(4)+2.9002426572104483*Ep(5)+  
     > 1.0606601717798214*Mm(2)+1.0606601717798212*Mm(3)+  
     > 1.988737822087165*Mm(4)+1.9887378220871632*Mm(5)-  
     > 1.0606601717798214*Mp(1)-1.060660171779821*Mp(2)-  
     > 1.988737822087165*Mp(3)-1.9887378220871632*Mp(4)-  
     > 2.9002426572104483*Mp(5))+  
     > dCos(4.5*th)*(-6.960582377305078*Em(5)-  
     > 6.960582377305078*Ep(4)+2.610218391489404*Ep(5)-  
     > 6.960582377305078*Mm(5)+6.960582377305078*Mp(4)-  
     > 2.610218391489404*Mp(5))+  
     > dCos(1.5*th)*(-1.0606601717798214*Em(2)+  
     > 1.590990257669732*Em(3)+0.6629126073623883*Em(4)+  
     > 2.6516504294495533*Em(5)-1.0606601717798214*Ep(1)+  
     > 1.590990257669732*Ep(2)+0.6629126073623883*Ep(3)+  
     > 2.6516504294495533*Ep(4)+1.7401455943262665*Ep(5)-  
     > 1.0606601717798214*Mm(2)+1.590990257669732*Mm(3)+  
     > 0.6629126073623883*Mm(4)+2.6516504294495533*Mm(5)+  
     > 1.0606601717798214*Mp(1)-1.5909902576697321*Mp(2)-  
     > 0.6629126073623883*Mp(3)-2.651650429449555*Mp(4)-  
     > 1.7401455943262647*Mp(5))+  
     > dCos(3.5*th)*(-4.640388251536718*Em(4)+2.32019412576836*Em(5)  
     > -4.640388251536718*Ep(3)+2.32019412576836*Ep(4)-  
     > 0.8700727971631341*Ep(5)-4.640388251536718*Mm(4)+  
     > 2.32019412576836*Mm(5)+4.640388251536718*Mp(3)-  
     > 2.32019412576836*Mp(4)+0.8700727971631359*Mp(5))+  
     > dCos(5.5*th)*(-9.570800768794482*Ep(5)+  
     > 9.570800768794482*Mp(5))
      
      H4=0.+(0.35355339059327373*Em(2)-1.0606601717798214*Em(3)+  
     > 1.1932426932522997*Em(4)-1.988737822087165*Em(5)+  
     > 1.4142135623730951*Ep(0)-1.0606601717798214*Ep(1)+  
     > 2.1213203435596424*Ep(2)-1.988737822087165*Ep(3)+  
     > 2.9831067331307484*Ep(4)-2.9002426572104554*Ep(5)+  
     > 1.4142135623730951*Mm(1)-1.0606601717798212*Mm(2)+  
     > 2.12132034355964*Mm(3)-1.9887378220871632*Mm(4)+  
     > 2.9831067331307395*Mm(5)-0.35355339059327395*Mp(1)+  
     > 1.0606601717798219*Mp(2)-1.1932426932522953*Mp(3)+  
     > 1.9887378220871597*Mp(4)-  
     > 2.071601898007465*Mp(5))*dSin(0.5*th)+  
     > (-1.0606601717798214*Em(2)+0.5303300858899109*Em(3)-  
     > 1.988737822087165*Em(4)+1.3258252147247749*Em(5)+  
     > 3.1819805153394642*Ep(1)-1.060660171779821*Ep(2)+  
     > 3.3145630368119416*Ep(3)-1.9887378220871632*Ep(4)+  
     > 4.060339720094632*Ep(5)+3.1819805153394642*Mm(2)-  
     > 1.060660171779821*Mm(3)+3.3145630368119416*Mm(4)-  
     > 1.9887378220871677*Mm(5)+1.0606601717798214*Mp(1)-  
     > 0.5303300858899105*Mp(2)+1.9887378220871632*Mp(3)-  
     > 1.3258252147247749*Mp(4)+  
     > 2.9002426572104767*Mp(5))*dSin(1.5*th)+  
     > (-2.6516504294495533*Em(3)+0.6629126073623883*Em(4)-  
     > 3.0935921676911464*Em(5)+5.303300858899107*Ep(2)-  
     > 1.1048543456039805*Ep(3)+4.64038825153672*Ep(4)-  
     > 2.030169860047316*Ep(5)+5.303300858899107*Mm(3)-  
     > 1.1048543456039805*Mm(4)+4.640388251536713*Mm(5)+  
     > 2.6516504294495533*Mp(2)-0.6629126073623883*Mp(3)+  
     > 3.0935921676911384*Mp(4)-  
     > 1.4501213286052366*Mp(5))*dSin(2.5*th)+  
     > (-4.640388251536718*Em(4)+0.7733980419227855*Em(5)+  
     > 7.733980419227864*Ep(3)-1.1600970628841782*Ep(4)+  
     > 6.090509580141946*Ep(5)+7.733980419227864*Mm(4)-  
     > 1.1600970628841782*Mm(5)+4.640388251536718*Mp(3)-  
     > 0.7733980419227873*Mp(4)+4.35036398581569*Mp(5))*dSin(3.5*th)  
     > +(-6.960582377305078*Em(5)+10.440873565957617*Ep(4)-  
     > 1.2181019160283881*Ep(5)+10.440873565957615*Mm(5)+  
     > 6.960582377305077*Mp(4)-  
     > 0.8700727971631341*Mp(5))*dSin(4.5*th)+  
     > (13.399121076312277*Ep(5)+  
     > 9.570800768794482*Mp(5))*dSin(5.5*th)
      
      H5=16.2421875*dCos(5.5*th)*Lp(5)+dCos(4.5*th)*(12.3046875*Lm(5)+  
     > 12.3046875*Lp(4)+1.4765625*Lp(5))+dCos(2.5*th)*(5.625*Lm(3)+  
     > 1.25*Lm(4)+5.46875*Lm(5)+5.625*Lp(2)+1.25*Lp(3)+  
     > 5.46875*Lp(4)+2.4609375*Lp(5))+dCos(0.5*th)*(1.*Lm(1)+  
     > 1.*Lm(2)+2.25*Lm(3)+2.25*Lm(4)+3.515625*Lm(5)+1.*Lp(0)+  
     > 1.*Lp(1)+2.25*Lp(2)+2.25*Lp(3)+3.515625*Lp(4)+  
     > 3.515625*Lp(5))+dCos(1.5*th)*(3.*Lm(2)+1.125*Lm(3)+3.75*Lm(4)  
     > +2.34375*Lm(5)+3.*Lp(1)+1.125*Lp(2)+3.75*Lp(3)+2.34375*Lp(4)  
     > +4.921875*Lp(5))+dCos(3.5*th)*(8.75*Lm(4)+1.3671875*Lm(5)+  
     > 8.75*Lp(3)+1.3671875*Lp(4)+7.3828125*Lp(5))
              
      H6=(1.*Lm(1)-1.*Lm(2)+2.25*Lm(3)-2.25*Lm(4)+3.515625*Lm(5)-  
     > 1.*Lp(0)+1.*Lp(1)-2.25*Lp(2)+2.25*Lp(3)-3.515625*Lp(4)+  
     > 3.515625*Lp(5))*dSin(0.5*th)+(3.*Lm(2)-1.125*Lm(3)+3.75*Lm(4)  
     > -2.34375*Lm(5)-3.*Lp(1)+1.125*Lp(2)-3.75*Lp(3)+2.34375*Lp(4)  
     > -4.921875*Lp(5))*dSin(1.5*th)+(5.625*Lm(3)-1.25*Lm(4)+  
     > 5.46875*Lm(5)-5.625*Lp(2)+1.25*Lp(3)-5.46875*Lp(4)+  
     > 2.4609375*Lp(5))*dSin(2.5*th)+(8.75*Lm(4)-1.3671875*Lm(5)-  
     > 8.75*Lp(3)+1.3671875*Lp(4)-7.3828125*Lp(5))*dSin(3.5*th)+  
     > (12.3046875*Lm(5)-12.3046875*Lp(4)+  
     > 1.4765625*Lp(5))*dSin(4.5*th)-16.2421875*Lp(5)*dSin(5.5*th)
     
      end
      
      !===============================================================
      
      
      subroutine gen_h(Ep,Em,Mp,Mm,Lp,Lm,thr,H1,H2,H3,H4,H5,H6,lmx)
      !--------------------------------------------------------------
      ! Generates helicity amplitudes from multipoles through steps
      ! highlighted in PRC 96, 025210 (2017)
      !--------------------------------------------------------------
      
      implicit none
      
      integer i,lmx,l
      complex*8 Ep(0:lmx),Em(0:lmx),Mp(0:lmx),Mm(0:lmx),
     >  Lp(0:lmx),Lm(0:lmx)
      complex*8 fc(6),H1,H2,H3,H4,H5,H6
      real*8 pl,x,y,xh,yh,sqr2,th,thr,rpi
      rpi=acos(-1.)
      
      x=cos(thr)
      y=sin(thr)
      xh=cos(thr/2.)
      yh=sin(thr/2.)
      sqr2=sqrt(2.)
      
      do i=1,6
      fc(i)=0.
      enddo
      
      do l=0,lmx
      fc(1)=fc(1)+(l*Mp(l)+Ep(l))*pl(x,l+1,1)
     > + ((l+1)*Mm(l)+Em(l))*pl(x,l-1,1)
     
      fc(5)=fc(5)+((l+1)*Lp(l)*pl(x,l+1,1)-l*Lm(l)*pl(x,l-1,1))
    
      fc(2)=fc(2)+((l+1)*Mp(l)+l*Mm(l))*pl(x,l,1)
      
      fc(3)=fc(3)+(Ep(l)-Mp(l))*pl(x,l+1,2)
     > + (Em(l)+Mm(l))*pl(x,l-1,2)
     
      fc(6)=fc(6)+(l*Lm(l)-(l+1)*Lp(l))*pl(x,l,1)
     
      fc(4)=fc(4)+(Mp(l)-Ep(l)-Mm(l)-Em(l))*pl(x,l,2)
      enddo
      !print*,fc(6)
      
      H1=-(1./sqr2)*y*xh*(fc(3)+fc(4))
      H2=sqr2*xh*(fc(2)-fc(1)+(fc(3)-fc(4))*yh**2)
      H3=(1./sqr2)*y*yh*(fc(3)-fc(4))
      H4=sqr2*yh*(fc(1)+fc(2)+(fc(3)+fc(4))*xh**2)
      H5=xh*(fc(5)+fc(6))
      H6=-yh*(fc(5)-fc(6))
      
      end
      
      !===============================================================
        
            
      function pl(x,n,nd)
      !---------------------------------------------------------------
      ! calculates Legendre polynomials Pn(x), as well as first and 
      ! second derivatives using the recurrence relation
      ! if n > 100 the function retuns 0.0
      !---------------------------------------------------------------
      implicit none
      double precision pl
      double precision x
      double precision pln(0:n),dpln(0:n),ddpln(0:n)
      integer n, k, nd
      
      pln(0) = 1.0
      pln(1) = x
      dpln(0)=0.
      dpln(1)=1.
      ddpln(0)=0.
      ddpln(1)=0.
      
      if (n <= 1) then
	if (nd.eq.0) pl = pln(n)
	if (nd.eq.1) pl = dpln(n)
	if (nd.eq.2) pl = ddpln(n)
	else
	do k=1,n-1
	pln(k+1) = ((2.0*k+1.0)*x*pln(k) -
     >       float(k)*pln(k-1))/(float(k+1))
	dpln(k+1)= float(k+1)*pln(k)+x*dpln(k)
	ddpln(k+1)=float(k+2)*dpln(k)+x*ddpln(k)
      end do
      if (nd.eq.0) pl = pln(n)
      if (nd.eq.1) pl = dpln(n)
	if (nd.eq.2) pl = ddpln(n)
	end if
	return
      end function
      
	!===============================================================
	  
	  
      subroutine get_M(q2,wt,parf,m_gam,lmx)
      !---------------------------------------------------------------
      ! Creates  a general matrix of all multipoles for a given 
      ! virtuality Q^2 ( in GeV^2) and center of mass energy ( in MeV) 
      !---------------------------------------------------------------
	implicit none
      complex*8 a_gami(3,5,260,4,13,2),TGa(13,13,5,4,260,3,2),
     > m_gam(4,5,13,3,2),ffit
      real*8 cme(260),q2,wt,parf(*),parv(100),chsq
      Integer Inp_indx(4),ic1,kjt,ivt,iemt,maxp,il,ilt,lmx
      integer ic,kj,iv,iem,iln,par_cnt(13,5,4,3),maxe,ilc 
      common /mtrx/a_gami,TGa,cme,par_cnt
      
Cf2py intent(out) m_gam      
          
	m_gam=0.
	parv=0.
	ffit=0.
          
	Inp_indx(1)=1
	
	ic1=Inp_indx(1)
              
	il=0
            
	do kjt=1,5
	Inp_indx(2)=kjt
	do ivt=1,4
	Inp_indx(3)=ivt
	ilt=ivt+kjt-3-(2*ivt-5)/abs(2*ivt-5)
	if(ilt.le.lmx)then
	do iemt=1,3
	Inp_indx(4)=iemt
            
	maxp=par_cnt(ic1,kjt,ivt,iemt)
            
	if(iemt.eq.1)maxe=maxp
            
	if((kjt.eq.1 .and. ivt.eq.1 .and. iemt.eq.2).or.
     >  (kjt.eq.1 .and. ivt.eq.2 .and. iemt.eq.1).or.
     >  (kjt.eq.1 .and. ivt.eq.3 .and. iemt.eq.2).or.
     >  (kjt.eq.1 .and. ivt.eq.4 .and. iemt.eq.1))then
     
	if(iemt.eq.1)parv(1:maxe)=0.
	goto 112
	endif
          
	if(iemt.eq.1)parv(1:maxe)=parf(il+1:il+maxp)
          
      if(iemt.eq.3)then
	if(ilt.eq.1 .and. kjt.eq.1 )then
	call funfit(Inp_indx,q2,wt,parf(il+1:il+maxp),ffit)
	else
	Inp_indx(4)=1
									
	call
     >  Mlong(Inp_indx,q2,wt,parv(1:maxe),
     >         parf(il+1:il+maxp),ffit,maxp,ilc)
      !Inp_indx,q2,wt,par,maxe,park,ffit,mxp,ilt
      endif
	else
	call funfit(Inp_indx,q2,wt,parf(il+1:il+maxp),ffit)
	endif
	m_gam(ivt,kjt,ic1,iemt,1)=ffit
	m_gam(ivt,kjt,ic1,iemt,2)=complex(maxp,0)
	il=il+maxp
	
112   enddo
	endif
	enddo
	enddo
				
666   end

	!===============================================================  
           
	subroutine funfit(Inp_indx,q2,wt,par,ffit)
        
	implicit none
	complex*8 a_gami(3,5,260,4,13,2),TGa(13,13,5,4,260,3,2), 
     >          Vg,TGVg,TGag(13,3),ag(3),ffit 
	real*8 q2,wt, Jt,Ist,cme(260),ed,Gem(13,3),conv,rpi,Gd,fitpd,
     >       fitfun,wtd,qpi,kg,wd,mp,mpi,Qpl,Qmn,rk,kgr,dk
	real*8 fv,fv2,dm,sm,v1,v2,kfac,rf,rfr 
	integer Inp_indx(4),ic,kj,iv,iem,ie,ict,kjt,ivt,iemt,iet,j,npar
      integer ipar,ir,jr,il,ilt,jlp,jlm,par_cnt(13,5,4,3) 
	real*8 par(*)
	real*8 mBW_F
        
      common /mtrx/a_gami,TGa,cme,par_cnt
        
	wtd=wt/1000.d0
        
	mp=0.938272d0
	mpi=0.1396d0
	sm=(wtd+mp)
	dm=(wtd-mp)
	Qpl=sqrt(sm**2+q2)
	Qmn=sqrt(dm**2+q2)
	qpi=(dm**2-mpi**2)*(sm**2-mpi**2)
	qpi=sqrt(qpi)/(2.*wtd)
	kgr=sm*dm/(2.*wtd)
	kg=Qpl*Qmn/(2.*wtd)
	rk=kg/kgr
	dk=(kg-kgr)/kgr
	v1=dm**2*(rk-1.)
	v2=q2/((rk*dm)**2/q2-1.)
	rf=sqrt(9. +3.*kg)
	rfr=sqrt(9. +3.*kgr)
      !---------------------------------------------------------------
        
	ict=Inp_indx(1)    ! Channel        
	kjt=Inp_indx(2)  ! Total angular momentum index         
	ivt=Inp_indx(3)       
	iemt=Inp_indx(4)  ! mode electric =1, magnetic =2

	ilt=ivt+kjt-3-(2*ivt-5)/abs(2*ivt-5)
        
        !jlp=-abs(1+2*ilt)+2*kjt
      jlm=-abs(1-2*ilt)+2*kjt
        
      il=ilt
	if(jlm.eq.1 .and. iemt.eq.1 .and. ilt.ge.2)then
	il=ilt-2
	endif
         
     
	!---------------------------------------------------------------
      !
      !       Scanning matrices for the requested indexes:
      !---------------------------------------------------------------
            
            
      ic=ict
      iet=0
      do ie=1,260
      ed=abs(Wt-cme(ie))
	if(ed.lt.5.)then
	iet=ie
	goto 34
	endif
	enddo
34    continue
            
      !---------------------------------------------------------------       
            
      TGag=0.
      do ic=1,13
      if(ic.lt.8 .or. ic.gt.11)then
      do j=1,3
	if(iemt.eq.3)then
	TGag(ic,j)=TGa(ict,ic,kjt,ivt,iet,j,1)
	if(kjt.eq.1 .and. (ivt.eq.2 .or. ivt.eq.4 ))
     >then
	TGag(ic,j)=TGa(ict,ic,kjt,ivt,iet,j,2)
	endif
	else
	TGag(ic,j)=TGa(ict,ic,kjt,ivt,iet,j,iemt)
	endif
	enddo
	endif
	enddo
      
      !---------------------------------------------------------------

      rpi=acos(-1.)
	conv=1000.*197.3!1./(8.883287/2./rpi)*1000.*197.3
       
      ag=0.      
      do j=1,3
	if (iemt.eq.3)then
	ag(j)=a_gami(j,kjt,iet,ivt,ict,1)
	if(kjt.eq.1 .and. (ivt.eq.2 .or. ivt.eq.4 ))
     >then
	ag(j)=a_gami(j,kjt,iet,ivt,ict,2)
	endif
	else
	ag(j)=a_gami(j,kjt,iet,ivt,ict,iemt)
	endif
	enddo
      !---------------------------------------------------------------
      !-       Q^2 parametrization
      !---------------------------------------------------------------
    
      fv=q2/(mp**2)
      fv2=v2!wt
        
      Gd=1.d0/(1.d0+q2/0.71d0)**2 !dipole formfactor
      npar=3
        
      Gem=0.
        
      do j=1,3 
	if(j.eq.1)then
	ipar=0
	do ic=1,13
	if(ic.eq.1)then
	fitpd=1.+(par(ipar+2)*fv+par(ipar+3)*fv**2)
	fitfun=exp(-par(ipar+1)*fv)*Gd*fitpd
	Gem(ic,j)=fitfun
	ipar=ipar+npar
	else
	Gem(ic,j)=Gd
	endif
      enddo
      else
      if(ag(j).ne.0.)then
      fitpd=1.+(par(ipar+2)*fv+par(ipar+3)*fv**2)
      fitfun=Gd*exp(-par(ipar+1)*fv)*fitpd
      do ic=1,13
      Gem(ic,j)=fitfun
      enddo
      ipar=ipar+npar
      endif
      endif
      enddo
        
      Vg=complex(0.d0,0.d0)
      TGVg=complex(0.d0,0.d0)
	do j=1,3
	Vg=Vg+ag(j)*Gem(1,j)
	do ic=1,13
      if(ic.lt.8 .or. ic.gt.11)then
      TGVg=TGVg+TGag(ic,j)*Gem(ic,j)
      endif
	enddo
	enddo
	kfac=rk**il*mBW_f(il,rk*par(ipar+1))
	kfac=kfac/mBW_f(il,par(ipar+1))
	
	if(iemt.eq.3)then
	ffit=(Vg+TGVg)*conv*kfac
	ffit=complex(par(ipar+2),par(ipar+3))*ffit
	else
	ffit=(Vg+TGVg)*conv*kfac
	endif
		
667   end

	!==============================================================
        
	subroutine Mlong(Inp_indx,q2,wt,par,park,ffit,mxp,ilt)
		
	implicit none
	complex*8 a_gami(3,5,260,4,13,2),TGa(13,13,5,4,260,3,2), 
     >Vg,TGVg,TGag(13,3),ag(3),ffit 
	real*8 q2,wt, Jt,Ist,cme(260),ed,Gem(13,3),conv,rpi,Gd,fitpd,
     >       fitfun,wtd,qpi,kg,wd,mp,mpi,Qpl,Qmn,rk,kgr,dk
	real*8 fv,fv2,dm,sm,v1,v2,kfac,rf,rfr 
	integer Inp_indx(4),ic,kj,iv,iem,ie,ict,kjt,ivt,iemt,iet,j,npar,
     >     ipar,ir,jr,il,ilt,jlp,jlm,par_cnt(13,5,4,3),ipark,npark,mxp
	integer maxe
	real*8 par(*),park(*)
	real*8 mBW_F
	
	common /mtrx/a_gami,TGa,cme,par_cnt
        
	wtd=wt/1000.d0
        
	mp=0.938272d0
	mpi=0.1396d0
	sm=(wtd+mp)
	dm=(wtd-mp)
	Qpl=sqrt(sm**2+q2)
	Qmn=sqrt(dm**2+q2)
	qpi=(dm**2-mpi**2)*(sm**2-mpi**2)
	qpi=sqrt(qpi)/(2.*wtd)
	kgr=sm*dm/(2.*wtd)
	kg=Qpl*Qmn/(2.*wtd)
	rk=kg/kgr
	dk=(kg-kgr)/kgr
	v1=dm**2*(rk-1.)
	v2=q2/((rk*dm)**2/q2-1.)
	rf=sqrt(9. +3.*kg)
      rfr=sqrt(9. +3.*kgr)
             
      !---------------------------------------------------------------
        
      ict=Inp_indx(1)    ! Channel        
      kjt=Inp_indx(2)  ! Total angular momentum          
	ivt=Inp_indx(3)     
	iemt=Inp_indx(4)  ! mode electric =1, magnetic =2
            
	ilt=ivt+kjt-3-(2*ivt-5)/abs(2*ivt-5)
        
	  !jlp=-abs(1+2*ilt)+2*kjt
	jlm=-abs(1-2*ilt)+2*kjt
        
      il=ilt
	if(jlm.eq.1 .and. iemt.eq.1 .and. ilt.ge.2)then
	il=ilt-2
	endif
        
	ic=ict
	iet=0
	do ie=1,260
	ed=abs(Wt-cme(ie))
	if(ed.lt.5.)then
	iet=ie
	goto 35
	endif
	enddo
35    continue
            
      !---------------------------------------------------------------
	TGag=0.     
      do ic=1,13
	if(ic.lt.8 .or. ic.gt.11)then
	do j=1,3
	TGag(ic,j)=TGa(ict,ic,kjt,ivt,iet,j,iemt)
	enddo
	endif
	enddo
      !---------------------------------------------------------------
            
      rpi=acos(-1.)
	conv=1000.*197.3!*1./(8.883287/2./rpi)!*1000.*197.3
       
	ag=0.     
	do j=1,3
	ag(j)=a_gami(j,kjt,iet,ivt,ict,iemt)
	enddo
		
	!---------------------------------------------------------------
      !-       Q^2 parametrization
      !---------------------------------------------------------------
    
      fv=q2/(mp**2)
	fv2=rk*mBW_f(1,rk)!(Qmn/mp)
        
	Gd=1.d0/(1.d0+q2/0.71d0)**2 !dipole formfactor
	npar=3
	npark=3
        
	Gem=0.
        
	do j=1,3 
	if(j.eq.1)then
	ipar=0
	ipark=0
	do ic=1,13
	!if(ic.eq.1 .or. ic.eq.5 .or. ic.eq.6 .or. ic.eq.12)then
	if(ic.eq.1)then
	fitpd=1.+(par(ipar+2)*fv+par(ipar+3)*fv**2) 
	fitfun=Gd*exp(-(par(ipar+1)*fv+park(ipark+1)*fv2))*fitpd
	Gem(ic,j)=fitfun*(1.+park(ipark+2)*fv2+park(ipark+3)*fv2**2)
	ipar=ipar+npar
	ipark=iparK+npark
	else
	Gem(ic,j)=Gd
	endif
	enddo
	else
	if(ag(j).ne.0.)then
	fitpd=1.+(par(ipar+2)*fv+par(ipar+3)*fv**2)
	fitfun=Gd*exp(-(par(ipar+1)*fv+park(ipark+1)*fv2))*fitpd
	do ic=1,13
	Gem(ic,j)=fitfun*(1.+park(ipark+2)*fv2+park(ipark+3)*fv2**2)
	enddo
	ipar=ipar+npar
	ipark=ipark+npark
	endif
	endif
	enddo
	mxp=ipark
	
	Vg=complex(0.d0,0.d0)
	TGVg=complex(0.d0,0.d0)
	do j=1,3
	Vg=Vg+ag(j)*Gem(1,j)
	if(kjt.eq.1 .and. (ivt.eq.2 .or. ivt.eq.4))then
	mxp=ipark+1
	Vg=Vg+park(mxp)*Gem(1,j)
	endif
	do ic=1,13
	if(ic.lt.8 .or. ic.gt.11)then
	TGVg=TGVg+TGag(ic,j)*Gem(ic,j)
	if(kjt.eq.1 .and. (ivt.eq.2 .or. ivt.eq.4))then
	mxp=ipark+3
	TGVg=TGVg+complex(park(mxp-1),park(mxp))*Gem(ic,j)
	endif
	endif
	enddo
	enddo

      kfac=rk**il*mBW_f(il,rk*par(ipar+1))
	kfac=kfac/mBW_f(il,par(ipar+1))
        
	if(jlm.eq.1 .and. iemt.eq.1 .and. ilt.ge.2)then
	kfac=float(1-ilt)/float(ilt)*kfac
	endif
	
      ffit=(Vg+TGVg)*conv*kfac
      ffit=ffit/(sm**2+q2)!*(kg**2-q2)/(sm**2+q2)
             
667   end

      !=============================================================
        
	function mBW_f(il,x) 
	!altered Blatt Weisskopff factors
	implicit none
	double precision mBW_f
	double precision x,ten
	integer il
	if(il.eq.0)mBW_f=1.
	if(il.eq.1)mBW_f=1./sqrt(1.+x**2)
	if(il.eq.2)mBW_f=3./sqrt(9.+3.*x**2+x**4)
	if(il.eq.3)mBW_f=15./sqrt(225.+45.*x**2+6.*x**4+x**6)
	if(il.eq.4)then
	ten=10.
	mBW_f=105./sqrt(11025.+1575.*x**2+135.*x**4+ten*x**6+x**8)
	endif
	end function 
	
      !===============================================================
          
	subroutine create(ict,a_gami,TGa,cme)
      implicit none
      integer ict
	integer i_c,k_j,i_v,i_en,i_j,ic,kj,iv,ien,iq2,j,ic1,ic2,ie_m
	real*8 cme(260),ez,ReTGa,AimTGa
	real*8 re_alph_NPe,ai_alph_NPe,re_alph_NPm,ai_alph_NPm
	complex*8 a_gami(3,5,260,4,13,2),TGa(13,13,5,4,260,3,2)

Cf2py intent(in) ict
Cf2py intent(out) a_gami,TGa,cme

	
	a_gami=0.
	TGa=0.
		
	!opening file with  photocouplings alpha_gamma
	open(unit=1,file='vgam_all4.txt')
            
	read(1,*)
	
	do i_c=1,13
	do k_j=1,5
	do i_v=1,4
	do i_en=1,260 ! Energy index
	do i_j=1,3 ! non-pole (j=1) and pole indexes(j=2,3)
                    
	read(1,*,end=111)ic,kj,iv,ien,iq2,ez,j 
     >    ,re_alph_NPe,ai_alph_NPe,re_alph_NPm,ai_alph_NPm
                 
	a_gami(j,kj,ien,iv,ic,1)=complex(re_alph_NPe,ai_alph_NPe)
	a_gami(j,kj,ien,iv,ic,2)=complex(re_alph_NPm,ai_alph_NPm)
	
	cme(ien)=ez
	
	enddo
	enddo
	enddo
	enddo
	enddo
            
111   close(unit=1)
    
      !---------------------------------------------------------------
    
      !Calling file with rescattering term, TGa=T_transport * G *alpha_gamma 
             
      open(unit=2, file='rescatt_mat2.txt')
            
	TGa=0.
            
	read(2,*)
            
	do 
	read(2,*,end=112)ic1,ic2,kj,iv,ien,j,ie_m,ReTGa,AimTGa
	!print*,ReTGa,AimTGa
	if(ic1.eq.ict+1)exit
	TGa(ic1,ic2,kj,iv,ien,j,ie_m)=complex(ReTGa,AimTGa)
	enddo
	
112   close(unit=2)

      !---------------------------------------------------------------
	end
	  
      !===============================================================
  
	!program Test_count
	!integer par_cnt(13,5,4,3)
	!
	!call par_count(par_cnt)
	!
	!
	!ic=1
	!do kj=1,3
	!do iv=1,4
	!do iem=1,2
	!print*,par_cnt(ic,kj,iv,iem)
	!enddo
	!enddo
	!enddo
	!
	!end
	
	subroutine par_count(par_cnt)
	
	implicit none
	real*8 parv(12),wt,chsq
	integer ic,kj,iv,iem,npar,par_cnt(13,5,4,3)
	  
Cf2py intent(out) par_cnt
	  
	parv=0.
	par_cnt=0
	
	open(unit=1, file='fill_par.txt')
	open(unit=2,file='fill_par2.txt')
	do
	read(1,*,end=100)ic,wt,kj,iv,iem,chsq,npar,parv(1:npar)
	
	par_cnt(ic,kj,iv,iem)=npar
	  
	enddo
	  
100   close(1)

	do
	read(2,*,end=113)ic,wt,kj,iv,iem,chsq,npar,parv(1:npar)
	if(kj.eq.1 .and. iv.eq.2)goto114
	if(kj.eq.1 .and. iv.eq.4)goto114
	par_cnt(ic,kj,iv,iem)=npar
114   enddo
        
113   close(2)

	end
	  
	!===============================================================
				
	subroutine creater(ict,ragami,iagami,rTGa,iTGa,cmet)
	implicit none
	integer ict
	integer i_c,k_j,i_v,i_en,i_j,ic,kj,iv,ien,iq2,j,ic1,ic2,ie_m
	real*8 cmet(260),ez,ReTGa,AimTGa
	real*8 re_alph_NPe,ai_alph_NPe,re_alph_NPm,ai_alph_NPm
	real*8 ragami(3,5,260,4,13,2),rTGa(13,13,5,4,260,3,2)
	real*8 iagami(3,5,260,4,13,2),iTGa(13,13,5,4,260,3,2)
	
Cf2py intent(in) ict
Cf2py intent(out) ragami,iagami,rTGa,iTGa,cmet
			
	ragami=0.
	iagami=0.
	rTGa=0.
	iTGa=0.
				
	!opening file with  photocouplings alpha_gamma
	open(unit=1,file='vgam_all4.txt')
				  
	read(1,*)
				  
	do i_c=1,13
	do k_j=1,5
	do i_v=1,4
	do i_en=1,260 ! Energy index
	do i_j=1,3 ! non-pole (j=1) and pole indexes(j=2,3)
						  
	read(1,*,end=211)ic,kj,iv,ien,iq2,ez,j 
     >   ,re_alph_NPe,ai_alph_NPe,re_alph_NPm,ai_alph_NPm
					   
	ragami(j,kj,ien,iv,ic,1)=re_alph_NPe
	ragami(j,kj,ien,iv,ic,2)=re_alph_NPm
	iagami(j,kj,ien,iv,ic,1)=ai_alph_NPe
	iagami(j,kj,ien,iv,ic,2)=ai_alph_NPm
				
	cmet(ien)=ez
				  
	enddo
	enddo
	enddo
	enddo
	enddo
				  
211   close(unit=1)
	  
	!---------------------------------------------------------------
		  
	!Calling file with rescattering term, TGa=T_transport * G *alpha_gamma 
				   
	open(unit=2, file='rescatt_mat2.txt')				  
				
	read(2,*)
				  
	do 
	read(2,*,end=212)ic1,ic2,kj,iv,ien,j,ie_m,ReTGa,AimTGa
	!print*,ReTGa,AimTGa
	if(ic1.eq.ict+1)exit
	rTGa(ic1,ic2,kj,iv,ien,j,ie_m)=ReTGa
	iTGa(ic1,ic2,kj,iv,ien,j,ie_m)=AimTGa
				
	enddo
				  
212   close(unit=2)
			
	!---------------------------------------------------------------
	end

	!===============================================================
C END OF FILE obs_py.f
	  





