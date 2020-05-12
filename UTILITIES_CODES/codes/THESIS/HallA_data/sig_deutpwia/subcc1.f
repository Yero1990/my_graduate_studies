c
c
c
c	subroutine subcc1
c
c	author:	R.W. Lourie
c	date:	V1.0	18-AUG-1986
c	modifications:	W. Boeglin
c	        Instead of angles , the four vectors of the incident 
c		beam, scattered electron and the outgoing proton are
c		entered in the argument list.
c		floor: x,z plane
c		       y-axis points out if the plane
c		       x-axis points to the electron side
c		       z-axis points in the beam direction
c	purpose:
C		Calculates elementary OFF-SHELL (e,p) cross-section
C		using the prescription given by DeForest in
C		Nucl. Phys.  A392 (1983)
C
c	variables:
c	input:
C		k_i : incident 4-vector
c		k_f : scattered 4-vector
c		p_f : outgoing proton 4-vector
C	OUTPUT:
C		sigma_ep: The off-shell cross section.
C		k_fact:   Kinematic factor ( = E_p*P_F)
C		pi_mag:   Magnitude of initial momentum (determined as
C			  PF - Q componentwise).
C
C	ROUTINES CALLED:  
C		form_fact
C
C
        subroutine subcc1 (k_i,k_f,p_f,  sigma_ep,k_fact,pi_mag,nucleon)
	implicit real*8 (a-h,k,m,o-z)
	real*8 k_i(4),k_f(4),q(4),p_i(4),p_f(4),n_e(4),n_p(4)
	character*7 fftype, nucleon*1
c Proton mass in MeV
	parameter (mp = 938.28)	
	parameter (pi = 3.14159)
	parameter (hbarc = 197.3286)
c fine structure constant
	parameter (alpha = 7.297e-3)		
c Proton anomalous moment
	parameter (kappa = 1.793)			
	data a,z,m_eff /1.,1.,1./
c	data fftype /'mainz  '/

c use the John Arrington fit as in SIMC
	data fftype /'jasimc '/
c
        dtr=pi/180.
        e_p_f=p_f(4)
        e_f=k_f(4)
c
C	Form components of momentum transfer 4-vector.
C	Form components of proton initial 4-momentum vector.
C
	do  10 i=1,4
	   q(i) = k_i(i) - k_f(i)
	   p_i(i) = p_f(i) - q(i)
10	continue
c
c vector normal to electron scattering plane
c
        call cross(k_i,k_f,n_e)
c
c vector normal to q p_f plane
c
        call cross(q,p_f,n_p)
c
c calculate out of plane angle phi_p_f
c
        rn_e=sqrt(dot3(n_e,n_e))
        rn_p=sqrt(dot3(n_p,n_p))
        cosphi=dot3(n_e,n_p)/(rn_e*rn_p)
        phi_p_f=acos(cosphi)
c
C	Form all dot products.
C
	qmu2  = -dot(q,q)
        pf=sqrt(p_f(4)**2-dot(p_f,p_f))
c
c calculate the electron scattering angle
c
        sin_phi=sqrt(qmu2/(4.*k_i(4)*k_f(4)))
        phi_e=2.*asin(sin_phi)
c
	pf_dot_q = 0.0
	pi_mag = 0.0
	q_mag = 0.0
	do 20 i=1,3
	   pi_mag = pi_mag + p_i(i)**2
	   q_mag = q_mag + q(i)**2
	   pf_dot_q = pf_dot_q + p_f(i)*q(i)
20	continue
	pi_mag = sqrt(pi_mag)
	q_mag = sqrt(q_mag)
c Angle between P_F and Q
	cosgam = (pf_dot_q/(pf*q_mag))    
	gamma = aacos(cosgam)
c
C	Set up kinematic factors etc.
C
c Initial proton energy
	e_bar = sqrt(pi_mag**2+mp**2)		
	omega_bar = e_p_f - e_bar
	q_mu_bar = sqrt(q_mag**2 - omega_bar**2)
	ebp = e_bar*e_p_f
	x = q_mu_bar**2/4./mp**2
	q_mu = (abs(qmu2))**.5
	qq = (q_mu/q_mag)**2
	tan_e2 = tan(phi_e/2.)**2
C
C	Put together various response functions
C
	call form_fact(fftype,a,z,q_mu,m_eff,ge2,
     #			     gm2,f1p,f2p,f1n,f2n)
c 
c select the proton or neutron form factors 
c
        if (nucleon.eq.'p') then
        
	   w_c = ((e_bar+e_p_f)**2*(f1p**2+kappa**2*x*f2p**2)-q_mag**2
     #          * (f1p+kappa*f2p)**2)/(4.*ebp)
c
	   w_t = (f1p+kappa*f2p)**2*q_mu_bar**2/2./ebp
c
	   w_s = (f1p**2+kappa**2*x*f2p**2)*pf**2*sin(gamma)**2/ebp
c
	   w_i = -(e_bar+e_p_f)*(f1p**2+kappa**2*x*f2p**2)*pf
     #           * sin(gamma)/ebp
        else
	   w_c = ((e_bar+e_p_f)**2*(f1n**2+kappa**2*x*f2n**2)-q_mag**2
     #          * (f1n+kappa*f2n)**2)/(4.*ebp)
c
	   w_t = (f1n+kappa*f2n)**2*q_mu_bar**2/2./ebp
c
	   w_s = (f1n**2+kappa**2*x*f2n**2)*pf**2*sin(gamma)**2/ebp
c
	   w_i = -(e_bar+e_p_f)*(f1n**2+kappa**2*x*f2n**2)*pf
     #           * sin(gamma)/ebp
        endif
C
C	Calculate Off-shell cross-section.
C
	sigma_mott = 4.*((alpha*cos(phi_e/2.)*e_f)/qmu2)**2
	sigma_ep = sigma_mott*(qq**2*w_c + (qq/2.+tan_e2)*w_t
     #           + qq*sqrt(qq+tan_e2)*w_i*cos(phi_p_f)
     #           + (qq*cos(phi_p_f)**2+tan_e2)*w_s)*hbarc**2.*10000.

	k_fact = e_p_f*pf
c
        return
	end
C
C	Function DOT(A,B)
C
C	Purpose:
C
C		Takes dot product of two 4-vectors.
C		Metric:  DOT(A,A) = A_time**2 - A_space .dot. A_space
C
C
	function dot(a,b)
	implicit real*8 (a-h,k,o-z)
	real*8 a(4),b(4)
	space = 0.
	do 10 i = 1,3
	   space = space + a(i)*b(i)
10	continue
	time = a(4)*b(4)
	dot = time - space
	return
	end
c
C	SUBROUTINE CROSS(A,B,C)
C
C	Purpose:
C
C		Takes the cross product of the space part of two 
C	        4-vectors.
C		CROSS(A,B,C): C = A X B
C
C
	subroutine cross(a,b,c)
	implicit real*8 (a-h,k,o-z)
	real*8 a(4),b(4),c(4)
c
        c(1)=a(2)*b(3)-a(3)*b(2)
        c(2)=a(3)*b(1)-a(1)*b(3)
        c(3)=a(1)*b(2)-a(2)*b(1)
        c(4)=0.
c
	return
	end
c
C	Function DOT3(A,B)
C
C	Purpose:
C
C		Takes 3-dot product of two 4-vectors.
C
C
	function dot3(a,b)
	implicit real*8 (a-h,k,o-z)
	real*8 a(4),b(4)
	space = 0.
	do 10 i = 1,3
	   space = space + a(i)*b(i)
10	continue
	dot3 = space
	return
	end
c
c
C
C -------------------------------------------------------------------
C     Deuteron momentum distribution.
C
C     Momenta are in GeV/c.
C     Constants are from Krautschneider's Ph.D. thesis
C -------------------------------------------------------------------
C
      double precision function h2spec(prmag)
      implicit real*8 (a-h,o-z)
      parameter (pi = 3.14159)
      if(prmag .le. 300.)then
c No rescattering.
        term = 0.         
      else
	term = 4.
      endif
c Convert to GeV/c
      pr = prmag/1000.  
      pr2 = pr**2
      t1 = pr2 + 0.002088
      t2 = pr2 + 0.0676
c
C     The constant TERM from K. Mueller Ph.D. thesis simulates
C     rescattering contribution for high
C     momenta to reach agreement with experimental data.
C     Normal Value:  TERM=4.
c
      h2spec = (1./t1-1./t2)**2 + term
      h2spec = 4.*pi*h2spec * 1.0e-12
      return
      end
c
c
c       subroutine form_fact
C
C
C	AUTHOR:	B.H.Cottman
C	DATE:	V1.0	23-SEP-83
C	MODIFICATIONS:	
C		Lots of additions by RWL: Janssens, Blatnik and Zovko,
C		and Gari/Krumpelmann formfactors.
C	PURPOSE:
C		Calculate dipole, MAINZ fit or FIT 8.2 mean formfactors.
C	VARIABLES:
C	INPUT:
C		a:	Number of nucleons
C		z:	Number of protons
C		q_mu:	4-momentum
C		m_eff:	Factor to scale nucleon mass.	
C		
C	OUTPUT:
C		ge2:	Electric formfactor squared	
C		gm2:	Magnetic formfactor squared	
C		f1p,n:	Dirac formfactor for proton and neutron
C		f2p,n:	Pauli formactor for proton and neutron
C	ROUTINES CALLED:  
C			None
C
C added parameterizations for John Arrington and P.Posted parameterization
C used in SIMC
C
	subroutine form_fact(fftype,a,z,q_mu,m_eff,ge2,gm2,
     #			     f1p,f2p,f1n,f2n)
c
	implicit real*8 (a-h,k,m,o-z)
	character*7 fftype
	parameter (pi = 3.1415927)
	parameter (mp = 938.2796)
	parameter (mn = 939.5655)
	parameter (hbarc = 197.3286)
c
C  dipole fit
	gdp(qq)=1./(1.+qq/0.71e+6)**2				
c
C mainz fit
	g_e_m(t)=0.312/(1+t/6.0)+				
     &		     1.312/(1+t/15.02)+
     &               (-0.709)/(1+t/44.08)+
     &               0.085/(1+t/154.2)
	g_m_m(t)=( 0.694/(1+t/8.5)+
     &		       0.719/(1+t/15.02)+
     &                 (-0.418)/(1+t/44.08)+
     &                 0.005/(1+t/355.4) )
c
	q2 = q_mu**2
	q2_f=q2/hbarc**2
C neutron mag. moment
	mun=-1.913					
C  proton mag. moment
	mup=2.793					
	kappan=mun
	kappap=mup-1.
	mpi = 0.1395
	x=1.+q2/4./mp**2
c
C	Dipole fit
C
	if(fftype .eq. 'dipole ') then
C  dipole chosen
	 f1p=(1.+q2*mup/4./mp**2)*gdp(q2)/x		
	 f1n=(q2*mun/2./mp**2)*gdp(q2)/x
	 f2p=gdp(q2)/x
	 f2n=(2.-x)*gdp(q2)/x
	elseif(fftype .eq. 'mainz  ') then
c
C mainz fit
c
	 f1p=(g_e_m(q2_f)+(x-1.)*mup*g_m_m(q2_f))/x		
	 f1n=(q2*mun/2./mp**2)*g_m_m(q2_f)/x
	 f2p=(mup*g_m_m(q2_f)-g_e_m(q2_f))/(kappap*x)
	 f2n=(2.-x)*g_m_m(q2_f)/x
	elseif(fftype .eq. 'fit8_2 ') then
c
C fit 8.2 chosen
c
	 t=-q2/1.e+6/MPI**2					
	 f1rho = (.4775+.045/(1-t*.01948/.355)**2) / (1-t*.01948/.536)
	 f2rho = (2.6675+.481/(1-t*.01948/.268)) / (1-t*.01948/.603)
	 f1v = f1rho + (0.03144/(1-t/75.736)) - (0.08575/(1-t/308.179))
     1   + (.03164/(1-t/446.642))
	 f1s = (-.61668/(1-t/53.37)) + (1.15702/(1-t/31.468))
     1   - (.040333/(1-t/166.648))
	 f2v = f2rho - (1.3512/(1-t/75.736)) + (.03295/(1-t/308.179))
     1   + (.02151/(1-t/446.642))
	 f2s = (.12097/(1-t/53.37)) + (-.17534/(1-t/31.468))
     1   - (.00363/(1-t/166.648))
	 f1p = (f1v+f1s)
	 f1n = f1s-f1v
	 f2p = (f2v+f2s)/kappap
	 f2n = (f2s-f2v)/kappan
	elseif(fftype .eq. 'gari   ')then
c
c Gari
c
	  call gari(q2/1.0e+06,f1p,f2p,f1n,f2n,gep,gmp,gen,gmn)
	  goto 777
	elseif(fftype .eq. 'blatnik')then
c
c Blatnik
c
	  call b_and_z(-q2/1.0e+06,f1p,f2p,f1n,f2n,gep,gmp,gen,gmn)
	  goto 777
	elseif(fftype .eq. 'janssen')then
c
c Janssen
c
	  call janssens(q2,f1p,f2p,f1n,f2n,gep,gmp,gen,gmn)
	  goto 777
	elseif(fftype .eq. 'jasimc ')then
c
c J. Arrington from SIMC
c
	  call jasimc(q2,f1p,f2p,f1n,f2n,gep,gmp,gen,gmn)
	  goto 777
	elseif(fftype .eq. 'bosted ')then
c
c P. Bosted from SIMC
c
	  call bosted(q2,f1p,f2p,f1n,f2n,gep,gmp,gen,gmn)
	  goto 777
	endif
c
	xx=  q2/4.0/mn/mn
C  neutron and proton...
	gen=f1n-xx*kappan*f2n			
C  elec. and mag. form..
	gep=f1p-xx*kappap*f2p			
C  factors
	gmn=f1n+kappan*f2n			
	gmp=f1p+kappap*f2p
C  nuclear elec. ff per nucleon
777	ge2=(z*gep**2+(a-z)*gen**2)		
C nuclear mag. ff per nucleon
	gm2=(z*gmp**2+(a-z)*gmn**2)		
c
	return
	end
c
	subroutine gari(q2,f1p,f2p,f1n,f2n,gep,gmp,gen,gmn)
c
c   program by gond yen, november 20, 1986
c
C   nucleon form factor parametrization by M. Gari and
C   W. Kruempelmann:
C   --- Z. Phys. A322(1985)689-693
C       Phys. Lett. B173(1986)10-14
C
C   the unit used in this code is GeV
C
C   dipole form: FD
C   FD=1.0/(Q2+0.71)**2
C
      implicit real*8 (a-h,l,k,m,o-z)
      data mn,mr,mw /0.938e0,0.776e0,0.784e0/
      data k,kv,ks,kr,kw /1.793e0,3.706e0,-0.12e0,6.62e0,0.163e0/
      data cr,cw /0.377e0,0.411e0/
      data lam1,lam2,lmqcd /0.795e0,2.27e0,0.29e0/
      mr2=mr**2
      mw2=mw**2
      mn2=mn**2
      lm1SQ=lam1**2
      lm2SQ=lam2**2
      lmqSQ=lmqcd**2
      qt2=q2*dlog((lm2SQ+q2)/lmqSQ)/dlog(lm2SQ/lmqSQ)
      f1=lm1SQ/(lm1SQ+qt2)*lm2SQ/(lm2SQ+qt2)
      f2=lm1SQ/(lm1SQ+qt2)*(lm2SQ/(lm2SQ+qt2))**2
      f1v=f1*(cr*mr2/(mr2+q2)+1.0-cr)
      f1s=f1*(cw*mw2/(mw2+q2)+1.0-cw)
      kf2v=f2*(cr*kr*mr2/(mr2+q2)+kv-cr*kr)
      kf2s=f2*(cw*kw*mw2/(mw2+q2)+ks-cw*kw)
C  Dirac and Pauli form factors: F1, F2's
      f1n=0.5*(f1s-f1v)
      f2n=0.5*(kf2s-kf2v)/((ks-kv)/2.)
      f1p=0.5*(f1s+f1v)
      f2p=0.5*(kf2s+kf2v)/k
C  Sachs form factors: GE, GM's
      gmn=f1n+(ks-kv)/2.*f2n
      gen=f1n-(ks-kv)/2.*q2/4.0/mn2*f2n
      gmp=f1p+k*f2p
      gep=f1p-k*q2/4.0/mn2*f2p
c
      return
      end

	subroutine b_and_z(q2,f1p,f2p,f1n,f2n,gep,gmp,gen,gmn)
	implicit real*8 (a-h,k,m,o-z)
	data bs,bv,mu_s,mu_v,mp /-0.91,-1.10,-0.06,1.853,0.93828/
	data t_r,t_o,t_p,t_rp,t_rpp,t_op /0.585,0.614,1.039,1.3,
     #	     2.1,1.4/

c
	r_s(t) = 1./((1.-t/t_o)*(1.-t/t_p)*(1.-t/t_op))	
	r_v(t) = 1./((1.-t/t_r)*(1.-t/t_rp)*(1.-t/t_rpp))	
c
	x = q2/4./mp**2
	mu_p = mu_s + mu_v + 1.
	mu_n = mu_s - mu_v
	gep = (0.5 + x*(mu_s + 2.*mp**2*bs))*r_s(q2)
     #	    + (0.5 + x*(mu_v + 2.*mp**2*bv))*r_v(q2)
	gen = (0.5 + x*(mu_s + 2.*mp**2*bs))*r_s(q2)
     #	    - (0.5 + x*(mu_v + 2.*mp**2*bv))*r_v(q2)
	gmp = (0.5 + mu_s + bs*q2/2.)*r_s(q2)
     #      + (0.5 + mu_v + bv*q2/2.)*r_v(q2)
	gmn = (0.5 + mu_s + bs*q2/2.)*r_s(q2)
     #      - (0.5 + mu_v + bv*q2/2.)*r_v(q2)
	f1p = (gep - x*gmp)/(1.-x)
	f2p = (gmp - gep)/(mu_p-1.)/(1.-x)
	f1n = (gen - x*gmn)/(1.-x)
	f2n = (gmn - gen)/mu_n/(1.-x)
	return
	end
c
	subroutine janssens(q2,f1p,f2p,f1n,f2n,gep,gmp,gen,gmn)
	implicit real*8 (a-h,k,m,o-z)
	data mp,mu_p,mu_n,hbarc /938.28,2.793,-1.913,197.3286/
c
	ge_s(t) = (2.5/(1.+t/15.7)-1.6/(1.+t/26.7)+0.1)*0.5
	gm_s(t) = (3.33/(1.+t/15.7)-2.77/(1.+t/26.7)+0.44)*0.44
	ge_v(t) = (1.16/(1.+t/8.19)-0.16)*0.5
	gm_v(t) = (1.11/(1.+t/8.19)-0.11)*2.353
c
	tt = q2/hbarc/hbarc
	gep = ge_s(tt) + ge_v(tt)
	gen = ge_s(tt) - ge_v(tt)
	gmp = gm_s(tt) + gm_v(tt)
	gmn = gm_s(tt) - gm_v(tt)
c
	x = q2/4./mp**2
	f1p = (gep+x*gmp)/(1.+x)
	f2p = (gmp-gep)/(1.+x)/(mu_p-1.)
	f1n = (gen+x*gmn)/(1.+x)
	f2n = (gmn-gen)/(1.+x)/mu_n
	return
	end
c

C -------------------------------------------------------------------
C     Now code in the John Arrington model from SIMC for comparison
C-------------------------------------------------------------------
      subroutine jasimc(q2,f1p,f2p,f1n,f2n,gep,gmp,gen,gmn)
      implicit none

      real*8 q2,f1p,f2p,f1n,f2n,gep,gmp,gen,gmn
      real *8 pi, mp, mn, hbarc, gdp, x, xx
      real*8 mun, mup, kappan, kappap, tau, qq
      real*8 q_mu_gev, q_mu2_gev

      parameter (pi = 3.1415927)
      parameter (mp = 938.2796)
      parameter (mn = 939.5655)
      parameter (hbarc = 197.3286)

      gdp(qq)=1./(1.+qq/0.71e+6)**2                           ! dipole fit

c neutron magnetic moment
      mun=-1.913    
c proton magnetic moment
      mup=2.793       

      kappan=mun
      kappap=mup-1.
      tau=q2/4./mp**2
      x=1.+tau
      xx=  q2/4.0/mn/mn

c qmu in gev/c
      q_mu2_gev = q2*1.d-6


c JohnA model from SIMC
        f1n=5.6*tau*tau*mun*gdp(q2)/((1.+tau)*(1.+5.6*tau))
        f2n=(1.+6.6*tau)*gdp(q2)/((1.+tau)*(1.+5.6*tau))
        gen=f1n-xx*kappan*f2n
        gmn=f1n+kappan*f2n      ! keep the neutron part same as above
c 
c       -- this calculates the form factors gep and gmp using
c       john arrington's fit to world data 
c       reference: phys. rev. c69, 022201, 2004
c
        gmp=mup/(1.d0
     &            + q_mu2_gev * (3.19)
     &            + q_mu2_gev**2 * (1.355)
     &            + q_mu2_gev**3 * (0.151)
     &            + q_mu2_gev**4 * (-0.114e-01)
     &            + q_mu2_gev**5 * (0.533e-03)
     &            + q_mu2_gev**6 * (-0.900e-05)  )
        gep=1.d0/(1.d0
     &            + q_mu2_gev * (3.226)
     &            + q_mu2_gev**2 * (1.508)
     &            + q_mu2_gev**3 * (-0.3773)
     &            + q_mu2_gev**4 * (0.611)
     &            + q_mu2_gev**5 * (-0.1853)
     &            + q_mu2_gev**6 * (0.1596e-01)  )
        f1p=(gep+xx*gmp)/x
        f2p=(gmp-gep)/(kappap*x)
        end

c ---------------------------------------------------------------
C   Add in the Peter Bosted model from SIMC as well
C --------------------------------------------------------------
      subroutine bosted()
      implicit none

      real*8 q2,f1p,f2p,f1n,f2n,gep,gmp,gen,gmn
      real *8 pi, mp, mn, hbarc, gdp, x, xx
      real*8 mun, mup, kappan, kappap, tau
      real*8 q_mu_gev, q_mu2_gev, qq

      parameter (pi = 3.1415927)
      parameter (mp = 938.2796)
      parameter (mn = 939.5655)
      parameter (hbarc = 197.3286)
C dipole fit
      gdp(qq)=1./(1.+qq/0.71e+6)**2 

c neutron magnetic moment
      mun=-1.913    
c proton magnetic moment
      mup=2.793       

      kappan=mun
      kappap=mup-1.
      tau=q2/4./mp**2
      x=1.+tau
      xx=  q2/4.0/mn/mn

c qmu in gev/c
      q_mu2_gev = q2*1.d-6
      q_mu_gev = sqrt(q2)*1.d-3


c Peter Bosted  model from SIMC
        f1n=5.6*tau*tau*mun*gdp(q2)/((1.+tau)*(1.+5.6*tau))
        f2n=(1.+6.6*tau)*gdp(q2)/((1.+tau)*(1.+5.6*tau))
        gen=f1n-xx*kappan*f2n
        gmn=f1n+kappan*f2n      ! keep the neutron part same as above
c
c       this calculates the gep and gmp form factors using the fit 
c       to world data by peter bosted.
c       reference: phys. rev. c51, 409, 1995
c
        gmp=mup/(1.d0
     &            + q_mu_gev * (0.35)
     &            + q_mu_gev**2 * (2.44)
     &            + q_mu_gev**3 * (0.5)
     &            + q_mu_gev**4 * (1.04)
     &            + q_mu_gev**5 * (0.34) )
        gep=1.d0/(1.d0
     &            + q_mu_gev * (0.62)
     &            + q_mu_gev**2 * (0.68)
     &            + q_mu_gev**3 * (2.8)
     &            + q_mu_gev**4 * (0.83) )
        f1p=(gep+xx*gmp)/x
        f2p=(gmp-gep)/(kappap*x)
c ..................................
        return 
        end


	function aacos(arg)
	implicit real*8 (a-h,o-z)
c
c this function calculates the acos and checkes whether the argument is
c bigger than 1. If it is the function is 1. 
	
	if (dabs(arg).gt.1.d0) then
	   aacos=0.d0
	else
	   aacos=acos(arg)
	endif
	return
	end



