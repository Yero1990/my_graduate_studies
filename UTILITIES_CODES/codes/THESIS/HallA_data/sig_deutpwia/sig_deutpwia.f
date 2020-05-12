c this program calculates the ee'p cross section in the PWIA for
c deuterium. 
c the parameterization of the paris potential momentum distribution
c
c this is the same basic code as the one for the deutpwia1 and deutpwia1_sub

      real*8 function sig_deutpwia(
c input
     >     e0, omega, thee,
     >     theta_p, phi_p, angle_type )


      implicit real*8 (a-h,o-z)
      real*8 md,mn,mp,mnp,klab,ml ,kfact,k1,k2
      real*8 k_i(4), k_f(4), p_f(4)
      logical first,lerr
      integer angle_type, spec_pos
      character nucleon*1,cread*1

c output
      common /results/
     >     rho, kfact, sigep, sigma, f_rec,
     >     thq, ppi, ppf, enp, qc, cthnpc, qlab
      
      data amu/931.48/,hbarc/197.33/,pi/3.1415927/,first/.true./
      data rpa/1.75d-5/

      data spec_pos/1/

      md=2.014*amu 
      mn=1.0087*amu 
      mp=1.0073*amu 
      dtr=pi/180. 
      rtd=180./pi 
c
      einc=e0
      thee=thee*dtr
      nucleon = 'p'

c
c kinematic for deuterium
c
      k1=einc 
      k2=(einc-omega) 
c 
      qnue2=4*k1*k2*(sin(thee/2))**2 
      qlab2=k1**2+k2**2-2*k1*k2*cos(thee) 
      qlab=sqrt(qlab2)
c theta-q
      cthq=(e0-(e0-omega)*cos(thee))/qlab
      thq=acos(cthq)

c handle angle types
      if ( angle_type .eq. spec_pos) then
         theta = theta_p
         phi = phi_p
         theta=theta*dtr
         phi=phi*dtr 
         cth=cos(phi)*cos(thq-theta)
      else
         thnp = theta_p
         phinp = phi_p
         thnp=thnp*dtr
         phinp=phinp*dtr
c
         cth=cos(thnp)
         phi=asin(sin(thnp)*sin(phinp))
         beta=acos(cth/cos(phi))
c
         diff=phinp-pi/2
c
         if (abs(diff).lt.rpa) then
            theta=thq
         elseif (diff.lt.0) then
            theta=thq-beta
         elseif (diff.gt.0) then
            theta=thq+beta
         endif
      endif
c
c
      call relkin(omega,qlab,cth, 
     &                 conv,cthnpc,enp,qc,ppf,ppi,lerr)
c
c
c form the 4-vectors k_i, k_f, p_f
c
c incident beam
      	    k_i(1)=0.
            k_i(2)=0.
            k_i(3)=e0
            k_i(4)=e0
c scattered electron
            k_f(1)=k2*sin(thee)
            k_f(2)=0.
            k_f(3)=k2*cos(thee)
            k_f(4)=k2
c outgoing proton
            p_f(1)=-ppf*cos(phi)*sin(theta)
            p_f(2)=ppf*sin(phi)
            p_f(3)=ppf*cos(phi)*cos(theta)
            p_f(4)=sqrt(ppf**2+mp**2)
c
            call subcc1(k_i,k_f,p_f,sigep,kfact,ppi,nucleon)
c
      sigep=sigep*1000.
c
      qlsq=4.*e0*(e0-omega)*sin(thee/2.)**2+omega**2
      er=sqrt(ppi**2+mn**2)
      ep=sqrt(ppf**2+mp**2)
      pipf=0.5*(qlsq-(ppf**2+ppi**2))

      rec=1.-ep/er*pipf/ppf**2

      er_h = sqrt(ppi**2+(15*mn)**2)
      rec_heavy=1.-ep/er_h*pipf/ppf**2

      f_rec = 1./rec

       
c get the paris momentum distribution
       
      rho=drho_paris( ppi ) 
c
      sigma=kfact*sigep*rho/rec
      sig_deutpwia = sigma
c
      thq = thq/dtr
      return

      end


