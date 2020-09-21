      subroutine relkin(omega,qlab,cth, 
     &                 conv,cthnpc,enp,qc,ppp,ppi,lerr)
      implicit double precision (a-h,o-z)
c
c 
c calculates the relevant kinematic variables in the center of 
c mass system and the final and initial proton momenta in the 
c laboratory frame. 
c 
c the input varable cth = cos (theta) where theta is the angle 
c between q and the out coming proton. 
c cthnpc is the cosine of the theta-np angle in the cms frame 
c 
      real*8 ml,kl2sq 
      logical lerr
      
      include 'constants.h'

      dtr=pi/180. 
      rtd=180./pi 
c 
c calcualate final momentum 
c 
      ql=qlab
      e0=md+omega 
      a0=ql**2-e0**2-mp**2+mn**2 
      a=4*e0**2-4*ql**2*cth**2 
      b=4*ql*a0*cth 
      c=4*e0**2*mp**2-a0**2 
      discr=b**2-4*a*c
      if (discr.lt.0) then
         lerr=.true.
         return
      endif
      ppp=(-b+sqrt(b**2-4*a*c))/2/a 
c other solution
      ppp2=(-b-sqrt(b**2-4*a*c))/2/a 
c      print *, ' sol. 1 :', ppp, ' sol2 : ', ppp2
c 
c calculate theta-np lab and the initial momentum 
c 
      fak1=-0.5*(ql-2*ppp*cth)
      fak2=sqrt(ppp**2+(ql/2)**2-ppp*ql*cth)
      cthnpl=fak1/fak2
      kl2sq=ppp**2+(ql/2)**2-ppp*ql*cth
      ppi=sqrt(kl2sq+(ql/2)**2-sqrt(kl2sq)*ql*cthnpl)
c
c make a Lorentz tranformation to the CMS frame
c
      beta=ql/(md+omega)
      gamma= 1/sqrt(1-beta**2)
c
      qc=gamma*(ql-beta*omega)
      wc=gamma*(omega-beta*ql)
c
      cthnpc=md*cthnpl/sqrt(md**2+qc**2-qc**2*cthnpl**2)
c 
c calculate invariant mass
c
      et=sqrt((omega+md)**2-ql**2)
      enp=et-(mn+mp)
c
c     cross section transformation factor cms->lab
c
      pcm=sqrt(enp*(enp+2.*(mn+mp)))
      ppcm= pcm/2.
      ml=enp+mn+mp
      el=sqrt(ml**2+ql**2)
      fact2=(ppp/ppcm)**3*ml/el
      conv=fact2/(1.+ql*ml*cthnpc/(2*ppcm*el))
c
      return
      end


