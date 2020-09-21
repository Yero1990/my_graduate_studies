      real*8 function drho_paris(pm)
c 
c return paris momentum distribution, does a polynomial interpolation
c between data points
c
c requires : polint and locate
c
      implicit none
      integer MAXDAT 
      parameter (MAXDAT=100)
      
      real*8 p(MAXDAT), rho(MAXDAT)
      real*8 pm, rho_, drho_

      integer j,k,m

      data m/4/  ! 4 point interpolation
      
c----------------------------------------------------------------------
c intialzie data: Paris momentum distribution
c
       data p(  1)/ 9.866E+00/,rho(  1)/ 1.540E+03/
       data p(  2)/ 1.973E+01/,rho(  2)/ 1.188E+03/
       data p(  3)/ 2.960E+01/,rho(  3)/ 8.177E+02/
       data p(  4)/ 3.947E+01/,rho(  4)/ 5.299E+02/
       data p(  5)/ 4.933E+01/,rho(  5)/ 3.358E+02/
       data p(  6)/ 5.920E+01/,rho(  6)/ 2.129E+02/
       data p(  7)/ 6.907E+01/,rho(  7)/ 1.366E+02/
       data p(  8)/ 7.893E+01/,rho(  8)/ 8.918E+01/
       data p(  9)/ 8.880E+01/,rho(  9)/ 5.933E+01/
       data p( 10)/ 9.867E+01/,rho( 10)/ 4.021E+01/
       data p( 11)/ 1.085E+02/,rho( 11)/ 2.773E+01/
       data p( 12)/ 1.184E+02/,rho( 12)/ 1.944E+01/
       data p( 13)/ 1.283E+02/,rho( 13)/ 1.383E+01/
       data p( 14)/ 1.381E+02/,rho( 14)/ 9.979E+00/
       data p( 15)/ 1.480E+02/,rho( 15)/ 7.289E+00/
       data p( 16)/ 1.579E+02/,rho( 16)/ 5.385E+00/
       data p( 17)/ 1.677E+02/,rho( 17)/ 4.021E+00/
       data p( 18)/ 1.776E+02/,rho( 18)/ 3.032E+00/
       data p( 19)/ 1.875E+02/,rho( 19)/ 2.308E+00/
       data p( 20)/ 1.973E+02/,rho( 20)/ 1.773E+00/
       data p( 21)/ 2.072E+02/,rho( 21)/ 1.374E+00/
       data p( 22)/ 2.171E+02/,rho( 22)/ 1.074E+00/
       data p( 23)/ 2.269E+02/,rho( 23)/ 8.461E-01/
       data p( 24)/ 2.368E+02/,rho( 24)/ 6.723E-01/
       data p( 25)/ 2.467E+02/,rho( 25)/ 5.387E-01/
       data p( 26)/ 2.565E+02/,rho( 26)/ 4.352E-01/
       data p( 27)/ 2.664E+02/,rho( 27)/ 3.546E-01/
       data p( 28)/ 2.763E+02/,rho( 28)/ 2.914E-01/
       data p( 29)/ 2.861E+02/,rho( 29)/ 2.416E-01/
       data p( 30)/ 2.960E+02/,rho( 30)/ 2.020E-01/
       data p( 31)/ 3.059E+02/,rho( 31)/ 1.703E-01/
       data p( 32)/ 3.157E+02/,rho( 32)/ 1.448E-01/
       data p( 33)/ 3.256E+02/,rho( 33)/ 1.242E-01/
       data p( 34)/ 3.355E+02/,rho( 34)/ 1.074E-01/
       data p( 35)/ 3.453E+02/,rho( 35)/ 9.357E-02/
       data p( 36)/ 3.552E+02/,rho( 36)/ 8.213E-02/
       data p( 37)/ 3.651E+02/,rho( 37)/ 7.259E-02/
       data p( 38)/ 3.749E+02/,rho( 38)/ 6.456E-02/
       data p( 39)/ 3.848E+02/,rho( 39)/ 5.777E-02/
       data p( 40)/ 3.947E+02/,rho( 40)/ 5.196E-02/
       data p( 41)/ 4.045E+02/,rho( 41)/ 4.697E-02/
       data p( 42)/ 4.144E+02/,rho( 42)/ 4.263E-02/
       data p( 43)/ 4.243E+02/,rho( 43)/ 3.884E-02/
       data p( 44)/ 4.341E+02/,rho( 44)/ 3.551E-02/
       data p( 45)/ 4.440E+02/,rho( 45)/ 3.255E-02/
       data p( 46)/ 4.539E+02/,rho( 46)/ 2.991E-02/
       data p( 47)/ 4.637E+02/,rho( 47)/ 2.754E-02/
       data p( 48)/ 4.736E+02/,rho( 48)/ 2.540E-02/
       data p( 49)/ 4.835E+02/,rho( 49)/ 2.346E-02/
       data p( 50)/ 4.933E+02/,rho( 50)/ 2.169E-02/
       data p( 51)/ 5.032E+02/,rho( 51)/ 2.007E-02/
       data p( 52)/ 5.131E+02/,rho( 52)/ 1.858E-02/
       data p( 53)/ 5.229E+02/,rho( 53)/ 1.722E-02/
       data p( 54)/ 5.328E+02/,rho( 54)/ 1.595E-02/
       data p( 55)/ 5.427E+02/,rho( 55)/ 1.478E-02/
       data p( 56)/ 5.525E+02/,rho( 56)/ 1.370E-02/
       data p( 57)/ 5.624E+02/,rho( 57)/ 1.270E-02/
       data p( 58)/ 5.723E+02/,rho( 58)/ 1.176E-02/
       data p( 59)/ 5.821E+02/,rho( 59)/ 1.089E-02/
       data p( 60)/ 5.920E+02/,rho( 60)/ 1.009E-02/
       data p( 61)/ 6.019E+02/,rho( 61)/ 9.333E-03/
       data p( 62)/ 6.117E+02/,rho( 62)/ 8.632E-03/
       data p( 63)/ 6.216E+02/,rho( 63)/ 7.979E-03/
       data p( 64)/ 6.315E+02/,rho( 64)/ 7.370E-03/
       data p( 65)/ 6.413E+02/,rho( 65)/ 6.804E-03/
       data p( 66)/ 6.512E+02/,rho( 66)/ 6.276E-03/
       data p( 67)/ 6.611E+02/,rho( 67)/ 5.785E-03/
       data p( 68)/ 6.709E+02/,rho( 68)/ 5.327E-03/
       data p( 69)/ 6.808E+02/,rho( 69)/ 4.902E-03/
       data p( 70)/ 6.907E+02/,rho( 70)/ 4.507E-03/
       data p( 71)/ 7.005E+02/,rho( 71)/ 4.140E-03/
       data p( 72)/ 7.104E+02/,rho( 72)/ 3.799E-03/
       data p( 73)/ 7.203E+02/,rho( 73)/ 3.483E-03/
       data p( 74)/ 7.301E+02/,rho( 74)/ 3.190E-03/
       data p( 75)/ 7.400E+02/,rho( 75)/ 2.918E-03/
       data p( 76)/ 7.499E+02/,rho( 76)/ 2.667E-03/
       data p( 77)/ 7.597E+02/,rho( 77)/ 2.434E-03/
       data p( 78)/ 7.696E+02/,rho( 78)/ 2.220E-03/
       data p( 79)/ 7.795E+02/,rho( 79)/ 2.021E-03/
       data p( 80)/ 7.893E+02/,rho( 80)/ 1.839E-03/
       data p( 81)/ 7.992E+02/,rho( 81)/ 1.670E-03/
       data p( 82)/ 8.091E+02/,rho( 82)/ 1.516E-03/
       data p( 83)/ 8.189E+02/,rho( 83)/ 1.373E-03/
       data p( 84)/ 8.288E+02/,rho( 84)/ 1.242E-03/
       data p( 85)/ 8.387E+02/,rho( 85)/ 1.123E-03/
       data p( 86)/ 8.485E+02/,rho( 86)/ 1.013E-03/
       data p( 87)/ 8.584E+02/,rho( 87)/ 9.123E-04/
       data p( 88)/ 8.683E+02/,rho( 88)/ 8.205E-04/
       data p( 89)/ 8.781E+02/,rho( 89)/ 7.367E-04/
       data p( 90)/ 8.880E+02/,rho( 90)/ 6.603E-04/
       data p( 91)/ 8.979E+02/,rho( 91)/ 5.908E-04/
       data p( 92)/ 9.077E+02/,rho( 92)/ 5.277E-04/
       data p( 93)/ 9.176E+02/,rho( 93)/ 4.704E-04/
       data p( 94)/ 9.275E+02/,rho( 94)/ 4.184E-04/
       data p( 95)/ 9.373E+02/,rho( 95)/ 3.715E-04/
       data p( 96)/ 9.472E+02/,rho( 96)/ 3.291E-04/
       data p( 97)/ 9.571E+02/,rho( 97)/ 2.908E-04/
       data p( 98)/ 9.669E+02/,rho( 98)/ 2.564E-04/
       data p( 99)/ 9.768E+02/,rho( 99)/ 2.255E-04/
       data p(100)/ 9.867E+02/,rho(100)/ 1.978E-04/
c----------------------------------------------------------------------
      
c from NR Vol2 p.113      

      call dlocate (p, MAXDAT, pm, j)
      k = min( max(j-(m-1)/2,1), MAXDAT+1-m)

c interpolate

      call dpolint(p(k), rho(k), m, pm, rho_, drho_)

      drho_paris = rho_*1.d-09     ! convert from GeV^-3 to MeV^-3

      return
      end




