!! HERE THE MC-VOLUME-MOVES !!!

CALL energysystem()
CALL lbox()
CALL checkindip(good)

! you save all old quatities:
eneold=global%energy
vold=global%volume
DO kd=1,global%dimension
   raold(kd)=global%ra(kd)
   rbold(kd)=global%rb(kd)       
   rcold(kd)=global%rc(kd)
ENDDO
DO ip=1,global%nr_particles
 DO kd=1,global%dimension
    cmold(ip,kd)=lattice%cm(ip,kd)
  DO ipv=1,2
     pvold(ip,kd,ipv)=lattice%patch_vector(ip,kd,ipv)
  ENDDO
 ENDDO
ENDDO

! you to to fractional coordinates
DO ip=1,global%nr_particles
 CALL realtofractional(ip,lattice%cm,lattice%fcm)
ENDDO


!you change the volume
ran=9*random()
IF(ran.LT.3)THEN
   IF(ran.LT.1)i=1
   IF(ran.GE.1.AND.ran.LT.2)i=2
   IF(ran.GE.2.AND.ran.LT.3)i=3
   global%ra(i)=global%ra(i)-global%deltav*(random()-0.5_dp)
ENDIF
IF(ran.GE.3.AND.ran.LT.6)THEN
  IF(ran.LT.4)i=1
  IF(ran.GE.4.AND.ran.LT.5)i=2
  IF(ran.GE.5.AND.ran.LT.6)i=3
  global%rb(i)=global%rb(i)-global%deltav*(random()-0.5_dp)
ENDIF
IF(ran.GE.6)THEN
  IF(ran.LT.7)i=1
  IF(ran.GE.7.AND.ran.LT.8)i=2
  IF(ran.GE.8.AND.ran.LT.9)i=3
  global%rc(i)=global%rc(i)-global%deltav*(random()-0.5_dp)
ENDIF

! here it comes the trick, you calculate the new box
CALL latticereduction3D()
CALL lbox()
CALL checkindip(good)
vnew=global%volume

!check if the box is too small before checking all overlaps (cutoff of the interaction)
CALL  shortestlength(shortest)
IF(shortest.LT.global%diameter*1.1_dp)THEN
  global%energy=eneold
  DO kd=1,global%dimension
     global%ra(kd)=raold(kd)
     global%rb(kd)=rbold(kd)
     global%rc(kd)=rcold(kd)
  ENDDO
  CALL lbox()
  CALL checkindip(good)
  RETURN
ENDIF
     
!you go back to real space
DO ip=1,global%nr_particles
   CALL fractionaltoreal(ip,lattice%fcm,lattice%cm)
ENDDO

!you calculate your overlaps
!keep in mind: calculate nvals=nr of boxes you need to consider
CALL shortestlength(mindist)
nvals=INT(global%range/mindist)+1
DO jpj=1,global%nr_particles
 DO ipi=1,jpj
    rijbare(1)=lattice%cm(ipi,1)-lattice%cm(jpj,1)
    rijbare(2)=lattice%cm(ipi,2)-lattice%cm(jpj,2)
    rijbare(3)=lattice%cm(ipi,3)-lattice%cm(jpj,3)
    DO k1=-nvals,nvals
     DO k2=-nvals,nvals
      DO k3=-nvals,nvals
        IF (ipi.NE.jpj.OR.k1.NE.0.OR.k2.NE.0.OR.k3.NE.0) THEN 
        rij(1)=rijbare(1)+k1*global%ra(1)+k2*global%rb(1)+k3*global%rc(1)
        rij(2)=rijbare(2)+k1*global%ra(2)+k2*global%rb(2)+k3*global%rc(2)
        rij(3)=rijbare(3)+k1*global%ra(3)+k2*global%rb(3)+k3*global%rc(3)  
        dist2=rij(1)**2+rij(2)**2+rij(3)**2
         IF(dist2.LT.sigma2)THEN
           overlap=.true.
           DO ip=1,global%nr_particles
            DO kd=1,global%dimension
               lattice%cm(ip,kd)=cmold(ip,kd)
             DO ipv=1,2
               lattice%patch_vector(ip,kd,ipv)=pvold(ip,kd,ipv)
             ENDDO
            ENDDO
            CALL patchposition(ip)
           ENDDO
           global%energy=eneold
           DO kd=1,global%dimension
              global%ra(kd)=raold(kd)
              global%rb(kd)=rbold(kd)
              global%rc(kd)=rcold(kd)
           ENDDO
           CALL lbox()
           CALL checkindip(good)
           RETURN
         ENDIF
        ENDIF
      ENDDO
     ENDDO
    ENDDO
 ENDDO
ENDDO


!!! HERE THE ROUTINES !!!

SUBROUTINE lbox()
INTEGER (kind=di) :: i

global%boxa=0._dp
global%boxb=0._dp
global%boxc=0._dp

DO i=1,global%dimension
   global%boxa=global%boxa+global%ra(i)**2._dp
   global%boxb=global%boxb+global%rb(i)**2._dp
   global%boxc=global%boxc+global%rc(i)**2._dp
ENDDO
global%boxa=SQRT(global%boxa)
global%boxb=SQRT(global%boxb)
global%boxc=SQRT(global%boxc)

END SUBROUTINE lbox

SUBROUTINE checkindip(answer)
INTEGER (kind=di) :: i
LOGICAL answer

answer=.true.

global%detbox=global%ra(1)*global%rb(2)*global%rc(3)-global%ra(1)*global%rb(3)*global%rc(2)+&
global%ra(3)*global%rb(1)*global%rc(2)-global%ra(2)*global%rb(1)*global%rc(3)+&
global%ra(2)*global%rb(3)*global%rc(1)-global%ra(3)*global%rb(2)*global%rc(1)

IF(ABS(global%detbox).LT.0.0000000000001_dp)THEN
answer=.false.
WRITE(*,*)'Box Vectors not Independent!'
STOP
ENDIF

global%volume=abs(global%detbox)

END SUBROUTINE checkindip

SUBROUTINE shortestlength(short)
INTEGER (kind=di) :: i,kd
REAL (kind=dp) :: short,short1,short2,short3,short4
REAL (kind=dp),DIMENSION(:) :: la(global%dimension),lb(global%dimension),lc(global%dimension)
REAL (kind=dp) :: lam,lbm,lcm
REAL (kind=dp),DIMENSION(:) :: shortest(4)

CALL lbox()
DO kd=1,global%dimension
   la(kd)=global%ra(kd)
   lb(kd)=global%rb(kd)
   lc(kd)=global%rc(kd)
ENDDO
lam=global%boxa
lbm=global%boxb
lcm=global%boxc

CALL shortestl2d(la,lb,lam,lbm,short1)
CALL shortestl2d(la,lc,lam,lcm,short2)
CALL shortestl2d(lb,lc,lbm,lcm,short3)

CALL shortestl3d(la,lb,lc,short4)

shortest(1)=short1
shortest(2)=short2
shortest(3)=short3
shortest(4)=short4

short=shortest(1)
DO i=2,4
   IF(shortest(i).LT.short) short=shortest(i)
ENDDO       

END SUBROUTINE shortestlength

SUBROUTINE shortestl2d(r1,r2,r1l,r2l,shortestl)
REAL (kind=dp), INTENT(IN) :: r1(global%dimension)
REAL (kind=dp), INTENT(IN) :: r2(global%dimension)
REAL (kind=dp), INTENT(IN):: r1l
REAL (kind=dp), INTENT(IN):: r2l
REAL (kind=dp), INTENT(OUT):: shortestl
REAL (kind=dp) :: theta,height1,height2
REAL (kind=dp) :: dotproduct12

dotproduct12=r1(1)*r2(1)+r1(2)*r2(2)+r1(3)*r2(3)

theta=ACOS(dotproduct12/(r1l*r2l))
height1=r1l*SIN(theta)
height2=r2l*SIN(theta)

shortestl=height1
IF(height2.LT.shortestl) shortestl= height2

END SUBROUTINE shortestl2d

SUBROUTINE shortestl3d(r1,r2,r3,shortestl)
INTEGER (kind=di) :: i
REAL (kind=dp), INTENT(IN) :: r1(global%dimension)
REAL (kind=dp), INTENT(IN) :: r2(global%dimension)
REAL (kind=dp), INTENT(IN) :: r3(global%dimension)
REAL (kind=dp), INTENT(OUT):: shortestl
REAL (kind=dp) :: height1,height2,height3,d1
REAL (kind=dp), DIMENSION(:) :: normal1(global%dimension),normal2(global%dimension),normal3(global%dimension)
REAL (kind=dp), DIMENSION(:) :: crossproduct(global%dimension)

crossproduct(1)=r1(2)*r2(3)-r1(3)*r2(2)
crossproduct(2)=-1._dp*(r1(1)*r2(3)-r1(3)*r2(1))
crossproduct(3)=r1(1)*r2(2)-r1(2)*r2(1)
DO i=1,global%dimension
   normal1(i)=crossproduct(i)
ENDDO

crossproduct(1)=r1(2)*r3(3)-r1(3)*r3(2)
crossproduct(2)=-1._dp*(r1(1)*r3(3)-r1(3)*r3(1))
crossproduct(3)=r1(1)*r3(2)-r1(2)*r3(1)
DO i=1,global%dimension
   normal2(i)=crossproduct(i)
ENDDO

crossproduct(1)=r3(2)*r2(3)-r3(3)*r2(2)
crossproduct(2)=-1._dp*(r3(1)*r2(3)-r3(3)*r2(1))
crossproduct(3)=r3(1)*r2(2)-r3(2)*r2(1) 
DO i=1,global%dimension
   normal3(i)=crossproduct(i)
ENDDO

height1 = ABS(distancefromplane(normal1,r1,r3))
height2 = ABS(distancefromplane(normal2,r3,r2))
height3 = ABS(distancefromplane(normal3,r2,r1))

shortestl = height1
IF (height2.LT.shortestl) shortestl = height2
IF (height3.LT.shortestl) shortestl = height3

END SUBROUTINE shortestl3d

FUNCTION distancefromplane(nn,r0,rr)
REAL (kind=dp) :: distancefromplane
INTEGER (kind=di) :: i
REAL (kind=dp), INTENT(IN) :: rr(global%dimension)
REAL (kind=dp), INTENT(IN) :: r0(global%dimension)
REAL (kind=dp), INTENT(IN) :: nn(global%dimension)
REAL (kind=dp), DIMENSION(:) :: bb(global%dimension)
REAL (kind=dp) :: dotproductbn,dotproductnn

DO i=1,3
   bb(i)=rr(i)-r0(i)
ENDDO

dotproductbn=bb(1)*nn(1)+bb(2)*nn(2)+bb(3)*nn(3)
dotproductnn=nn(1)*nn(1)+nn(2)*nn(2)+nn(3)*nn(3) 

distancefromplane= dotproductbn/SQRT(dotproductnn)

END FUNCTION distancefromplane

SUBROUTINE latticereduction3D
INTEGER (kind=di) :: counter,i
REAL (kind=dp) :: rasafe(global%dimension),rbsafe(global%dimension),rcsafe(global%dimension)
REAL (kind=dp) :: raold(global%dimension),rbold(global%dimension),rcold(global%dimension)
REAL (kind=dp) :: rtest(global%dimension)
REAL (kind=dp) :: surface,surfaceold,surfacenew
LOGICAL good

DO i=1,global%dimension
   rasafe(i)=global%ra(i)
   rbsafe(i)=global%rb(i) 
   rcsafe(i)=global%rc(i)
ENDDO
surface=surface3d(rasafe,rbsafe,rcsafe)

DO i=1,global%dimension
   raold(i)=global%ra(i)
   rbold(i)=global%rb(i) 
   rcold(i)=global%rc(i)
   rasafe(i)=global%ra(i)
   rbsafe(i)=global%rb(i) 
   rcsafe(i)=global%rc(i)
ENDDO
   surfaceold=surface
        
counter=0
DO
  counter=counter+1
  DO i=1,global%dimension
     rtest(i)=0._dp
     rtest(i)=global%ra(i)+global%rb(i)
  ENDDO
  surfacenew = surface3D(rtest,rbsafe,rcsafe)
  IF(surfacenew.LT.surfaceold)THEN
     surfaceold = surfacenew
     DO i=1,global%dimension
        raold(i)=rtest(i)
        rbold(i)=global%rb(i) 
        rcold(i)=global%rc(i)
     ENDDO
  ENDIF
  DO i=1,global%dimension
     rtest(i)=0._dp
     rtest(i)=global%ra(i)-global%rb(i)
  ENDDO
  surfacenew = surface3D(rtest,rbsafe,rcsafe)
  IF(surfacenew.LT.surfaceold)THEN
     surfaceold = surfacenew
     DO i=1,global%dimension
        raold(i)=rtest(i)
        rbold(i)=global%rb(i) 
        rcold(i)=global%rc(i)
     ENDDO
  ENDIF
  DO i=1,global%dimension
     rtest(i)=0._dp
      rtest(i)=global%ra(i)+global%rc(i)
  ENDDO
  surfacenew = surface3D(rtest,rbsafe,rcsafe)
  IF(surfacenew.LT.surfaceold)THEN
        surfaceold = surfacenew
        DO i=1,global%dimension
           raold(i)=rtest(i)
           rbold(i)=global%rb(i) 
           rcold(i)=global%rc(i)
        ENDDO
  ENDIF
  DO i=1,global%dimension
     rtest(i)=0._dp
     rtest(i)=global%ra(i)-global%rc(i)
  ENDDO
  surfacenew = surface3D(rtest,rbsafe,rcsafe)
  IF(surfacenew.LT.surfaceold)THEN
     surfaceold = surfacenew
     DO i=1,global%dimension
        raold(i)=rtest(i)
        rbold(i)=global%rb(i) 
        rcold(i)=global%rc(i)
     ENDDO
  ENDIF
  DO i=1,global%dimension
     rtest(i)=0._dp
     rtest(i)=global%rb(i)+global%rc(i)
  ENDDO
  surfacenew = surface3D(rasafe,rtest,rcsafe)
  IF(surfacenew.LT.surfaceold)THEN
     surfaceold = surfacenew
     DO i=1,global%dimension
        raold(i)=global%ra(i)
        rbold(i)=rtest(i) 
        rcold(i)=global%rc(i)
     ENDDO
  ENDIF
  DO i=1,global%dimension
     rtest(i)=0._dp
     rtest(i)=global%rb(i)-global%rc(i)
  ENDDO
  surfacenew = surface3D(rasafe,rtest,rcsafe)
  IF(surfacenew.LT.surfaceold)THEN
     surfaceold = surfacenew
     DO i=1,global%dimension
        raold(i)=global%ra(i)
        rbold(i)=rtest(i) 
        rcold(i)=global%rc(i)
     ENDDO
  ENDIF
  DO i=1,global%dimension
     rtest(i)=0._dp
     rtest(i)=global%rb(i)+global%ra(i)
  ENDDO
  surfacenew = surface3D(rasafe,rtest,rcsafe)
  IF(surfacenew.LT.surfaceold)THEN
     surfaceold = surfacenew
     DO i=1,global%dimension
        raold(i)=global%ra(i)
        rbold(i)=rtest(i) 
        rcold(i)=global%rc(i)
     ENDDO
  ENDIF
  DO i=1,global%dimension
      rtest(i)=0._dp
      rtest(i)=global%rb(i)-global%ra(i)
  ENDDO
  surfacenew = surface3D(rasafe,rtest,rcsafe)
  IF(surfacenew.LT.surfaceold)THEN
     surfaceold = surfacenew
     DO i=1,global%dimension
        raold(i)=global%ra(i)
        rbold(i)=rtest(i) 
        rcold(i)=global%rc(i)
     ENDDO
  ENDIF
  DO i=1,global%dimension
     rtest(i)=0._dp
     rtest(i)=global%rc(i)+global%ra(i)
  ENDDO
  surfacenew = surface3D(rasafe,rbsafe,rtest)
  IF(surfacenew.LT.surfaceold)THEN
     surfaceold = surfacenew
     DO i=1,global%dimension
        raold(i)=global%ra(i)
        rbold(i)=global%rb(i) 
        rcold(i)=rtest(i)
     ENDDO
  ENDIF
  DO i=1,global%dimension
     rtest(i)=0._dp
     rtest(i)=global%rc(i)-global%ra(i)
  ENDDO
  surfacenew = surface3D(rasafe,rbsafe,rtest)
  IF(surfacenew.LT.surfaceold)THEN
     surfaceold = surfacenew
     DO i=1,global%dimension
        raold(i)=global%ra(i)
        rbold(i)=global%rb(i) 
        rcold(i)=rtest(i)
     ENDDO
  ENDIF
  DO i=1,global%dimension
     rtest(i)=0._dp
     rtest(i)=global%rc(i)+global%rb(i)
  ENDDO
  surfacenew = surface3D(rasafe,rbsafe,rtest)
  IF(surfacenew.LT.surfaceold)THEN
     surfaceold = surfacenew
     DO i=1,global%dimension
        raold(i)=global%ra(i)
        rbold(i)=global%rb(i) 
        rcold(i)=rtest(i)
     ENDDO
  ENDIF
  DO i=1,global%dimension
     rtest(i)=0._dp
     rtest(i)=global%rc(i)-global%rb(i)
  ENDDO
  surfacenew = surface3D(rasafe,rbsafe,rtest)
  IF(surfacenew.LT.surfaceold)THEN
     surfaceold = surfacenew
     DO i=1,global%dimension
        raold(i)=global%ra(i)
        rbold(i)=global%rb(i) 
        rcold(i)=rtest(i)
     ENDDO
  ENDIF

  IF(surface-surfaceold.LT.0.0000001_dp)THEN
     EXIT
  ELSE
     surface = surfaceold
     DO i=1,global%dimension
        global%ra(i) = raold(i)
        global%rb(i) = rbold(i)
        global%rc(i) = rcold(i)
     ENDDO
  ENDIF
  
  IF(counter.EQ.100)EXIT
ENDDO

CALL checkindip(good)

IF(.NOT.good)THEN
   print*,'Back to old vectors in Lattice Reduction ',good
   DO i=1,global%dimension
      global%ra(i)=rasafe(i)
      global%rb(i)=rbsafe(i)
      global%rc(i)=rcsafe(i)
    ENDDO
ENDIF

END SUBROUTINE latticereduction3D


FUNCTION surface3d(raa,rbb,rcc)
REAL (kind=dp) :: surface3d
REAL (kind=dp) :: raa(global%dimension),rbb(global%dimension),rcc(global%dimension)

surface3d=2._dp*((raa(2)*rbb(3)-raa(3)*rbb(2))*(raa(2)*rbb(3)-raa(3)*rbb(2))+(raa(3)*rbb(1)-raa(1)*rbb(3))*(raa(3)*rbb(1)-raa(1)*rbb(3))+(raa(1)*rbb(2)-raa(2)*rbb(1))*(raa(1)*rbb(2)-raa(2)*rbb(1))+(raa(2)*rcc(3)-raa(3)*rcc(2))*(raa(2)*rcc(3)-raa(3)*rcc(2))+(raa(3)*rcc(1)-raa(1)*rcc(3))*(raa(3)*rcc(1)-raa(1)*rcc(3))+(raa(1)*rcc(2)-raa(2)*rcc(1))*(raa(1)*rcc(2)-raa(2)*rcc(1))+(rbb(2)*rcc(3)-rbb(3)*rcc(2))*(rbb(2)*rcc(3)-rbb(3)*rcc(2))+(rbb(3)*rcc(1)-rbb(1)*rcc(3))*(rbb(3)*rcc(1)-rbb(1)*rcc(3))+(rbb(1)*rcc(2)-rbb(2)*rcc(1))*(rbb(1)*rcc(2)-rbb(2)*rcc(1)))

END FUNCTION surface3d
