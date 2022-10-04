!===========================================================================================!
!===========================================================================================!
! Numerical code In Fortran 90 Written by Ghasem Safaei                                     !
!                                                                                           !
! Created on 2015.04.25 with Version 1.0                                                    !
!                                                                                           !
!===========================================================================================!
!===========================================================================================!

Program orbit
Implicit none
Real x0,y0,z0,vx0,vy0,vz0,vr0,r0,h,x,xi,xf   
real pc,Mcore,G,select,yr,time
real a,b,q,r,v0,pi,vc,rho
integer ,parameter :: n=6
integer i,j,nstep
real y(n),dydx(n),yout(n)
CHARACTER(LEN=100) :: output
CHARACTER(LEN=100) :: aaa
CHARACTER(LEN=100) :: bbb
CHARACTER(LEN=100) :: ccc
common  /constant/	a,b,q,r,v0,pi,G,Mcore,vc,rho
external derivspoint
external derivsplumer
external derivsisochrone
external derivslogarithmic
external derivsmiyamoto
external derivskuzmin
external derivshernquist
external derivsjaffe
external derivsnfw
external derivssatoh

call selector(select,x0,y0,z0,r0,vx0,vy0,vz0,vr0,yr,time,xf,xi,h,nstep,x,y,a,b,q,r,rho,pi,G,Mcore)
call GM(G,Mcore)
call initial(x0,y0,z0,r0,vx0,vy0,vz0,vr0,yr,time,xf,xi,h,nstep,x,y)
!---------------------------------1----------------------------------------------------Point-Mass
if (select.eq.1) then
    open(1,file="Point-Mass.txt")
    write(1,*)y(1),y(3),y(5)
        do i=1,nstep
   	        call rk4(y,dydx,n,x,h,yout,derivspoint)
            x=x+h
           do j=1,n
	         call rk4(y,dydx,n,x,h,yout,derivspoint)
             y(j)=yout(j)
           end do
	      write(1,*)yout(1),yout(3),yout(5)	  
       end do

   !******
    output = 'Point-Mass'
    aaa='---'
    bbb='---'
    ccc='---'
    call script(output,aaa,bbb,ccc)
!--------------------------------------------2----------------------------------------- Plummer
else if (select.eq.2) then
   call constplumer(b)
   open(3,file="Plummer.txt")
   write(3,*)y(1),y(3),y(5)
  do i=1,nstep

   	 call rk4(y,dydx,n,x,h,yout,derivsplumer)
     x=x+h
       do j=1,n
	     call rk4(y,dydx,n,x,h,yout,derivsplumer)
         y(j)=yout(j)
       end do
	      write(3,*)yout(1),yout(3),yout(5)	  
  end do
!***
output = 'Plummer'
aaa='---'
bbb='b=0.05'
ccc='---'
call script(output,aaa,bbb,ccc)
!---------------------------------------------3----------------------------------------Isochrone
else if (select.eq.3) then
  call constisochrone(b)
  open(1,file="Isochrone.txt")
  write(1,*)y(1),y(3),y(5)
 do i=1,nstep
   	 call rk4(y,dydx,n,x,h,yout,derivsisochrone)
     x=x+h
       do j=1,n
	     call rk4(y,dydx,n,x,h,yout,derivsisochrone)
         y(j)=yout(j)
       end do
	      write(1,*)yout(1),yout(3),yout(5)	  
end do
output = 'Isochrone'
aaa='---'
bbb='b=0.05'
ccc='---'
call script(output,aaa,bbb,ccc)
!--------------------------------------4-----------------------------------------------Logarithmic
else if (select.eq.4) then
call constlogarithmic(v0,r,q)
  open(1,file="Logarithmic.txt")
  write(1,*)y(1),y(3),y(5)
 do i=1,nstep
   	 call rk4(y,dydx,n,x,h,yout,derivslogarithmic)
     x=x+h
       do j=1,n
	    call rk4(y,dydx,n,x,h,yout,derivslogarithmic)
        y(j)=yout(j)
       end do
	      write(1,*)yout(1),yout(3),yout(5)	  
 end do
output = 'Logarithmic'
aaa='-V_0=0.7351-'
bbb='-R_c=0.1-'
ccc='-q_phi=0.1-'
call script(output,aaa,bbb,ccc)
!--------------------------------------5----------------------------------------Miyamoto
else if (select.eq.5) then
call constmiyamoto(a,b)
  open(1,file="Miyamoto.txt")
  write(1,*)y(1),y(3),y(5)
 do i=1,nstep
   	 call rk4(y,dydx,n,x,h,yout,derivsmiyamoto)
     x=x+h
       do j=1,n
	   call rk4(y,dydx,n,x,h,yout,derivsmiyamoto)
         y(j)=yout(j)
         end do
	      write(1,*)yout(1),yout(3),yout(5)	  
 end do
output = 'Miyamoto'
aaa='-a=0.1-'
bbb='-b=0.1-'
ccc='----'
call script(output,aaa,bbb,ccc)
!------------------------------------6----------------------- Kuzmin
else if (select.eq.6) then
call constkuzmin(a)
  open(1,file="Kuzmin.txt")
  write(1,*)y(1),y(3),y(5)
 do i=1,nstep
   	 call rk4(y,dydx,n,x,h,yout,derivskuzmin)
     x=x+h
       do j=1,n
	     call rk4(y,dydx,n,x,h,yout,derivskuzmin)
         y(j)=yout(j)
       end do
	      write(1,*)yout(1),yout(3),yout(5)	  
 end do
output = 'Kuzmin'
aaa='-a=0.1-'
bbb='-----'
ccc='----'
call script(output,aaa,bbb,ccc)
!---------------------------------7------------------------------- Hernquist
else if (select.eq.7) then
call consthernquist(a)
   open(1,file="Hernquist.txt")
   write(1,*)y(1),y(3),y(5)
 do i=1,nstep
   	 call rk4(y,dydx,n,x,h,yout,derivshernquist)
     x=x+h
       do j=1,n
	     call rk4(y,dydx,n,x,h,yout,derivshernquist)
         y(j)=yout(j)
       end do
	      write(1,*)yout(1),yout(3),yout(5)	  
 end do
output = 'Hernquist'
aaa='-a=0.1-'
bbb='-----'
ccc='----'
call script(output,aaa,bbb,ccc)
!--------------------------------------------8--------------------------------------	NFW
else if (select.eq.8) then
  call constnfw(a)
  open(1,file="NFW.txt")
  write(1,*)y(1),y(3),y(5)
 do i=1,nstep
   	 call rk4(y,dydx,n,x,h,yout,derivsnfw)
     x=x+h
       do j=1,n
	     call rk4(y,dydx,n,x,h,yout,derivsnfw)
         y(j)=yout(j)
       end do
	      write(1,*)yout(1),yout(3),yout(5)	  
 end do
output = 'NFW'
aaa='-a=0.1-'
bbb='-----'
ccc='----'
call script(output,aaa,bbb,ccc)
!-----------------------------------------------9------------------------------------Jaffe
else if (select.eq.9) then
call constjaffe(a)
open(1,file="Jaffe.txt")
write(1,*)y(1),y(3),y(5)
do i=1,nstep
   	 call rk4(y,dydx,n,x,h,yout,derivsjaffe)
     x=x+h
       do j=1,n
	   call rk4(y,dydx,n,x,h,yout,derivsjaffe)
       y(j)=yout(j)
       end do
	      write(1,*)yout(1),yout(3),yout(5)	  
end do
output = 'Jaffe'
aaa='-a=0.1-'
bbb='-----'
ccc='----'
call script(output,aaa,bbb,ccc)
!----------------------------------10---------------------------------------Satoh
else if (select.eq.10) then
  call constsatoh(a,b)
  open(1,file="Satoh.txt")
  write(1,*)y(1),y(3),y(5)
 do i=1,nstep
   	 call rk4(y,dydx,n,x,h,yout,derivssatoh)
     x=x+h
       do j=1,n
	     call rk4(y,dydx,n,x,h,yout,derivssatoh)
         y(j)=yout(j)
       end do
	      write(1,*)yout(1),yout(3),yout(5)	  
 end do
output = 'Satoh'
aaa='-a=0.1-'
bbb='-b=0.1-'
ccc='----'
call script(output,aaa,bbb,ccc)
!else if (select.gt.10) then
!print*,'wrong input'
!stop
!else if (select.lt.0) then
!print*,'wrong input'

else
print*,'***wrong input---!-please enter a number-1-2-3-4-5-6-7-8-9-10-!----***'
stop


end if
close(1)
End Program orbit
!##################################################################################################################################################################
!***********************DERIVS***********************************
subroutine derivspoint(x,y,dydx)
real :: x,y(*),dydx(*),G,Mcore
common /constant/	a,b,q,r,v0,pi,G,Mcore,vc,rho
!***************x*************
  dydx(1)=y(2)
  dydx(2)=-G*Mcore*y(1)/((y(1)**2.+y(3)**2.+y(5)**2.)**1.5)
!***************y*************
  dydx(3)=y(4)
  dydx(4)=-G*Mcore*y(3)/((y(1)**2.+y(3)**2.+y(5)**2.)**1.5)
!***************z*************
  dydx(5)=y(6)
  dydx(6)=-G*Mcore*y(5)/((y(1)**2.+y(3)**2.+y(5)**2.)**1.5)
return
end
!********************
subroutine derivsplumer(x,y,dydx)
real :: x,y(*),dydx(*),G,Mcore,b
common /constant/	a,b,q,r,v0,pi,G,Mcore,vc,rho
!***************x*************
  dydx(1)=y(2)
  dydx(2)=-(G*Mcore*y(1))/((b**2.+y(1)**2.+y(3)**2.+y(5)**2.)**1.5)
!***************y*************
  dydx(3)=y(4)
  dydx(4)=-(G*Mcore*y(3))/((b**2.+y(1)**2.+y(3)**2.+y(5)**2.)**1.5)
!***************z*************
  dydx(5)=y(6)
  dydx(6)=-(G*Mcore*y(5))/((b**2.+y(1)**2.+y(3)**2.+y(5)**2.)**1.5)
return
end
!*****************
subroutine derivsisochrone(x,y,dydx)
real :: x,y(*),dydx(*),G,Mcore,b
common /constant/	a,b,q,r,v0,pi,G,Mcore,vc,rho
  !***************x*************
  dydx(1)=y(2)
  dydx(2)=(- G*Mcore*y(1)) / (((b**2.+ y(1)**2. + y(3)**2.  + y(5)**2. )**(1./2.))* ( b + (b**2.+ y(1)**2. + y(3)**2.+ y(5)**2.)**(1./2.))**2.)
  !***************y*************
  dydx(3)=y(4)
  dydx(4)=(- G*Mcore*y(3)) / (((b**2.+ y(1)**2. + y(3)**2.  + y(5)**2. )**(1./2.))* ( b + (b**2.+ y(1)**2. + y(3)**2.+ y(5)**2.)**(1./2.))**2.)
  !***************z*************
  dydx(5)=y(6)
  dydx(6)=(- G*Mcore*y(5)) / (((b**2.+ y(1)**2. + y(3)**2.  + y(5)**2. )**(1./2.))* ( b + (b**2.+ y(1)**2. + y(3)**2.+ y(5)**2.)**(1./2.))**2.) 
return
end
!********************
subroutine derivslogarithmic(x,y,dydx)
real :: x,y(*),dydx(*),G,Mcore,q,v0,r
common /constant/	a,b,q,r,v0,pi,G,Mcore,vc,rho
!***************x*************-(( x Subscript[v^2, 0])/(x^2+y^2+z^2/Subscript[q^2, \[Phi]]+Subscript[R^2, c]))
  dydx(1)=y(2)
  dydx(2)=-(v0**2.*y(1))/((y(1)**2.+y(3)**2.+(y(5)**2./q**2.)+r**2.))
!***************y*************
  dydx(3)=y(4)
  dydx(4)=-(v0**2.*y(3))/((y(1)**2.+y(3)**2.+(y(5)**2./q**2.)+r**2.))
!***************z*************
  dydx(5)=y(6)
  dydx(6)=-(v0**2.*y(5))/(((y(1)**2.+y(3)**2.+(y(5)**2./q**2.)+r**2.)))*q**2.
return
end
!***************
subroutine derivsmiyamoto(x,y,dydx)
real :: x,y(*),dydx(*),G,Mcore,b,a
common /constant/	a,b,q,r,v0,pi,G,Mcore,vc,rho
!***************x*************
  dydx(1)=y(2)
  dydx(2)=-(G*Mcore*y(1))/(  (y(1)**2.+y(3)**2.+(a+sqrt(b**2.+y(5)**2.))**2.)**1.5)
!***************y*************
  dydx(3)=y(4)
  dydx(4)=-(G*Mcore*y(3))/(  (y(1)**2.+y(3)**2.+(a+sqrt(b**2.+y(5)**2.))**2.)**1.5)
!***************z*************
  dydx(5)=y(6)
  dydx(6)=-(G*Mcore*y(5)*(a+sqrt(b**2.+y(5)**2.)))/(sqrt(b**2.+y(5)**2.)*(  (y(1)**2.+y(3)**2.+(a+sqrt(b**2.+y(5)**2.))**2.)**1.5)) 
return
end
!*****************************
subroutine derivskuzmin(x,y,dydx)
real :: x,y(*),dydx(*),G,Mcore,a
common /constant/	a,b,q,r,v0,pi,G,Mcore,vc,rho
!***************x*************-((G M x)/(x^2+y^2+(a+z)^2)^(3/2))
  dydx(1)=y(2)
  dydx(2)=-((G*Mcore*y(1))/(y(1)**2.+y(3)**2.+(a+y(5))**2.)**1.5)
!***************y*************
  dydx(3)=y(4)
  dydx(4)=-((G*Mcore*y(3))/(y(1)**2.+y(3)**2.+(a+y(5))**2.)**1.5)
!***************z*************
  dydx(5)=y(6)
  dydx(6)=-((G*Mcore*y(5))/(y(1)**2.+y(3)**2.+(a+y(5))**2.)**1.5) 
return
end
!*********************************************hernquist
subroutine derivshernquist(x,y,dydx)
real :: x,y(*),dydx(*),G,Mcore,pi,a,rho
common /constant/	a,b,q,r,v0,pi,G,Mcore,vc,rho
rho=(3.*Mcore)/(4.*pi*a**3.)
 !***************x*************-((2 a G \[Pi]  x[t] Subscript[\[Rho], 0])/(Sqrt[x[t]^2+y[t]^2+z[t]^2] (1+Sqrt[x[t]^2+y[t]^2+z[t]^2]/a)^2))
  dydx(1)=y(2)
  dydx(2)=-2.*a*G*pi*y(1)*rho/(sqrt(y(1)**2+y(3)**2+y(5)**2)*(1.+ sqrt(y(1)**2+y(3)**2+y(5)**2) / a )**2)
!***************y*************
  dydx(3)=y(4)
  dydx(4)=-2.*a*G*pi*y(3)*rho/(sqrt(y(1)**2+y(3)**2+y(5)**2)*(1.+ sqrt(y(1)**2+y(3)**2+y(5)**2) / a )**2)
!***************z*************
  dydx(5)=y(6)
  dydx(6)=-2.*a*G*pi*y(5)*rho/(sqrt(y(1)**2+y(3)**2+y(5)**2)*(1.+ sqrt(y(1)**2+y(3)**2+y(5)**2) / a )**2)
return
end
!***********************************jaffe
subroutine derivsjaffe(x,y,dydx)
real :: x,y(*),dydx(*),G,Mcore,pi,a,rho
common /constant/	a,b,q,r,v0,pi,G,Mcore,vc,rho
rho=(3.*Mcore)/(4.*pi*a**3.)
!***************x************* -((4 a^3 G \[Pi] x Subscript[\[Rho], 0])/((x^2+y^2+z^2)^(3/2) (1. +a/Sqrt[x^2+y^2+z^2])))
  dydx(1)=y(2)
  dydx(2)=-(4.*a**3.*G*pi*rho*y(1))/((y(1)**2.+y(3)**2.+y(5)**2.)**1.5 *(1.+a/(sqrt(y(1)**2.+y(3)**2.+y(5)**2.))))
!***************y*************
  dydx(3)=y(4)
  dydx(4)=-(4.*a**3.*G*pi*rho*y(3))/((y(1)**2.+y(3)**2.+y(5)**2.)**1.5 *(1.+a/(sqrt(y(1)**2.+y(3)**2.+y(5)**2.))))
!***************z*************
  dydx(5)=y(6)
  dydx(6)=-(4.*a**3.*G*pi*rho*y(5))/((y(1)**2+y(3)**2.+y(5)**2.)**1.5 *(1.+a/(sqrt(y(1)**2.+y(3)**2.+y(5)**2.))))
return
end
!********************************************************NFW
subroutine derivsnfw(x,y,dydx)
real :: x,y(*),dydx(*),G,Mcore,a,pi,rho
common /constant/	a,b,q,r,v0,pi,G,Mcore,vc,rho
rho=(3.*Mcore)/(4.*pi*a**3.)
!***************x*************(2 z+(2 a z)/Sqrt[b^2+z^2])
  dydx(1)=y(2)
  dydx(2)=(4.*a**2.*G*pi*y(1)*rho)/((y(1)**2.+y(3)**2.+y(5)**2.)*(1.+sqrt(y(1)**2.+y(3)**2.+y(5)**2.)/a))-& 
  &(4.*a**3.*G*pi*y(1)*rho*log((1.+sqrt(y(1)**2.+y(3)**2.+y(5)**2.)/a))/((y(1)**2.+y(3)**2.+y(5)**2.)**1.5))
!***************y*************
  dydx(3)=y(4)
  dydx(4)=(4.*a**2.*G*pi*y(3)*rho)/((y(1)**2.+y(3)**2.+y(5)**2.)*(1.+sqrt(y(1)**2.+y(3)**2.+y(5)**2.)/a))-& 
  &(4.*a**3.*G*pi*y(3)*rho*log((1.+sqrt(y(1)**2.+y(3)**2.+y(5)**2.)/a))/((y(1)**2.+y(3)**2.+y(5)**2.)**1.5))
!***************z*************
  dydx(5)=y(6)
  dydx(6)=(4.*a**2.*G*pi*y(5)*rho)/((y(1)**2.+y(3)**2.+y(5)**2.)*(1.+sqrt(y(1)**2.+y(3)**2.+y(5)**2.)/a))-& 
  &(4.*a**3.*G*pi*y(5)*rho*log((1.+sqrt(y(1)**2.+y(3)**2.+y(5)**2.)/a))/((y(1)**2.+y(3)**2.+y(5)**2.)**1.5))
return
end
!!!****************************************************************satoh
subroutine derivssatoh(x,y,dydx)
real :: x,y(*),dydx(*),G,Mcore,a,b
common /constant/	a,b,q,r,v0,pi,G,Mcore,vc,rho
!***************x*************(2 z+(2 a z)/Sqrt[b^2+z^2])
  dydx(1)=y(2)
  dydx(2)=-(G*Mcore*y(1))/(y(1)**2.+y(3)**2.+y(5)**2.+a*(a+2.*sqrt(b**2.+y(5)**2.)))**1.5
!***************y*************
  dydx(3)=y(4)
  dydx(4)=-(G*Mcore*y(3))/(y(1)**2.+y(3)**2.+y(5)**2.+a*(a+2.*sqrt(b**2.+y(5)**2.)))**1.5
!***************z*************
  dydx(5)=y(6)
  dydx(6)= -(G*Mcore*(y(5)+a*y(5)/sqrt(b**2.+y(5)**2.)) )/(y(1)**2.+y(3)**2.+y(5)**2.+a*(a+2.*sqrt(b**2.+y(5)**2.)))**1.5
return
end
!*******************************rk4**************************************
SUBROUTINE rk4(y,dydx,n,x,h,yout,derivs)
INTEGER n,NMAX
REAL h,x,dydx(n),y(n),yout(n)
EXTERNAL derivs
PARAMETER (NMAX=50)
INTEGER i
REAL h6,hh,xh,dym(NMAX),dyt(NMAX),yt(NMAX)
hh=h*0.5
h6=h/6.
xh=x+hh
do 11 i=1,n
yt(i)=y(i)+hh*dydx(i)
11 continue
call derivs(xh,yt,dyt)
do 12 i=1,n
yt(i)=y(i)+hh*dyt(i)
12 continue
call derivs(xh,yt,dym)
do 13 i=1,n
yt(i)=y(i)+h*dym(i)
dym(i)=dyt(i)+dym(i)
13 continue
call derivs(x+h,yt,dyt)
do 14 i=1,n
yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i))
	14 continue
return
END
!***************************************************************
 subroutine script(output,aaa,bbb,ccc)
CHARACTER(LEN=100) :: output
CHARACTER(LEN=100) :: aaa
CHARACTER(LEN=100) :: bbb
CHARACTER(LEN=100) :: ccc
!_________
open (4,file=TRIM(ADJUSTL(output))//'-script.txt')
write(4,*)"set term png"
write(4,*)'set output "', TRIM(ADJUSTL(output)), '.png"'
write(4,*)"set pm3d"
write(4,*)"set isosample 50"
write(4,*)"set samples 1000"
write(4,*)"set multiplot"
write(4,*)"set nokey"
write(4,*)'set xlabel "X"'
write(4,*)'set ylabel "Y"'
write(4,*)'set zlabel "Z"'
write(4,*)'set title   "',TRIM(ADJUSTL(output)),'--', TRIM(ADJUSTL(aaa)),'--',TRIM(ADJUSTL(bbb)),'-t=10000.0-',TRIM(ADJUSTL(ccc)),'--V_0=0.7351--"'
write(4,*)'set label "Galaxy"'
write(4,*)"set grid"
write(4,*)"set style data lines"
write(4,*)"sp 'C:\Users\ghasem\Desktop\exercise\",TRIM(ADJUSTL(output)),".txt' w l"
write(4,*)"quit"
!------
return
end
!*******************************************************************************
subroutine initial(x0,y0,z0,r0,vx0,vy0,vz0,vr0,yr,time,xf,xi,h,nstep,x,y)
real x0,y0,z0,r0,vx0,vy0,vz0,vr0,yr,time,xf,xi,h,x,y(*)
integer nstep
print*,"------------------------initial position------------------------------"
print*,"----------------------X_0------------------------------"
print*,'please enter X_0 as initial value for position align X axis '
read*,x0
print*,"----------------------Y_0------------------------------"
print*,'please enter Y_0 as initial value for position align Y axis '
read*,y0
print*,"----------------------Z_0------------------------------"
print*,'please enter Z_0 as initial value for position align Z axis '
read*,z0
r0=sqrt(x0**2+y0**2+z0**2)
print*,"------------------------initial velocity------------------------------"
print*,"----------------------Vx_0------------------------------"
print*,'please enter Vx_0 as initial value for velocity align X axis '
read*,vx0
print*,"----------------------Vy_0------------------------------"
print*,'please enter Vy_0 as initial value for velocity align Y axis '
read*,vy0
print*,"----------------------Vz_0------------------------------"
print*,'please enter Vz_0 as initial value for velocity align Z axis '
read*,vz0
vr0=sqrt(vx0**2+vy0**2+vz0**2)
print*,"----------------------time------------------------------"
print*,'please enter time(according to yr)'
read*,time
yr=365.25*24.*60.*60.
time=time*yr
xf=time/1.e7
print*,"----------------------step------------------------------"
print*,'please enter step(according to yr)'
read *,h

xi=0.
nstep=(xf-xi)/h
x=xi

y(1)=x0
y(2)=vx0
y(3)=y0
y(4)=vy0
y(5)=z0
y(6)=vz0

return
end
!**********************************************************
subroutine GM(G,Mcore)
real G,Mcore
print*,'please enter "G" gravitation constant'
read*,G
print*,'please enter "M" Mass of inside'
read*,Mcore
return
end
!***************
subroutine constplumer(b)
real b
print*,'please enter b for Plummer model'
read*,b
return
end
!****************************
subroutine constisochrone(b)
real b
print*,'please enter b for Isochrone model'
read*,b
return
end
!****************************
subroutine constmiyamoto(a,b)
real a,b
print*,'please enter "a" for Miyamoto model'
read*,a
print*,'please enter "b" for Miyamoto model'
read*,b
return
end
!****************************
subroutine constlogarithmic(v0,r,q)
real v0,r,q
print*,'please enter "v0" for Logarithmic model'
read*,v0
print*,'please enter "r" for Logarithmic model'
read*,r
print*,'please enter "q" for Logarithmic model'
read*,q
return
end
!****************************
subroutine constjaffe(a)
real a
print*,'please enter "a" for Jaffe model'
read*,a
return
end
!****************************
subroutine consthernquist(a)
real a
print*,'please enter "a" for Hernquist model'
read*,a
return
end
!****************************
subroutine constnfw(a)
real a
print*,'please enter "a" for NFW model'
read*,a
return
end
!****************************
subroutine constsatoh(a,b)
real a,b
print*,'please enter "a" for Satoh model'
read*,a
print*,'please enter "b" for Satoh model'
read*,b
return
end
!****************************
subroutine constkuzmin(a)
real a
print*,'please enter "a" for Kuzmin model'
read*,a
return
end
!****************************
subroutine selector(select,x0,y0,z0,r0,vx0,vy0,vz0,vr0,yr,time,xf,xi,h,nstep,x,y,a,b,q,r,rho,pi,G,Mcore)
real G,Mcore,pi,rho,a,b,q,r
real x0,y0,z0,r0,vx0,vy0,vz0,vr0,yr,time,xf,xi,h,x,y(*),pc
integer nstep

!*******************------initial-value----****************************

 print*,"please enter a number betwen (1-10) this number determine potential type"
print*,"------------------------------------------------------"
print*,"1=Point mass************** (for Point mass potential )"
print*,"------------------------------------------------------"
print*,"2=Plummer******************** (for Plummer potential )"
print*,"------------------------------------------------------"
print*,"3=Isochrone**************** (for Isochrone potential )"
print*,"------------------------------------------------------"
print*,"4=Logarithmic************ (for Logarithmic potential )"
print*,"------------------------------------------------------"
print*,"5=Miyamoto****************** (for Miyamoto potential )"
print*,"------------------------------------------------------"
print*,"6=Kuzmin********************** (for Kuzmin potential )"
print*,"------------------------------------------------------"
print*,"7=Hernquist**************** (for Hernquist potential )"
print*,"------------------------------------------------------"
print*,"8=NFW ****************************(for NFW potential )"
print*,"------------------------------------------------------"
print*,"9=Jaffe************************ (for Jaffe potential )"
print*,"------------------------------------------------------"
print*,"10=Satoh ***********************(for Satoh potential )"
read*,select

!********
a=0.1
b=0.1
q=0.1
r=0.1
pi=3.141592
pc=3.085
x0=8.0*pc 
y0=0.
z0=0.
G=6.67
Mcore=2.
rho=(3.*Mcore)/(4.*pi*a**3.)
r0=sqrt(x0**2+y0**2+z0**2)
vc=sqrt(G*Mcore/r0)
v0=vc
vx0=0.
vy0=-vc
vz0=0.
vr0=sqrt(vx0**2+vy0**2+vz0**2)
yr=365.25*24.*60.*60.
time=time*yr
xf=time/1.e7
xf=1000.
xi=0.
h=0.01
nstep=(xf-xi)/h
time=sqrt(1.e31)*xf/yr
x=xi
y(1)=x0
y(2)=vx0
y(3)=y0
y(4)=vy0
y(5)=z0
y(6)=vz0
!************

return
end
!****************************
!****************************
!****************************
!****************************
!****************************
!****************************
!****************************
