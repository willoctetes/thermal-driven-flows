program TPD_tube_cl
!
!    *******************************************************************************************************
!    *******************************************************************************************************
!    **                                                                                                   **
!    **  Program by           :    Giorgos Tatsios														  **
!    **  First release        :    20 February 2020				                                          **
!	 **  Last modified		  :	   06 April 2020	(Guillermo López Quesada)							  **
!    **  Problem              :    Single Gas Flow through a circular tube                                **
!    **  Collision model      :    Shakhov model                                                          **
!    **  Components           :    Single gas                                                             **
!    **                                                                                                   **
!    *******************************************************************************************************
!    *******************************************************************************************************
!-----------------------------------------------------------------------------------------------------------
!
!Declaration of the variables to be used in the code
!
implicit none
integer::i,n,j,k,valsd,valsa,casesol,iter
real(8)::dz,error, pi,delin,delout,mu_ref,Pin,rgas, len, g, gav, dav, omega, tc, th, t_ref, a, convergencecrit, Tin, Tout, Rin, Rout
real(8), allocatable::delta(:),deldat(:),adat(:),qpdat(:,:), qtdat(:,:), P(:), T(:), u0(:), mu(:)
Real(8), allocatable:: spdat(:),stdat(:),asp(:),ast(:), rad(:), z(:)
Real(8):: gain, m_init, m,pout
Integer:: valsasp,valsast,status
namelist /physical_parameters/ mu_ref, t_ref, pin, pout, rgas, omega, a
namelist /numerical_parameters/ gain, m_init,casesol, convergencecrit

! This code is used for calculating 3 different cases:
	!casesol = 0 -> Pin and Pout are known, m is calculated (maximum m)
	!casesol = 1 -> Pin is known and m=0, Pout is calculated (maximum TPD)
	!casesol = 2 -> Pin and m are known, Pout is calculated (intermediate case)



pi=dacos(-1.0d0)

!OPEN and read the physical and numerical parameters
open(11,file='parameters')
read(11,NML=physical_parameters)
read(11,NML=numerical_parameters)
close(11)

!ASSOCIATE the geometry file to the code
 open(11,file='geometry')				!Open the geometry file
 read(11,*) n,len						!read and associate the number of nodes of the channel in the code
 allocate(delta(n),P(n), T(n), u0(n), mu(n), rad(n),z(n))
 do i = 1,n
	 read(11,*) z(i),rad(i),t(i)			!read and associate the length, radius and temperature of the channel for each node
 enddo
 close(11)

 do i = 1,n	
	 u0(i)=sqrt(2.0d0*rgas*t(i))				!calculate the most probable velocity in each node
	 mu(i) = mu_ref*(t(i)/t_ref)**omega		!calculate the viscosity in each node
 enddo

! NEXT LINES FOR SIMPLE AND SINGLE GEOMETRY CREATION WITHOUT USING THE SCRIPT

! casesol=0
! m_init=1e-17

! n=10001				!number of nodes
! len=5e-5		!length of the channel
! dz=len/(n-1.0d0)	!differential increasae for the nodes in the length direction
! Tin=300				!inlet temperature
! Tout=400			!outlet temperature
! Rin=5e-7			!inlet radius
! Rout=5e-7			!outlet radius

!allocate(delta(n),P(n), T(n), u0(n), mu(n), rad(n),z(n))	!allocation of the vector size to the number of nodes

! do i=1,n
	! z(i)=0+dz*(i-1)								!length direction
	! T(i)=Tin + (i-1)*dz*(Tout-Tin)/len			!linear distribution of temperature (it can be changed to other T distributions)
	! rad(i)=Rin + (i-1)*dz*(Rout-Rin)/len		!radius of each node (normally straight circular channel)
! enddo

!write(*,*) 'GEOMETRY completed'

! do i = 1,n	
	! u0(i)=sqrt(2.0d0*rgas*t(i))				!calculate the most probable velocity in each node
	! mu(i) = mu_ref*(t(i)/t_ref)**omega		!calculate the viscosity in each node
! enddo
	

delin=pin*rad(1)/(sqrt(2*rgas*t(1))*mu_ref*(t(1)/t_ref)**omega)			!delta at inlet
delout=pout*rad(n)/(sqrt(2*rgas*t(n))*mu_ref*(t(n)/t_ref)**omega)		!delta at outlet

!number of different values of At and An at the Sp database
valsASp = 10

!number of different values of At and An at the St database
valsASt = 5

!number of different values of delta, At and An at the qp and qt database
valsd = 17
valsa = 11

allocate(deldat(valsd),adat(valsa),qpdat(valsd,valsa), qtdat(valsd,valsa))
allocate(spdat(valsASp),stdat(valsASt),asp(valsASp),ast(valsASt))


!read database for qp and qt
open(11,file='qpqtdata')
do j = 1,valsa
	read(11,*) adat(j)
	do i = 1,valsd
		read(11,*) deldat(i),qpdat(i,j),qtdat(i,j)
	enddo
enddo
close(11)

!read database for sp
open(11,file='spdata')
do i = 1,valsASp
	read(11,*) asp(i),spdat(i)
enddo
close(11)

!read database for st
open(11,file='stdata')
do i = 1,valsASt
	read(11,*) ast(i),stdat(i)
enddo
close(11)

m=m_init									!m_init is the initial mass flow rate to iterate for casesol=0 or the given m for casesol=2
if(casesol == 1) m=0d0						!casesol=1, m=0 and Pin are given to calculate Pout

P(1)=Pin									!Pin is the inlet PRESSURE (always known for the 3 cases)
delta(1)=Pin*rad(1)/(mu(1)*u0(1))
iter = 0
do 
	iter = iter + 1
	!solve the equation to find the pressure at the other end
	do i = 2,n
		P(i) = P(i-1) - (z(i)-z(i-1))*u0(i-1)/pi/rad(i)**3/qp(delta(i-1),a)*m + (qt(delta(i-1),a)/qp(delta(i-1),a))*(P(i-1)/T(i-1))*(t(i)-t(i-1))
		delta(i) = P(i)*rad(i)/(mu(i)*u0(i))
	enddo
	
	if(casesol==1 .or. casesol==2) exit								!exit the loop, no iteration required.
	if(mod(iter,1)==0) write(*,*) iter, m, p(n), pout, p(n)-pout
	
	if(abs(p(n)-pout)<convergencecrit) exit							!iterate upon convergence
	m=m*(1d0+gain*(p(n)-pout)/p(1))									!recalculation of m for the next iteration
enddo

!WRITE of the output data
open(11,file='presdist.dat')
write(11,*) 'VARIABLES = z, r, rtop, rbot, t, P, delta'
write(11,*) 'ZONE T="1", I=',n/100 + 1
do i = 1,n,100
	write(11,*) z(i),rad(i),rad(i)/2d0,-rad(i)/2d0,t(i),p(i),delta(i)
enddo
close(11)

!WRITE the output depending on the case selected
if(casesol == 0) then
	write(*,*) iter, p(n), p(n)-p(1),m
else if (casesol == 1) then
	write(*,*) iter, p(n), p(n)-p(1),m, log(p(n)/p(1))/log(t(n)/t(1))
else if (casesol == 2) then
	write(*,*) iter, p(n), p(n)-p(1),m
endif

contains


!interpolation for sp
function sp(aa)
real(8)::sp,aa
integer:: i1, i2

if(aa>aSp(valsaSp) .or. aa<aSp(1)) then
	write(*,*) 'a out of bound for sp',aa,asp(1),asp(valsaSp)
	stop
else
	do i2 = 2,valsaSp
		if(asp(i2)>aa .or. iseq(asp(i2),aa)) exit
	enddo
	i1 = i2 - 1
	
	sp = (spdat(i1)*(asp(i2)-aa) + spdat(i2)*(aa-asp(i1)))/(asp(i2)-asp(i1))
endif

end function sp

!interpolation for st
function st(aa)
real(8)::st,aa
integer:: i1, i2

if(aa>aSt(valsaSt) .or. aa<aSt(1)) then
	write(*,*) 'a out of bound for st',aa,ast(1),ast(valsaSt)
	stop
else
	do i2 = 2,valsaSt
		if(ast(i2)>aa .or. iseq(ast(i2),aa)) exit
	enddo
	i1 = i2 - 1
	
	st = (stdat(i1)*(ast(i2)-aa) + stdat(i2)*(aa-ast(i1)))/(ast(i2)-ast(i1))
endif

end function st

!trilinear interpolation for qp
function qp(del,aa)
real(8)::qp,del,aa
integer::i1d,i2d,i1a,i2a

if(del>deldat(valsd)) then
    qp = del/4.0d0+sp(aa)

else if(del<deldat(1)) then
    write(*,*) 'delta smaller than',deldat(1),del
    stop
else if (aa> adat(valsa) .or. aa<adat(1)) then
	write(*,*) 'a out of range',adat(1),adat(valsa),aa
	stop
else

	do i2d=1,valsd
		if(deldat(i2d)>=del) exit
	enddo
	i1d=i2d-1

	do i2a=1,valsa
		if(adat(i2a)>=aa) exit
	enddo
	i1a=i2a-1

	qp = (&
	(&
	(qpdat(i1d,i1a)*(adat(i2a)-aa) + qpdat(i1d,i2a)*(aa-adat(i1a)))/(adat(i2a)-adat(i1a))&
	 )*(deldat(i2d)-del) +&
	(&
	(qpdat(i2d,i1a)*(adat(i2a)-aa) + qpdat(i2d,i2a)*(aa-adat(i1a)))/(adat(i2a)-adat(i1a))&
	)*(del-deldat(i1d)) &
	 )/(deldat(i2d)-deldat(i1d))
endif
end function qp

!trilinear interpolation for qt
function qt(del,aa)
real(8)::qt,del,aa
integer::i1d,i2d,i1a,i2a

if(del>deldat(valsd)) then
    qt = st(aa)/del

else if(del<deldat(1)) then
    write(*,*) 'delta smaller than',deldat(1),del
    stop
else if (aa> adat(valsa) .or. aa<adat(1)) then
	write(*,*) 'a out of range',adat(1),adat(valsa),aa
	stop
else

	do i2d=1,valsd
		if(deldat(i2d)>=del) exit
	enddo
	i1d=i2d-1

	do i2a=1,valsa
		if(adat(i2a)>=aa) exit
	enddo
	i1a=i2a-1

	qt = (&
	(&
	(qtdat(i1d,i1a)*(adat(i2a)-aa) + qtdat(i1d,i2a)*(aa-adat(i1a)))/(adat(i2a)-adat(i1a))&
	 )*(deldat(i2d)-del) +&
	(&
	(qtdat(i2d,i1a)*(adat(i2a)-aa) + qtdat(i2d,i2a)*(aa-adat(i1a)))/(adat(i2a)-adat(i1a))&
	)*(del-deldat(i1d)) &
	 )/(deldat(i2d)-deldat(i1d))
endif
end function qt

!checks if two float numbers are equal to some tolerance 1d-8
function iseq(x,y)
logical:: iseq
real(8)::x,y
	iseq=.FALSE.
	if(abs(x-y)<1d-8) then
		iseq=.TRUE.
	endif
end function iseq

end program TPD_tube_cl
