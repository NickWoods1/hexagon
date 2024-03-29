MODULE Armchair_Geometry
IMPLICIT NONE

CONTAINS

!=========================================================!
!=		     Armchair Subroutine		 =!
!=========================================================!

SUBROUTINE Armchair(wmax,lreal,lreal1,b,a,pos,nn1,nn2,nn3,L,lmax,Lwidth,a1,shear,ex,Poisson)
IMPLICIT NONE

INTEGER :: i,j, counter
DOUBLE PRECISION, INTENT(INOUT) :: shear,ex,Poisson
INTEGER, INTENT(INOUT) :: wmax,lreal,lreal1,b,a,L,lmax,Lwidth
INTEGER, DIMENSION(L,L), INTENT(INOUT) :: nn1,nn2,nn3
DOUBLE PRECISION, DIMENSION(2), INTENT(INOUT) :: a1
DOUBLE PRECISION, DIMENSION(L,2), INTENT(INOUT) :: pos


!=========================================================!
!=		Generating AC Coordinates		 =!
!=========================================================!

DO j=1,lmax


	b=(j-1)*Lwidth

	pos(1+b,1) = 0
	pos(1+b,2) = 0 - a*(j-1)*3

	counter = 2

	pos(counter+b,1) = pos(counter+b-1,1) - a*0.5*sqrt(3d0)
	pos(counter+b,2) = pos(counter+b-1,2) - a*0.5
	
	counter = counter + 1

	pos(counter+b,1) = pos(counter+b-1,1)
	pos(counter+b,2) = pos(counter+b-1,2) - 1*a
	
	counter = counter + 1

	pos(counter+b,1) = pos(counter+b-1,1) + a*0.5*sqrt(3d0)
	pos(counter+b,2) = pos(counter+b-1,2) - a*0.5

	counter = counter + 1

	pos(counter+b,1) = pos(counter+b-1,1) + a*0.5*sqrt(3d0)
	pos(counter+b,2) = pos(counter+b-1,2) + a*0.5

	counter = counter + 1

	pos(counter+b,1) = pos(counter+b-1,1)
	pos(counter+b,2) = pos(counter+b-1,2) + a*1


IF (wmax/=1)THEN

counter = counter + 1

	pos(counter+b,1) = pos(counter+b-1,1) + a*0.5*sqrt(3d0)
	pos(counter+b,2) = pos(counter+b-1,2) + a*0.5

	counter = counter + 1

	pos(counter+b,1) = pos(counter+b-3,1) + a*0.5*sqrt(3d0)
	pos(counter+b,2) = pos(counter+b-3,2) - a*0.5

IF (wmax/=2)THEN
DO i=1,(int(0.5*(wmax-1)))


	counter = counter + 1

	pos(counter+b,1) = pos(counter+b-2,1) + a*0.5*sqrt(3d0)
	pos(counter+b,2) = pos(counter+b-2,2) - a*0.5

	counter = counter + 1

	pos(counter+b,1) = pos(counter+b-2,1) + a*0.5*sqrt(3d0)
	pos(counter+b,2) = pos(counter+b-2,2) + a*0.5


	IF ((i==(int(0.5*(wmax-1)))).AND.(mod(wmax,2)/=0))THEN
	EXIT
	END IF

	counter = counter + 1

	pos(counter+b,1) = pos(counter+b-2,1) + a*0.5*sqrt(3d0)
	pos(counter+b,2) = pos(counter+b-2,2) + a*0.5

	counter = counter + 1

	pos(counter+b,1) = pos(counter+b-2,1) + a*0.5*sqrt(3d0)
	pos(counter+b,2) = pos(counter+b-2,2) - a*0.5

	IF ((i==(int(0.5*(wmax-1)))).AND.(mod(wmax,2)==0))THEN
	EXIT
	END IF


END DO
END IF
END IF
END DO


!=========================================================!
!=	   Creating the nearest neighbor lists   	 =!
!=========================================================!

DO i=(1+(lreal*Lwidth)),(L-(lreal*Lwidth))
	DO j=1,L
		IF (i/=j)THEN

			!=================================================================================================== First nn
			
			IF (sqrt((abs(pos(i,1)-pos(j,1)))**2 + (abs(pos(i,2)-pos(j,2)))**2)<(a+0.05))THEN 
			IF (j<(1+(lreal*Lwidth)))THEN
		
				nn1(i,j) = - j

			ELSEIF (j>(L-(lreal*Lwidth)))THEN

				nn1(i,j) = - (j - 2*Lwidth*lreal)

			ELSE

				nn1(i,j) = j - Lwidth*lreal

			END IF
			END IF

			!=================================================================================================== Second nn
		
			IF (sqrt((abs(pos(i,1)-pos(j,1)))**2 + (abs(pos(i,2)-pos(j,2)))**2)<(sqrt(3d0)*a+0.01))THEN
			IF (sqrt((abs(pos(i,1)-pos(j,1)))**2 + (abs(pos(i,2)-pos(j,2)))**2)>(a+0.05))THEN
			IF (j<(1+lreal*Lwidth))THEN

				nn2(i,j) = - j

			ELSEIF (j>(L-(lreal*Lwidth)))THEN

				nn2(i,j) = -(j - lreal*Lwidth*2)
	
			ELSE

				nn2(i,j) = j - lreal*Lwidth

			END IF
			END IF
			END IF	


			!================================================================================================= Third nn

			IF (sqrt((abs(pos(i,1)-pos(j,1)))**2 + (abs(pos(i,2)-pos(j,2)))**2)<(2*a+0.01))THEN 
			IF (sqrt((abs(pos(i,1)-pos(j,1)))**2 + (abs(pos(i,2)-pos(j,2)))**2)>(sqrt(3d0)*a+0.02))THEN 
			IF (j<(1+lreal*Lwidth))THEN

				nn3(i,j) = - j

			ELSEIF (j>(L-(lreal*Lwidth)))THEN

				nn3(i,j) = -(j - lreal*Lwidth*2)
	
			ELSE

				nn3(i,j) = j - lreal*Lwidth

			END IF
			END IF
			END IF


			!================================================================================================= 


		END IF
	END DO
END DO

!=========================================================!
!=	  	 End Nearest Neighbour list	    	 =!
!=========================================================!





!============= Smith algorithm for Shear Strain ==============!
DO i=1,L
	pos(i,2) = pos(i,2) + shear*pos(i,1)
END DO
!=============================================================!



!=========== Smith algorithm for Uniaxial Strain =============!
DO i=1,L
	pos(i,1) = pos(i,1)*(-Poisson*ex + 1)*a
	pos(i,2) = pos(i,2)*(ex + 1)*a 
END DO
!=============================================================!


!Write strained positions to positions.dat
DO i=1,L
	WRITE(2,*)pos(i,1),pos(i,2)
END DO




END SUBROUTINE

END MODULE Armchair_Geometry
