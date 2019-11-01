MODULE Zigzag_Geometry
IMPLICIT NONE

CONTAINS



!=========================================================!
!=		     ZigZag Subroutine			 =!
!=========================================================!


SUBROUTINE Zigzag(wmax,lreal,lreal1,b,a,pos,nn1,nn2,nn3,L,lmax,Lwidth,a1,shear,ex,Poisson)
IMPLICIT NONE

INTEGER :: i,j, counter
DOUBLE PRECISION, INTENT(INOUT) :: shear,ex,Poisson
INTEGER, INTENT(INOUT) :: wmax,lreal,lreal1,b,a,L,lmax,Lwidth
INTEGER, DIMENSION(L,L), INTENT(INOUT) :: nn1,nn2,nn3
DOUBLE PRECISION, DIMENSION(L,2), INTENT(INOUT) :: pos
DOUBLE PRECISION, DIMENSION(2), INTENT(INOUT) :: a1



!=========================================================!
!=		Generating ZZ Coordinates		 =!
!=========================================================!


DO j=1,lmax

	b=(j-1)*Lwidth

	pos(1+b,1) = 0
	pos(1+b,2) = 0 - a*(j-1)*(sqrt(3d0))

	counter = 2

DO i=1,int((0.5*wmax))+1

	pos(counter+b,1) = (pos(counter+b-1,1) + 0.5)*a
	pos(counter+b,2) = (pos(counter+b-1,2) - 0.5*sqrt(3d0))*a

		IF ((i==int((0.5*wmax))+1).AND.(mod(wmax,2)==0))THEN !i=wmax and wmax is even
			EXIT
		END IF

	counter = counter + 1

	pos(counter+b,1) = (pos(counter+b-1,1) + 1)*a
	pos(counter+b,2) = (pos(counter+b-1,2) + 0 )*a

	counter = counter + 1

	pos(counter+b,1) = (pos(counter+b-1,1) + 0.5)*a
	pos(counter+b,2) = (pos(counter+b-1,2) + 0.5*sqrt(3d0) )*a

		IF ((i==int((0.5*wmax))+1).AND.(mod(wmax,2)/=0))THEN !i=wmax and wmax is odd
			EXIT
		END IF

	counter = counter + 1

	pos(counter+b,1) = (pos(counter+b-1,1) + 1)*a
	pos(counter+b,2) = (pos(counter+b-1,2) + 0)*a

	counter = counter + 1

END DO
END DO

!=========================================================!
!=	   Creating the nearest neighbor lists   	 =!
!=========================================================!



DO i=(1+(2+(wmax*2))*lreal1),(L-(2+(wmax*2))*lreal1)
	DO j=1,L
		IF (i/=j)THEN

			!=================================================================================================== First nn

			IF (sqrt((abs(pos(i,1)-pos(j,1)))**2 + (abs(pos(i,2)-pos(j,2)))**2)<(a+0.1))THEN 
	
			IF (j<=(2+(wmax*2))*lreal1)THEN
		
				IF (lreal==1)THEN
					nn1(i,j) = - (j - Lwidth)
				
				ELSEIF (lreal/=1)THEN
					nn1(i,j) = - j
				END IF

			ELSEIF (j>(L-(2+(wmax*2))*lreal1))THEN


				IF (lreal==1)THEN
					nn1(i,j) = -(j - 3*Lwidth)
				ELSEIF (lreal/=1)THEN
					nn1(i,j) = -(j - lreal*Lwidth*2)
				END IF
	
			ELSE

			IF (lreal==1)THEN
			nn1(i,j) = j - 2*Lwidth
			ELSE
			nn1(i,j) = j - lreal*Lwidth
			END IF



			END IF	
			END IF


			!=================================================================================================== Second nn

			IF (sqrt((abs(pos(i,1)-pos(j,1)))**2 + (abs(pos(i,2)-pos(j,2)))**2)<(sqrt(3d0)*a+0.01))THEN
			IF (sqrt((abs(pos(i,1)-pos(j,1)))**2 + (abs(pos(i,2)-pos(j,2)))**2)>(a+0.1))THEN
			IF (j<=(2+(wmax*2))*lreal1)THEN

				IF (lreal==1)THEN
					nn2(i,j) = - (j - Lwidth)
				ELSEIF (lreal/=1)THEN
					nn2(i,j) = - j
				END IF

			ELSEIF (j>(L-(2+(wmax*2))*lreal1))THEN

				IF (lreal==1)THEN
					nn2(i,j) = -(j - 3*Lwidth)
				ELSEIF (lreal/=1)THEN
					nn2(i,j) = -(j - 2*lreal*Lwidth)
				END IF
	
			ELSE

			IF (lreal==1)THEN
			nn2(i,j) = j - 2*Lwidth
			ELSE
			nn2(i,j) = j - lreal*Lwidth
			END IF


			END IF	
			END IF
			END IF

			!================================================================================================= Third nn

			IF (sqrt((abs(pos(i,1)-pos(j,1)))**2 + (abs(pos(i,2)-pos(j,2)))**2)<(2*a+0.01))THEN 
			IF (sqrt((abs(pos(i,1)-pos(j,1)))**2 + (abs(pos(i,2)-pos(j,2)))**2)>(sqrt(3d0)*a+0.01))THEN 

			IF (j<=(2+(wmax*2))*lreal1)THEN

				IF (lreal==1)THEN
					nn3(i,j) = - (j - Lwidth)
				ELSEIF (lreal/=1)THEN
					nn3(i,j) = - j
				END IF

			ELSEIF (j>(L-(2+(wmax*2))*lreal1))THEN

				IF (lreal==1)THEN
					nn3(i,j) = -(j - 3*Lwidth)
				ELSEIF (lreal/=1)THEN
					nn3(i,j) = -(j - 2*lreal*Lwidth)
				END IF
	
			ELSE
			
			IF (lreal==1)THEN
			nn3(i,j) = j - 2*Lwidth
			ELSE
			nn3(i,j) = j - lreal*Lwidth
			END IF


			END IF	
			END IF
			END IF

			!================================================================================================= 

		END IF
	END DO
END DO


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


END MODULE Zigzag_Geometry
