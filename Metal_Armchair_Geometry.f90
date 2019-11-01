MODULE Metal_Armchair_Geometry
IMPLICIT NONE

CONTAINS

!=========================================================!
!=	      Metal Doped Armchair Subroutine     	 =!
!=========================================================!


SUBROUTINE Metal(wmax,lreal,lreal1,b,a,pos,nn1,nn2,nn3,L,lmax,Lwidth,a1,N)
IMPLICIT NONE

INTEGER :: i,j, counter
INTEGER, INTENT(INOUT) :: wmax,lreal,lreal1,b,a,L,lmax,Lwidth,N
INTEGER, DIMENSION(L,L), INTENT(INOUT) :: nn1,nn2,nn3
DOUBLE PRECISION, DIMENSION(2), INTENT(INOUT) :: a1
DOUBLE PRECISION, DIMENSION(L,2), INTENT(INOUT) :: pos


!=========================================================!
!=		Generating MAC Coordinates		 =!
!=========================================================!

DO j=1,lmax

	b=(j-1)*N

	pos(1+b,1) = 0
	pos(1+b,2) = 0 - (j-1)*4

	counter = 2

	pos(counter+b,1) = pos(counter+b-1,1) - 0.5*sqrt(3d0)
	pos(counter+b,2) = pos(counter+b-1,2) - 0.5
	
	counter = counter + 1

	pos(counter+b,1) = pos(counter+b-1,1)
	pos(counter+b,2) = pos(counter+b-1,2) - 1
	
	counter = counter + 1

	pos(counter+b,1) = pos(counter+b-1,1) + 0.5*sqrt(3d0)
	pos(counter+b,2) = pos(counter+b-1,2) - 0.5

	counter = counter + 1

	pos(counter+b,1) = pos(counter+b-1,1) + 0.5*sqrt(3d0)
	pos(counter+b,2) = pos(counter+b-1,2) + 0.5

	counter = counter + 1

	pos(counter+b,1) = pos(counter+b-1,1)
	pos(counter+b,2) = pos(counter+b-1,2) + 1

	counter = counter + 1

	pos(counter+b,1) = pos(counter+b-6,1) 
	pos(counter+b,2) = pos(counter+b-6,2) + 1

IF (wmax/=1)THEN

	counter = counter + 1

	pos(counter+b,1) = pos(counter+b-2,1) + 0.5*sqrt(3d0)
	pos(counter+b,2) = pos(counter+b-2,2) + 0.5

	counter = counter + 1

	pos(counter+b,1) = pos(counter+b-4,1) + 0.5*sqrt(3d0)
	pos(counter+b,2) = pos(counter+b-4,2) - 0.5

	counter = counter + 1

	pos(counter+b,1) = pos(counter+b-3,1) + sqrt(3d0)
	pos(counter+b,2) = pos(counter+b-3,2)

IF (wmax/=2)THEN
DO i=1,(int(0.5*(wmax-1)))


	counter = counter + 1

	pos(counter+b,1) = pos(counter+b-3,1) + 0.5*sqrt(3d0)
	pos(counter+b,2) = pos(counter+b-3,2) - 0.5

	counter = counter + 1

	pos(counter+b,1) = pos(counter+b-3,1) + 0.5*sqrt(3d0)
	pos(counter+b,2) = pos(counter+b-3,2) + 0.5


	IF ((i==(int(0.5*(wmax-1)))).AND.(mod(wmax,2)/=0))THEN
	EXIT
	END IF

	counter = counter + 1

	pos(counter+b,1) = pos(counter+b-2,1) + 0.5*sqrt(3d0)
	pos(counter+b,2) = pos(counter+b-2,2) + 0.5

	counter = counter + 1

	pos(counter+b,1) = pos(counter+b-2,1) + 0.5*sqrt(3d0)
	pos(counter+b,2) = pos(counter+b-2,2) - 0.5

	counter = counter + 1

	pos(counter+b,1) = pos(counter+b-5,1) + sqrt(3d0)
	pos(counter+b,2) = pos(counter+b-5,2)

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

DO i=1+N,2*N
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


DO i=(1+(lreal*Lwidth)),(L-(lreal*Lwidth))
	DO j=1,L

	if (nn1(i,j)/=0)THEN
	!!PRINT*, nn1(i,j)
	end if

	end do

!PRINT*, "=="
end do



!=========================================================!
!=	   Creating DFT Cartesian Coordinate set   	 =!
!=========================================================!


DO j=1,lmax

	b=(j-1)*N

	pos(1+b,1) = 0
	pos(1+b,2) = 0 - (j-1)*4.601335

	counter = 2

	pos(counter+b,1) = pos(counter+b-1,1) - 0.836328
	pos(counter+b,2) = pos(counter+b-1,2) - 0.4916316
	
	counter = counter + 1

	pos(counter+b,1) = pos(counter+b-1,1)
	pos(counter+b,2) = pos(counter+b-1,2) - 0.9849997
	
	counter = counter + 1

	pos(counter+b,1) = pos(counter+b-1,1) + 0.836328
	pos(counter+b,2) = pos(counter+b-1,2) - 0.4916316

	counter = counter + 1

	pos(counter+b,1) = pos(counter+b-1,1) + 0.87129789
	pos(counter+b,2) = pos(counter+b-1,2) + 0.4699066

	counter = counter + 1

	pos(counter+b,1) = pos(counter+b-1,1)
	pos(counter+b,2) = pos(counter+b-1,2) + 1.0299747

	counter = counter + 1

	pos(counter+b,1) = pos(counter+b-6,1) 
	pos(counter+b,2) = pos(counter+b-6,2) + 1.3155552

IF (wmax/=1)THEN

	counter = counter + 1

	pos(counter+b,1) = pos(counter+b-2,1) + 0.87129789
	pos(counter+b,2) = pos(counter+b-2,2) + 0.4699341

	counter = counter + 1

	pos(counter+b,1) = pos(counter+b-4,1) + 0.8711303
	pos(counter+b,2) = pos(counter+b-4,2) - 0.469612

	counter = counter + 1

	pos(counter+b,1) = pos(counter+b-3,1) + 1.692307
	pos(counter+b,2) = pos(counter+b-3,2)

IF (wmax/=2)THEN
DO i=1,(int(0.5*(wmax-1)))


	counter = counter + 1

	pos(counter+b,1) = pos(counter+b-3,1) + 0.8356697
	pos(counter+b,2) = pos(counter+b-3,2) - 0.4928334

	counter = counter + 1

	pos(counter+b,1) = pos(counter+b-3,1) + 0.8356697
	pos(counter+b,2) = pos(counter+b-3,2) + 0.4928334


	IF ((i==(int(0.5*(wmax-1)))).AND.(mod(wmax,2)/=0))THEN
	EXIT
	END IF

	counter = counter + 1

	pos(counter+b,1) = pos(counter+b-2,1) + 0.87129789
	pos(counter+b,2) = pos(counter+b-2,2) + 0.4699341

	counter = counter + 1

	pos(counter+b,1) = pos(counter+b-2,1) + 0.8711303
	pos(counter+b,2) = pos(counter+b-2,2) - 0.469612

	counter = counter + 1

	pos(counter+b,1) = pos(counter+b-5,1) + 1.692307
	pos(counter+b,2) = pos(counter+b-5,2)

	IF ((i==(int(0.5*(wmax-1)))).AND.(mod(wmax,2)==0))THEN
	EXIT
	END IF


END DO
END IF
END IF
END DO


DO i=1,L

WRITE(2,*)pos(i,1),pos(i,2)

END DO




END SUBROUTINE








END MODULE
