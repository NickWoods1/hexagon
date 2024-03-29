MODULE Construct_AC_Hamiltonian
IMPLICIT NONE

CONTAINS




!=========================================================!
!=	         Hamiltonian Subroutine	AC		 =!
!=========================================================!

SUBROUTINE HamiltonianAC(wmax,lreal,lreal1,pos,nn1,nn2,nn3,t1,t2,t3,N,counter,Im,k,a,Lwidth,L,Ham,rij,E,nup,Uhub,EP)
IMPLICIT NONE

INTEGER :: i,j,EP
DOUBLE PRECISION :: strain,aedge,b,t1_edge
INTEGER, INTENT(INOUT) :: wmax,lreal,lreal1,Lwidth,a,L,N,counter
INTEGER, DIMENSION(L,L), INTENT(INOUT) :: nn1,nn2,nn3
DOUBLE PRECISION, DIMENSION(L,2), INTENT(INOUT) :: pos
COMPLEX*16, DIMENSION(L), INTENT(INOUT) :: nup
COMPLEX*16, DIMENSION(N,N), INTENT(INOUT) :: Ham,t1,t2,t3
DOUBLE PRECISION, DIMENSION(N), INTENT(INOUT) :: E,Uhub
COMPLEX*16, INTENT(INOUT) :: Im,rij,k


!=========================================================!
!=		    Hamiltonian Loop AC			 =!
!=========================================================!

	!==============================================================================================================
	!=====================================================================================First Nearest Neighbor==
	!==============================================================================================================

	DO j=(1+(lreal*Lwidth)),(L-(lreal*Lwidth)) !loop over atoms per unit

		DO i=1,L !constructing the hamiltonian


		!================================================================================================
		!============================================================================Edge Pertubations==
		!================================================================================================
		
		IF (EP==1)THEN

		t1_edge = 3.4

		t1(2,3) = t1_edge
		t1(3,2) = t1_edge

		t1(N,N-1) = t1_edge
		t1(N-1,N) = t1_edge

		END IF

		!===================================================================================================
		!===========================================================================End Edge Pertubations==
		!===================================================================================================

		
		IF ((i/=j).AND.(nn1(j,i)/=0))THEN


		!========================================================================================
		!===============================================================================Strain==
		!========================================================================================

		rij = sqrt((pos(i,1)-pos(j,1))**2 + (pos(i,2)-pos(j,2))**2) !atomic displacement 

		b = 1d0 !beta parameter

		!strain = (1/(real(rij)))**2 !Harrison
		
		strain = exp( -b*(real(rij)-1) )  !Exponential		

		!========================================================================================
		!===========================================================================End Strain==
		!========================================================================================


		IF (nn1(j,i)<0)THEN !Intra-cell hopping
			IF ((pos(j,2)-pos(i,2))<0)THEN
				Ham(counter,-nn1(j,i))=Ham(counter,-nn1(j,i))-strain*t1(counter,-nn1(j,i))*exp(Im*k*a)
			ELSE
				Ham(counter,-nn1(j,i))=Ham(counter,-nn1(j,i))-strain*t1(counter,-nn1(j,i))*exp(-Im*k*a)
			END IF

		ELSEIF (nn1(j,i)>0)THEN
			Ham(counter,nn1(j,i))=Ham(counter,nn1(j,i))-strain*t1(counter,nn1(j,i)) !Inter-cell hopping
		END IF
		
		END IF

	


	END DO

	Ham(counter,counter)=Ham(counter,counter)+E(counter)



	counter = counter + 1


	END DO

	!===============================================================================================================
	!=====================================================================================Second Nearest Neighbor==
	!===============================================================================================================

	counter = 1

	DO j=(1+(lreal*Lwidth)),(L-(lreal*Lwidth))  !loop over atoms per unit

		DO i=1,L !constructing the hamiltonian
		
		IF ((i/=j).AND.(nn2(j,i)/=0))THEN

		rij = a*sqrt((pos(i,1)-pos(j,1))**2 + (pos(i,2)-pos(j,2))**2) !displacement for the exponent

		!strain = ((sqrt(3d0)*a)/(real(rij)))**2 !Harrison
		strain = exp( -b*((real(rij)/sqrt(3d0))-1) )  !Exponential	

		IF (nn2(j,i)<0)THEN !Intra-cell hopping
			IF ((pos(j,2)-pos(i,2))<0)THEN
				Ham(counter,-nn2(j,i))=Ham(counter,-nn2(j,i))-strain*t2(counter,-nn2(j,i))*exp(Im*k*a)
			ELSE		
				Ham(counter,-nn2(j,i))=Ham(counter,-nn2(j,i))-strain*t2(counter,-nn2(j,i))*exp(-Im*k*a)
			END IF

		ELSEIF (nn2(j,i)>0)THEN
			Ham(counter,nn2(j,i))=Ham(counter,nn2(j,i))-strain*t2(counter,nn2(j,i))!Inter-cell hopping
		END IF

		END IF

	END DO

	counter = counter + 1

	END DO

	!=============================================================================================================
	!=====================================================================================Third Nearest Neighbor==
	!=============================================================================================================
	
	counter = 1

	DO j=(1+(lreal*Lwidth)),(L-(lreal*Lwidth))  !loop over atoms per unit

		DO i=1,L !constructing the hamiltonian
		
		IF ((i/=j).AND.(nn3(j,i)/=0))THEN

		rij = a*sqrt((pos(i,1)-pos(j,1))**2 + (pos(i,2)-pos(j,2))**2) !displacement for the exponent

		!strain = ((2*a)/(real(rij)))**2 !Harrison
		strain = exp( -b*((real(rij)/2d0)-1) )  !Exponential	


		IF (nn3(j,i)<0)THEN !Intra-cell hopping
			IF ((pos(j,2)-pos(i,2))<0)THEN		
				Ham(counter,-nn3(j,i))=Ham(counter,-nn3(j,i))-strain*t3(counter,-nn3(j,i))*exp(Im*k*a)
			ELSE
				Ham(counter,-nn3(j,i))=Ham(counter,-nn3(j,i))-strain*t3(counter,-nn3(j,i))*exp(-Im*k*a)
			END IF

		ELSEIF (nn3(j,i)>0)THEN
			Ham(counter,nn3(j,i))=Ham(counter,nn3(j,i))-strain*t3(counter,nn3(j,i)) !Inter-cell hopping
		END IF

		END IF

	END DO

	counter = counter + 1

	END DO


!=========================================================!
!=		   End Hamiltonian Loop			 =!
!=========================================================!

END SUBROUTINE


!=========================================================!
!=		  Random Number Generator		 =!
!=========================================================!

SUBROUTINE sub_seed()
IMPLICIT NONE

INTEGER :: i,n,ck
INTEGER, ALLOCATABLE :: seed(:)

CALL RANDOM_SEED(size=n)
ALLOCATE(seed(n))
CALL SYSTEM_CLOCK(COUNT=ck)
seed = ck + 37*(/(i-1,i=1,n)/)
CALL RANDOM_SEED(PUT=seed)
DEALLOCATE(seed)

END SUBROUTINE




END MODULE
