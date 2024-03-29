MODULE Construct_ZZ_Hamiltonian
IMPLICIT NONE

CONTAINS


!=========================================================!
!=	         Hamiltonian Subroutine	ZZ		 =!
!=========================================================!

SUBROUTINE HamiltonianZZ(wmax,lreal1,pos,nn1,nn2,nn3,t1,t2,t3,N,counter,k,a,Lwidth,L,Ham,E ,spin_density,Hubbard_U,n_spin_old,g,EP)
IMPLICIT NONE

INTEGER :: i,j,rr
DOUBLE PRECISION :: harrison,r, t_edge
INTEGER, INTENT(INOUT) :: wmax,lreal1,Lwidth,a,L,N,counter
INTEGER, INTENT(IN) :: g,EP
INTEGER, DIMENSION(L,L), INTENT(INOUT) :: nn1,nn2,nn3
DOUBLE PRECISION, DIMENSION(L,2), INTENT(INOUT) :: pos
DOUBLE PRECISION, DIMENSION(N), INTENT(INOUT) :: spin_density,Hubbard_U,n_spin_old
COMPLEX*16, DIMENSION(N,N), INTENT(INOUT) :: Ham,t1,t2,t3
DOUBLE PRECISION, DIMENSION(N), INTENT(INOUT) :: E
COMPLEX*16, INTENT(INOUT) :: k
COMPLEX*16 :: Im,rij


Im = (0.0,1.0)


!=========================================================!
!=		    Hamiltonian Loop ZZ			 =!
!=========================================================!

	IF (EP==1)THEN!edge pertubations
		t_edge = 2.808
	
		t1(1,2) = t_edge
		t1(2,1) = t_edge

		t1(N,N-1) = t_edge
		t1(N-1,N) = t_edge
	END IF
	

	!===================================================================================== First Nearest Neighbor
	
	DO j=(1+(2+(wmax*2))*lreal1),(L-(2+(wmax*2))*lreal1) !loop over atoms per unit

		DO i=1,L !constructing the hamiltonian
		
		IF ((i/=j).AND.(nn1(j,i)/=0))THEN

		rij = sqrt((pos(i,1)-pos(j,1))**2 + (pos(i,2)-pos(j,2))**2) !displacement for the exponent

		harrison = (a/(real(rij)))**2

		IF (nn1(j,i)<0)THEN !Intra-cell hopping
			IF ((pos(j,2)-pos(i,2))<0)THEN
				Ham(counter,-nn1(j,i))=Ham(counter,-nn1(j,i))-harrison*t1(counter,-nn1(j,i))*exp(Im*k*a)
			ELSE
				Ham(counter,-nn1(j,i))=Ham(counter,-nn1(j,i))-harrison*t1(counter,-nn1(j,i))*exp(-Im*k*a)
			END IF

		ELSEIF (nn1(j,i)>0)THEN
			Ham(counter,nn1(j,i))=Ham(counter,nn1(j,i))-harrison*t1(counter,nn1(j,i)) !Inter-cell hopping
		END IF

		END IF

	END DO


	IF (g==1)THEN
	Ham(counter,counter) = Ham(counter,counter) + Hubbard_U(counter)*spin_density(counter)
	ELSE
	Ham(counter,counter) = Ham(counter,counter) + Hubbard_U(counter)*n_spin_old(counter)
	END IF

	counter = counter + 1

	END DO

	!===================================================================================== Second Nearest Neighbor

	counter = 1

	DO j=(1+(2+(wmax*2))*lreal1),(L-(2+(wmax*2))*lreal1) !loop over atoms per unit

		DO i=1,L !constructing the hamiltonian
		
		IF ((i/=j).AND.(nn2(j,i)/=0))THEN

		rij = a*sqrt((pos(i,1)-pos(j,1))**2 + (pos(i,2)-pos(j,2))**2) !displacement for the exponent

		harrison = ((sqrt(3d0)*a)/(real(rij)))**2

		IF (nn2(j,i)<0)THEN !Intra-cell hopping
			IF ((pos(j,2)-pos(i,2))<0)THEN

				IF (real(rij)>(2.6*a))THEN
				Ham(counter,-nn2(j,i))=Ham(counter,-nn2(j,i))-harrison*t2(counter,-nn2(j,i))*exp(2*Im*k*a)
				ELSE
				Ham(counter,-nn2(j,i))=Ham(counter,-nn2(j,i))-harrison*t2(counter,-nn2(j,i))*exp(Im*k*a)
				END IF

			ELSE
				
				IF (real(rij)>(2.6*a))THEN
				Ham(counter,-nn2(j,i))=Ham(counter,-nn2(j,i))-harrison*t2(counter,-nn2(j,i))*exp(-2*Im*k*a)
				ELSE
				Ham(counter,-nn2(j,i))=Ham(counter,-nn2(j,i))-harrison*t2(counter,-nn2(j,i))*exp(-Im*k*a)
				END IF


			END IF

		ELSEIF (nn2(j,i)>0)THEN
			Ham(counter,nn2(j,i))=Ham(counter,nn2(j,i))-t2(counter,nn2(j,i))!Inter-cell hopping
		END IF

		END IF

	END DO

	counter = counter + 1

	END DO


	!===================================================================================== Third Nearest Neighbor

	counter = 1

	DO j=(1+(2+(wmax*2))*lreal1),(L-(2+(wmax*2))*lreal1) !loop over atoms per unit

		DO i=1,L !constructing the hamiltonian
		
		IF ((i/=j).AND.(nn3(j,i)/=0))THEN

		rij = a*sqrt((pos(i,1)-pos(j,1))**2 + (pos(i,2)-pos(j,2))**2) !displacement for the exponent

		harrison = ((2*a)/(real(rij)))**2

		IF (nn3(j,i)<0)THEN !Intra-cell hopping
			IF ((pos(j,2)-pos(i,2))<0)THEN


				IF (real(rij)>(2.6*a))THEN
				Ham(counter,-nn3(j,i))=Ham(counter,-nn3(j,i))-harrison*t3(counter,-nn3(j,i))*exp(2*Im*k*a)
				ELSE
				Ham(counter,-nn3(j,i))=Ham(counter,-nn3(j,i))-harrison*t3(counter,-nn3(j,i))*exp(Im*k*a)
				END IF

			ELSE

				IF (real(rij)>(2.6*a))THEN
				Ham(counter,-nn3(j,i))=Ham(counter,-nn3(j,i))-harrison*t3(counter,-nn3(j,i))*exp(-2*Im*k*a)
				ELSE
				Ham(counter,-nn3(j,i))=Ham(counter,-nn3(j,i))-harrison*t3(counter,-nn3(j,i))*exp(-Im*k*a)
				END IF

			END IF

		ELSEIF (nn3(j,i)>0)THEN
			Ham(counter,nn3(j,i))=Ham(counter,nn3(j,i))-t3(counter,nn3(j,i)) !Inter-cell hopping
		END IF

		END IF

	END DO

	counter = counter + 1

	END DO

	!===================================================================================== 

!=========================================================!
!=		   End Hamiltonian Loop			 =!
!=========================================================!

END SUBROUTINE








END MODULE
