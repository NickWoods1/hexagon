PROGRAM GrapheneTBM

!import the modules
USE Armchair_Geometry
USE Zigzag_Geometry
USE Metal_Armchair_Geometry
USE Construct_ZZ_Hamiltonian
USE Construct_AC_Hamiltonian
USE Construct_MAC_Hamiltonian

IMPLICIT NONE

!=========================================================================================================!
!          ____       							        	     ____         !
!         /    \        ██╗  ██╗███████╗██╗  ██╗ █████╗  ██████╗  ██████╗ ███╗   ██╗        /    \        !
!    ____/      \____	██║  ██║██╔════╝╚██╗██╔╝██╔══██╗██╔════╝ ██╔═══██╗████╗  ██║    ___/      \____	  !	
!   /    \      /    \  ███████║█████╗   ╚███╔╝ ███████║██║  ███╗██║   ██║██╔██╗ ██║  /    \      /    \  !
!  /      \____/      \ ██╔══██║██╔══╝   ██╔██╗ ██╔══██║██║   ██║██║   ██║██║╚██╗██║ /      \____/      \ !	
!  \      /    \      / ██║  ██║███████╗██╔╝ ██╗██║  ██║╚██████╔╝╚██████╔╝██║ ╚████║ \      /    \      / !
!   \____/      \____/  ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝  ╚═════╝ ╚═╝  ╚═══╝  \____/      \____/  !
!   /    \      /    \	 							      /    \      /    \  !
!  /      \____/      \			             By		  		     /      \____/      \ !			
!  \      /    \      /							             \      /    \      / !
!   \____/      \____/		                   N. D. Woods	                      \____/      \____/  !
!        \      /								           \      /       !
!         \____/								  	    \____/        !
!=========================================================================================================!


INTEGER :: i,j,L,m,ok,g,atm,nrepeat,p,f,b,counter,ii,wmax,lmax,Lwidth,lreal,N,a,lreal1,inpt1,inpt2,inpt3,hub,rr,counter1,EP,gg

DOUBLE PRECISION :: r,ex,ey,varswitch1,varswitch2,shear,Metal_energy_hopping3,Metal_energy_onsite,Metal_energy_hopping1

DOUBLE PRECISION :: Metal_energy_hopping2,Poisson

INTEGER, DIMENSION(:,:),ALLOCATABLE :: atom,nn,nn1,nn2,nn3

DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: pos

DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: E,Uhub,a1,bubble,spin_density,Hubbard_U,n_spin_old,n_spin,n_spin_compare

COMPLEX*16, DIMENSION(:,:),ALLOCATABLE :: Ham,EigVec,EigVec1,t1,t2,t3

COMPLEX*16, DIMENSION(:), ALLOCATABLE :: Eigen,WORK,nup,RWORK

COMPLEX*16, DIMENSION(1,1) :: DUMMY

COMPLEX*16 :: Im,z,U,rij,k


OPEN(UNIT=1,FILE="Edisp.dat")
OPEN(UNIT=2,FILE="positions.dat")
OPEN(UNIT=3,FILE="EigVec.dat")
OPEN(UNIT=20, FILE="TBM.in")
OPEN(UNIT=40,FILE="Strain.plot",ACCESS = 'append',STATUS='old')

!=========================================================!
!=		    Defining Code Parameters		 =!
!=========================================================!

Im = (0.0,1.0)
k = (0.0,0.0)
a = 1.0d0 !Lattice Length
counter = 1
counter1 = 1
Poisson = 0.186 !Poisson ratio



!=========================================================!Parse the input file

READ(20,*)
READ(20,*) inpt1

READ(20,*)
READ(20,*)
READ(20,*) wmax

READ(20,*)
READ(20,*)
READ(20,*) lreal

READ(20,*)
READ(20,*)
READ(20,*) inpt2

READ(20,*)
READ(20,*)
READ(20,*) ex

READ(20,*)
READ(20,*)
READ(20,*) EP

READ(20,*)
READ(20,*)
READ(20,*) shear
!=========================================================



!=========================================================!
!=		    User Selected Metal AC		 =!
!=========================================================!

IF (inpt1==3)THEN

	lmax = lreal + 2

	IF (mod(wmax,2)==0)THEN
		L = 7 + (0.5*(wmax))*3 + (0.5*(wmax)-1)*2
		ELSEIF (mod(wmax,2)/=0)THEN
		L = 7 + (0.5*(wmax)-1)*3 + (wmax+1)
	END IF

	Lwidth = L
 	N = L

	PRINT*, "Number of Atoms in the Unit Cell:", N

	ALLOCATE(nn(f,L),Ham(L,L),E(L),Eigen(L),WORK(2*L),EigVec(L,L),bubble(201),EigVec1(L,201*L),nup(L),Uhub(L),t1(L,L))
	ALLOCATE(RWORK(2*L),t2(L,L),t3(L,L),a1(2))

	L = N*3

	ALLOCATE(pos(L,2),nn1(L,L),nn2(L,L),nn3(L,L))

	CALL Metal(wmax,lreal,lreal1,b,a,pos,nn1,nn2,nn3,L,lmax,Lwidth,a1,N)

END IF


!=========================================================!
!=		    	 End Metal AC			 =!
!=========================================================!


!=========================================================!
!=		    User Selected ZigZag		 =!
!=========================================================!

IF (inpt1==2)THEN

	IF (lreal==1)THEN
		lreal1 = 2
		ELSE
		lreal1 = lreal
	END IF

	N = ((2+2*wmax)*lreal)
	L = N

	PRINT*, "Number of Atoms in the Unit Cell:", N

	ALLOCATE(nn(f,L),Ham(L,L),E(L),Eigen(L),WORK(2*L),EigVec(L,L),bubble(201),EigVec1(L,201*L),nup(L),Uhub(L),t1(L,L))
	ALLOCATE(RWORK(2*L),t2(L,L),t3(L,L),a1(2),spin_density(L),Hubbard_U(L),n_spin(L),n_spin_old(L))
	ALLOCATE(n_spin_compare(L))

	IF (lreal==1)THEN
		lmax = lreal +  4
		ELSE
		lmax = 3*lreal
	END IF

	Lwidth = 4+(wmax-1)*2
	L = lmax*Lwidth

	ALLOCATE(pos(L,2),nn1(L,L),nn2(L,L),nn3(L,L))

CALL Zigzag(wmax,lreal,lreal1,b,a,pos,nn1,nn2,nn3,L,lmax,Lwidth,a1,shear,ex,Poisson)

END IF

!=========================================================!
!=		    	 End ZigZag			 =!
!=========================================================!


!=========================================================!
!=		    User Selected Armchair		 =!
!=========================================================!

IF (inpt1==1)THEN

	lmax = 3*lreal
	N = ((4+2*wmax)*lreal)
	L = N

	PRINT*, "Number of Atoms in the Unit Cell:", N

	ALLOCATE(nn(f,L),Ham(L,L),E(L),Eigen(L),WORK(2*L),EigVec(L,L),bubble(201),EigVec1(L,201*L),nup(L),Uhub(L),t1(L,L))
	ALLOCATE(RWORK(2*L),t2(L,L),t3(L,L),a1(2))


	Lwidth = 6+(wmax-1)*2
	L = lmax*Lwidth

	ALLOCATE(pos(L,2),nn1(L,L),nn2(L,L),nn3(L,L))

CALL Armchair(wmax,lreal,lreal1,b,a,pos,nn1,nn2,nn3,L,lmax,Lwidth,a1,shear,ex,Poisson)

END IF

!=========================================================!
!=		    	 End Armchair			 =!
!=========================================================!


!=========================================================!
!=             		Parameterization		 =!
!=========================================================!


!Unperturbed Hopping Energies
IF (inpt2==3)THEN
	t1 = (2.7,0.0)
	t2 = (0.2,0.0)
	t3 = (0.18,0.0)
ELSEIF (inpt2==2)THEN
	t1 = (2.7,0.0)
	t2 = (0.2,0.0)
ELSE
	t1 = (2.7,0.0)
END IF

!=========================================================!
!=             	Metal Doped Parameterization		 =!
!=========================================================!

Metal_energy_onsite = 3.3
E(:) = 0

IF (inpt1==3)THEN!if the metal-doped structure has been chosen

	E(7)=Metal_energy_onsite

	DO i=7,N
		IF (mod(i,5)==0)THEN
			E(i)=Metal_energy_onsite
		END IF

	END DO

END IF


!Perturb the hopping energies
Metal_energy_hopping3 = 0.05
Metal_energy_hopping1 = 1.25
Metal_energy_hopping2 = 0.1

IF (inpt1==3)THEN!if the metal-doped structure has been chosen

	t3(7,:)= Metal_energy_hopping3
	t3(:,7)= Metal_energy_hopping3

	t1(7,:)= Metal_energy_hopping1
	t1(:,7)= Metal_energy_hopping1

	t2(7,:)= Metal_energy_hopping2
	t2(:,7)= Metal_energy_hopping2

	DO i=7,N
		IF (mod(i,5)==0)THEN
			t3(i,:)= Metal_energy_hopping3
			t3(:,i)= Metal_energy_hopping3

			t1(i,:)= Metal_energy_hopping1
			t1(:,i)= Metal_energy_hopping1

			t2(i,:)= Metal_energy_hopping2
			t2(:,i)= Metal_energy_hopping2

		END IF

	END DO

END IF


!=========================================================!
!=                 End Parameterization			 =!
!=========================================================!



!=========================================================!
!=		    	k-space Loop			 =!
!=========================================================!

IF (inpt1/=2)THEN

DO g=0,200 !kloop

	counter = 1
	Ham=(0.0,0.0)


	IF (inpt1==1)THEN
		CALL HamiltonianAC(wmax,lreal,lreal1,pos,nn1,nn2,nn3,t1,t2,t3,N,counter,Im,k,a,Lwidth,L,Ham,rij,E,nup,Uhub,EP)!Generate Hamiltonian
	ELSEIF (inpt1==3)THEN
		CALL HamiltonianMAC(wmax,lreal,lreal1,pos,nn1,nn2,nn3,t1,t2,t3,N,counter,Im,k,a,Lwidth,L,Ham,rij,E,nup,Uhub,EP)
	END IF


	CALL ZGEEV('N', 'V', N, Ham, N, Eigen, DUMMY, 1, EigVec, N, WORK, 2*L, RWORK, ok)!Diagonalize H


DO j=1,(N**2)
	DO i=1,N-1

		IF (real(Eigen(i))<real(Eigen(i+1)))THEN
			varswitch1=Eigen(i)
			varswitch2=Eigen(i+1)

			Eigen(i)=varswitch2
			Eigen(i+1)=varswitch1

		END IF

	END DO
END DO

	bubble(counter1) = real(Eigen(N/2)) - real(Eigen((N/2)+1))

	counter1=counter1+1

	DO i=1,N

		WRITE(1,*)(0.255/3.14159)*Real(k),(Real(Eigen(i))-2.85672)
		WRITE(3,*)real(EigVec(:,i))

	END DO

	k=k+(0.01570795,0.0)

END DO

PRINT*, "Band Gap (eV) =", minval(real(bubble))!bandgap

WRITE(40,*) wmax+2, minval(real(bubble))

END IF


IF (inpt1==2)THEN!User chose zigzag, need hubbard

!=========================================================!
!=             Self-Consistant Hubbard Loop	         =!
!=========================================================!

spin_density = 0
Hubbard_U = 0

DO g=1,1000

k=(0.0,0.0)

DO i=1,50!Adiabatic U implementation
	CALL random_number(r)
	rr=int(((N))*r)+1
	IF (Hubbard_U(rr)<2)THEN
		Hubbard_U(rr)=Hubbard_U(rr)+0.2
		EXIT
	END IF
END DO


DO gg=0,200 !kloop

counter = 1

Ham=(0.0,0.0)

CALL HamiltonianZZ(wmax,lreal1,pos,nn1,nn2,nn3,t1,t2,t3,N,counter,k,a,Lwidth,L,Ham,E,spin_density,Hubbard_U,n_spin_old,g,EP)!Create Hamiltonian

CALL ZGEEV('N', 'V', N, Ham, N, Eigen, DUMMY, 1, EigVec, N, WORK, 2*L, RWORK, ok)!Diagonalize H

spin_density = 0

DO j=1,(N/2)
	DO i=1,N
		spin_density(i) = spin_density(i) + EigVec(i,minloc(real(Eigen),1))*CONJG(EigVec(i,minloc(real(Eigen),1)))!calculate average occupancies
	END DO
	Eigen(minloc(real(Eigen))) = Eigen(minloc(real(Eigen))) + 100
END DO

DO i=1,N
	n_spin(i) = n_spin(i) + (1d0/201d0)*spin_density(i)!normalize over k-space
END DO

k=k+(0.01570795,0.0)

END DO

!=========================================================!
!=      	   Convergence Condition	         =!
!=========================================================!

IF (mod(g,2)==0)THEN

	IF (ALL( abs( (real(n_spin_compare) - real(n_spin_old)) ) < 0.0000001 ))THEN
		PRINT*,"Average spin polarization on edge atom:", abs(real(n_spin(1)) - real(n_spin_old(1)))
		EXIT
	END IF

	n_spin_compare(:) = n_spin_old(:)

END IF

	n_spin_old(:) = n_spin(:)!Recycle old arrays
	n_spin = 0

END DO

!=========================================================!
!=      	     End Hubbard Loop	  	         =!
!=========================================================!



!=========================================================!
!=      	  Plot the Converged Output	         =!
!=========================================================!

k=(0.0,0.0)

DO gg=0,200

counter = 1
Ham=(0.0,0.0)

CALL HamiltonianZZ(wmax,lreal1,pos,nn1,nn2,nn3,t1,t2,t3,N,counter,k,a,Lwidth,L,Ham,E,spin_density,Hubbard_U,n_spin_old,g,EP)

CALL ZGEEV('N', 'V', N, Ham, N, Eigen, DUMMY, 1, EigVec, N, WORK, 2*L, RWORK, ok)!Diagonalize H

DO j=1,(N**2)
	DO i=1,N-1

		IF (real(Eigen(i))<real(Eigen(i+1)))THEN
			varswitch1=Eigen(i)
			varswitch2=Eigen(i+1)

			Eigen(i)=varswitch2
			Eigen(i+1)=varswitch1

		END IF

	END DO
END DO

bubble(counter1) = real(Eigen(N/2)) - real(Eigen((N/2)+1))

counter1 = counter1 + 1


DO i=1,N
	WRITE(1,*)Real(k),Real(Eigen(i))
END DO

k=k+(0.01570795,0.0)

END DO

PRINT*, "Direct Band Gap (eV) =", minval(real(bubble))!bandgap

WRITE(40,*) ex, minval(real(bubble))

END IF


!=========================================================!
!=		   	   End				 =!
!=========================================================!


END PROGRAM
