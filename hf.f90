program hartree_fock
	
	use utils
	use hfmod

	IMPLICIT NONE

	real*8, allocatable :: h(:,:), s(:,:), f(:,:), p(:,:), dipint(:,:,:)
	logical havedip



	!================================================================================
	!                  	MY DECLARATIONS HERE
	!================================================================================
	
	
	integer 							     :: i, j, k, kk, INFO, LWORK, lendiis
	integer , dimension(:),   allocatable    :: IPIV
	real(dp) 							     :: error, dlange, nrg, ddot
	real(dp), dimension(:),   allocatable    :: W, WORK, diisv
	real(dp), dimension(:,:), allocatable    :: U, S_half, Shalf, temp1, temp2, h_til, &
											  & C, C_k, C_k1, P_k, P_k1, Fk_1, B, M
	real(dp), dimension(:,:,:), allocatable  :: F_tens, E_tens
	
	logical  ::  damping 
    logical  ::  diis
  
	
	! USEFUL PARAMETERS: all in utils.f90
	!
    ! tol_scf -->  Convergence threshold
	! itmax   -->  Max SCF iterations
	! f_damp  -->  Damping factor            ( if = 0 : Damping skipped )
	! ndiis   -->  Max dimension DIIS buffer ( if = 0 : DIIS skipped    )
    !
    !	
	!================================================================================
	!================================================================================
	!
	call init
	!
	! print some information
	!
	call print_header
	!
	! allocate memory and then read the one-electron integrals.
	!
	allocate (h(nbas,nbas), s(nbas,nbas), f(nbas,nbas), p(nbas,nbas), dipint(nbas,nbas,3))
	!
	call read_h(h,s,dipint,havedip)
	!
	call mkfock(h,p,f)
	!
	!===============================================================================
	! 			MY PROGRAM FROM HERE
	!===============================================================================

	
	allocate( W(nbas), U(nbas,nbas), S_half(nbas, nbas), Shalf(nbas, nbas), &
            & temp1(nbas, nbas), temp2(nbas, nbas), h_til(nbas, nbas), &
            & C(nbas, nbas), C_k(nbas,nbas), C_k1(nbas,nbas), P_k(nbas,nbas),&
			& P_k1(nbas,nbas), Fk_1(nbas, nbas), &
			& F_tens(nbas, nbas, ndiis ), E_tens(nbas, nbas, ndiis), M(nbas,nbas) )

	!  Variable description:
	!  S_half --> matrix S^-1/2
	!  Shalf  --> matrix S^1/2
	!  h_til  --> S^-1/2*h*S^1/2

	
	
	! ###  MAKING THE S^-1/2 MATRIX  	###
	
	! Copying the S matrix
	U = S


	! WRITING S
	!write(6,*) '---  MATRICE S ---'
	!call output(S, 1, nbas, 1, nbas, nbas, nbas, 1)


	! Dummy call of dsyev
	allocate (WORK(1))
	LWORK = -1

	call dsyev('V', 'U', nbas, U, nbas, W, WORK, LWORK, INFO)

	LWORK = int(WORK(1))

	deallocate(WORK)
	allocate( WORK(LWORK) )


	! Diagonalizing S
	call dsyev('V', 'U', nbas, U, nbas, W, WORK, LWORK, INFO)
	deallocate(WORK)


	! Making S_half = U*L^-1/2*U^t
	do i = 1, nbas

  		temp1(:,i) = (1/sqrt(W(i)))*U(:,i)
  		temp2(:,i) = sqrt(W(i))*U(:,i)
	
	end do

	call dgemm('N', 'T', nbas, nbas, nbas, 1.0_dp, temp1, nbas, U, nbas, 0.0_dp, S_half, nbas)
	call dgemm('N', 'T', nbas, nbas, nbas, 1.0_dp, temp2, nbas, U, nbas, 0.0_dp, Shalf,  nbas)


	! WRITING S^-1/2
	!write(6,*) '---  MATRICE S^-1/2 ---'
	!call output(S_half, 1, nbas, 1, nbas, nbas, nbas, 1)
	
	!###############################

	
	! ### MAKING THE INITIAL GUESS ###

	! WAY 1:
	! Solving hC = SCE
	! h already contains the monoelectronic hamiltonian
	! (subroutine read_h do so)
	call gendiag(nbas, h, s_half, C)

	!WAY 2:
	! Solving MC = SCE
	!do i = 1, nbas
	!	do j = 1, nbas
	!		
	!		if (i .eq. j ) then 
	!			M(i,i) = h(i,i)
	!		else
	!			M(j,i) = 0.5_dp*1.75_dp*S(j,i)*( h(i,i) + h(j,j) )
	!		end if
	!					
	!	end do
	!end do
	!
	!call gendiag(nbas, M, s_half, C)



	! WRITING C (guess mo)
	!write(6,*) '---  MATRICE C ---'
	!call output(C, 1, nbas, 1, nbas, nbas, nbas, 1)

  	!################################


    ! INITILIZATION SOME VARIABLES:

	k       = 0       		 ! Index SCF iterations
	error   = 1.0_dp  		 ! Error
	F       = 0.0_dp  		 ! Fock matrix before entering iterations
	Fk_1    = h       		 ! Core hamiltonian for first damping
    C_k     = C       		 ! MO coefficents: from initial guess
	P_k1    = 0.0_dp  		 ! Density matrix for first error
    nocc    = nele/2  		 ! Number of occupied MOs

	diis    = ndiis .gt. 0   ! Decide whether DIIS will be done or not
	damping = f_damp .gt. 0  !  Decide whether Damping will be done or not
    kk      = 0              ! Index DIIS iterations
    lendiis = 0              ! DIIS buffer dimension
	E_tens  = 0.0_dp         ! DIIS error tensor
	F_tens  = 0.0_dp         ! DIIS Fock tensor

	


	! USEFUL OUTPUT 
	write(6,*) ''
	write(6,*) ''
	write(6,*) '====  ENTERING SCF ITERATIONS ===='
	write(6,*) ''
	write(6,'(a, ES7.1)') ' SCF threshold  = ', tol_scf
	write(6,'(a, I3)')    ' Max iterations = ', itmax
	write(6,'(a, L1)')    ' DAMPING        = ', damping
	write(6,'(a, F5.3)')  ' Damping factor = ', F_DAMP
	write(6,'(a, L1)')    ' DIIS           = ', diis
	write(6,'(a, I1)')    ' DIIS start     = ', diis_start
	write(6,'(a, I2)')    ' DIIS dimension = ', ndiis
	write(6,*) ''
	

	
	
	do while ( (error .gt. tol_scf) .and. (k .lt. itmax) )
    	
		k = k+1
		
		write(6,*) '================================================================='
		write(6,*) '			 ITERATION ', k
		write(6,*) '================================================================='
     	write(6,*) ''
		
		
		! Making density matrix
		call dgemm('N', 'T', nbas, nbas, nocc, 2.0_dp, C_k, nbas, C_k, nbas, 0.0_dp, P_k, nbas) 
		
	
		! Printing density matrix
		!write(6,*) '--- DENSITY MATRIX ---'
		!call output(P_k, 1, nbas, 1, nbas, nbas, nbas, 1)
		!write(6,*) '----------------------'	
		!write(6,*) ''
		
		
		! Calculating the error
		P_k1 = P_k - P_k1
		error = dlange('F', nbas, nbas, P_k1, nbas)
		
		
		! Building Fock matrix
		call mkfock(h, P_k, F)
		

		! Calculating SCF Energy
		temp1  = F + h
        nrg    = 0.50_dp * ddot(nbas*nbas, temp1, 1, P_k, 1) + repnuc
		

			

		! DAMPING
		if( (damping .eqv. .true.) .and. (k .le. diis_start) .and. (diis .eqv. .true.) ) then
	
			F = (1.0_dp - F_DAMP)*F + F_DAMP*Fk_1

		end if
		
		Fk_1 = F
		! DIIS
		if ((diis .eqv. .true.) .and. (k .gt. diis_start) ) then
			
			
			!NEW WAY: WORKING
			if ( kk .eq. ndiis) then
				kk = 0
				write(6,*) 'DIIS BUFFER FULL: RESTART'
			end if

			lendiis = lendiis + 1
			lendiis = min(lendiis, ndiis)
			kk = kk + 1
			
			write(6,'(X,a,i4,a,i4,a,i4)') 'k = ',k,' kk = ', kk, ' lendiis = ', lendiis
			
			
			! F matrix in F tensor
			F_tens(:, :, kk) = F

			
			! Error = FPS - SPF in E tensor
			call dgemm('N', 'N', nbas, nbas, nbas, 1.0_dp, P_k, nbas, S,     nbas, 0.0_dp, temp1, nbas) 
			call dgemm('N', 'N', nbas, nbas, nbas, 1.0_dp, F,   nbas, temp1, nbas, 0.0_dp, temp2, nbas) 
			E_tens(:, :, kk) = temp2		
			call dgemm('N', 'N', nbas, nbas, nbas, 1.0_dp, P_k, nbas, F,     nbas, 0.0_dp, temp1, nbas) 
			call dgemm('N', 'N', nbas, nbas, nbas, 1.0_dp, S,   nbas, temp1, nbas, 0.0_dp, temp2, nbas) 

			E_tens(:, :, kk) = E_tens(:, :, kk) - temp2

			
			! Allocate DIIS variables
			allocate( B(lendiis+1, lendiis+1), diisv(lendiis+1), IPIV(lendiis+1) )

			
			! Building B matrix and diisv vector
			diisv = 0.0_dp
			diisv(lendiis+1) = 1.0_dp
			B = 1.0_dp
			do i = 1, lendiis
				do j = 1, i
					B(i,j) = ddot(nbas*nbas, E_tens(:,:,i), 1, E_tens(:,:,j), 1)
					B(j,i) = B(i,j)
				end do
			end do
			B(lendiis+1,lendiis+1) = 0.0_dp 
	
			
			! solving B*E = diisv system
            call dgesv( lendiis+1, 1, B, lendiis+1, IPIV, diisv, lendiis+1, info)
			write(6,*) 'EXIT DGESV =', info	
           ! write(6,'(20f8.4)') diisv(1:lendiis)
           ! write(6,*) sum(diisv(1:lendiis))
			
			F=0.0_dp
 			do i=1,lendiis
				F= F+diisv(i)*F_tens(:,:,i)
			end do
 
			deallocate ( B,diisv,IPIV )

		end if
		
		! Printing Fock matrix
		!write(6,*) '--- FOCK MATRIX ------'
		!call output(F, 1, nbas, 1, nbas, nbas, nbas, 1)
       	!write(6,*) '----------------------'
		!write(6,*) '' 
 
		
		! k-th iteration output		
		write(6,*) 'SCF ENERGY (Hartree) =', nrg
		write(6,*) 'Error                =', error
		write(6,*) ''
		write(6,*) ''

		! clean output
		write(2,*) k, nrg, error
		

		! Solving Roothan equations
		call gendiag(nbas, F, s_half, C_k)
		

		! Updating the density matrix
		P_k1 = P_k

	end do

	
	! FINAL OUTPUT
		write(6,*) ''
		write(6,*) '================================================================='
	if (error .le. tol_scf) then
		write(6,*) 'SCF CONVERGED IN', k, 'ITERATIONS'
	else
		write(6,*) 'SCF NOT CONVERGED IN', k, 'ITERATIONS'
	end if
		write(6,*) 'SCF ENERGY (Hartree) = ', nrg
		write(6,*) 'Error                = ', error
		write(6,*) '================================================================='
		write(6,*) ''

	!================================================================================
	!
	!free the memory:
	!
	deallocate( W, U, S_half, Shalf, temp1, temp2, h_til, C, C_k, C_k1, P_k,&
			 &  P_k1, Fk_1, F_tens, E_tens)
	deallocate (h,s,f,p)
	!
end program hartree_fock
