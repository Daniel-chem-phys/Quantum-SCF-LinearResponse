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

 !declaring some additional variables (polarizability)
	
	integer 							     :: a, d, ia, jb, r, t, cnt, nvirt
	integer, dimension(:), allocatable       :: IPIV_p
	real(dp) 							     :: norm, P_DAMP
	real(dp), dimension(:),   allocatable    :: sigma_mo_vec, C_vi, C_va, temp4, diisv_p
	real(dp), dimension(:,:), allocatable    :: X_0, X, X_k, X_ao, Y, Y_ao, G, G_mo, temp3, sigma_mo_mat, incr, B_p, pol									
	real(dp), dimension(:,:,:), allocatable  :: Q_mo, Q_ov, X_tens, incr_tens
	logical  ::  damping_p 
    logical  ::  diis_p

  
	
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


    !-----------------------------------------------------------------------------------------------------------------------
    !                                     POLARIZABILITY
	!-----------------------------------------------------------------------------------------------------------------------

	! define some useful variables
	nvirt = nbas-nocc

	! initializing variables
	P_DAMP = 0.1_dp
	damping_p = .false.
	diis_p = .true.


    ! allocate polarizability variables
	allocate( sigma_mo_vec(nocc*nvirt), C_vi(nbas), C_va(nbas), temp4(nocc*nvirt) )
    allocate( X_0(nocc,nvirt), X(nocc,nvirt), X_k(nocc,nvirt), Y(nbas,nbas), Y_ao(nbas,nbas), &
	& G(nbas,nbas), G_mo(nbas,nbas), temp3(nbas,nbas), sigma_mo_mat(nocc,nvirt), incr(nocc,nvirt), &
	& pol(3,3))
	allocate(Q_mo(nbas,nbas,3), Q_ov(nocc,nvirt,3), X_tens(nocc,nvirt,300), incr_tens(nocc,nvirt,300))

    Q_ov = 0.0_dp
	Q_mo = 0.0_dp


	! write(6,*) 'la matrice degli integrali di dipolo è'
	! do i = 1, nbas
	!	 write (6,*) dipint(i,:,1)
	! end do 
	! write (6,*) '-------------------------------------------------------------------------------------'

	! trasforming Q(x,y,z): AO-->MO (C^T * Q * C)

	do j=1,3
		call dgemm('T', 'N', nbas, nbas, nbas, 1.0_dp, C_k, nbas, dipint(:,:,j), nbas, 0.0_dp, temp3, nbas)
    	call dgemm('N', 'N', nbas, nbas, nbas, 1.0_dp, temp3, nbas, C_k, nbas, 0.0_dp, Q_mo(:,:,j) , nbas)

	 	! building matrix Q_ov (occupied-virtual block)
		do i=1, nocc
			do a=1, nvirt
				q_ov(i,a,j) = Q_mo(i,a+nocc,j)
			end do
		end do

        !write(6,*) 'q_ov ', j

		!do i = 1, nocc
		!  write (6,*) q_ov(i,:,j)
		!end do 

		!write (6,*) '--------------------------------------------------------------------------------------'

	end do

	!write(6,*) 'la q_ov è'

	!do i = 1, nocc
	!  write (6,*) q_ov(i,:,1)
	!end do 

	!write (6,*) '--------------------------------------------------------------------------------------'

	
	pol= 0.0_dp

	!calculating Jacobi's guess (X_0)  	
  	do k=1,3

		X_0= 0.0_dp
	  	X_tens(:,:,:) = 0.0_dp

      	do i=1,nocc
	    	do a=1,nvirt
			  	X_0(i,a) = Q_ov(i,a,k)/(w(nocc+a)-w(i))
		  	end do
	  	end do

	  	!write(6,*) 'il guess X_0 è', k
	
		!do i = 1, nocc
	  		!write (6,*) X_0(i,:)
	   	!end do 

	 	! write(6,*) '----------------------------------------------------------------------------------------------------'			
	  
	  	!------entering Jacobi--------------
	  	cnt = 0
	  	norm = 1.0_dp 
      	do while (norm .ge. 1.0E-06 .and. cnt .le. 200)

        	cnt = cnt+1

         	if (cnt == 1) then
	        	X_k = X_0
         	else
	        	X_k = X
         	end if 

         	! costruzione matrice Y a partire da X all'iterazione precedente
		 	Y = 0.0_dp

         	do i=1,nocc
       	    	do a=1,nvirt
       		    	Y(i,a+nocc)=X_k(i,a)                   
       		    	Y(a+nocc,i)=X_k(i,a)
       	    	end do
         	end do

			!if (cnt == 28) then
			!	write (6,*)'la matrice Y all iterazione 28', k
			!	do i = 1, nbas 
			!   	write (6,*) Y(i,:)
			!	end do 
			!end if 

         	! trasforming Y: MO --> AO (C * Y * C^T)
         	call dgemm('N', 'N', nbas, nbas, nbas, 1.0_dp, C_k, nbas, Y, nbas, 0.0_dp, temp3, nbas)
         	call dgemm('N', 'T', nbas, nbas, nbas, 1.0_dp, temp3, nbas, C_k, nbas, 0.0_dp, Y_ao , nbas)

         	! building bielectronic term with Y_ao,instead of P as in HF
         	call mkfock(h,Y_ao,f)

         	G=f-h

	     	! write (6,*) '-------------------------------------------------------------------------------------------------------'
	     	! write (6,*) 'la matrice G è'
	     	! do i=1,nbas
	     	!    write (6,*) G(:,i)
	     	! end do

	     	! tranforming G: AO-->MO
	     	call dgemm('T', 'N', nbas, nbas, nbas, 1.0_dp, C_k, nbas, G, nbas, 0.0_dp, temp3, nbas)
         	call dgemm('N', 'T', nbas, nbas, nbas, 1.0_dp, temp3, nbas, C_k, nbas, 0.0_dp, G_mo , nbas)

	     	! building sigma_mo_mat from G_mo (OV) 
	     	do i=1,nocc
		    	do a=1,nvirt
               		sigma_mo_mat(i,a) = G_mo(i,a+nocc)
		    	end do
	     	end do

	     	! write (6,*) '-------------------------------------------------------------------------------------------------------'
	     	! do i=1,nocc
	     	! write(6,*) sigma_mo_mat(i,:)
	     	! end do

	     	! sigma(MO) = sigma(MO) - Q(x)
	     	sigma_mo_mat = sigma_mo_mat - Q_ov(:,:,k)

		 	!if (cnt == 28) then
			!	write (6,*)'la matrice sigma all iterazione 28', k
			!	do i = 1, nocc
		 	!   	write (6,*) sigma_mo_mat(i,:)
			!	end do 
			!end if

		 	! write (6,*) '-------------------------------------------------------------------------------------------------------'
		 	! write(6,*) 'sigma_mo_mat- Q è'
		 	! do i=1,nocc
		 	!	   write(6,*) sigma_mo_mat(i,:)
		 	! end do

		 	X = 0.0_dp

	     	! product D^-1 * (sigma(MO) - Q)
	     	do i=1,nocc
		    	do a=1,nvirt
	            	X(i,a)= -sigma_mo_mat(i,a)/(W(nocc+a)-W(i))
		    	end do
	     	end do		

		 	! write (6,*) '-------------------------------------------------------------------------------------------------------'
		 	! write(6,*) 'la matrice X è'
		 	! do i=1,nocc
		 	!	 write(6,*) X(i,:)
		 	! end do

		 	X_tens(:,:, cnt) = X 

         	! DAMPING
	     	if( (cnt .ge. 2) .and. (damping_p .eqv. .true.) ) then
		    	X = (1.0_dp - P_DAMP)*X + P_DAMP*X_tens(:,:,cnt-1)
	     	end if

		 	incr = 0.0_dp
		 	norm = 0.0_dp

		 	if (cnt == 1) then
				incr = X_tens(:,:,cnt) - X_0
		 	else
				incr = X_tens(:,:,cnt) - X_tens(:,:,cnt-1)
		 	end if 

		 	incr_tens(:,:,cnt) = incr 

		 	! write(6,*) '------------------------------------------------------------------------------------------'
		 	! write(6,*) 'la matrice degli incrementi è'
		 	! do i = 1, nocc
		 	!    write (6,*) incr(i,:)
		 	! end do 
		
		 	norm = dlange('F', nocc, nvirt, incr, nocc)

		 	write(6,*) 'la norma all iteraz', cnt , 'è', norm

         	! DIIS
		 	if ((cnt .ge. 20) .and. (diis_p .eqv. .true.)) then
	     		!calculating increment

	     		! if (cnt == 1) then 
		     		!incr = X - X_0
         		! else 
	      			! incr = X_tens(:,:,cnt) - X_tens(:,:,cnt-1)
	      		!end if 
	
	     		! incr_tens(:,:,cnt) = incr 

         		! calculating error (Frobenius norm of incr)
	

		 		! if ( kk .eq. ndiis) then
		 		!	kk = 0
		 		!	write(6,*) 'DIIS BUFFER FULL: RESTART'
		 		! end if

		 		! lendiis = lendiis + 1
		 		! lendiis = min(lendiis, ndiis)
		 		! kk = kk + 1

		
		 		! Allocate DIIS variables

		 		allocate( B_p(cnt+1, cnt+1), diisv_p(cnt+1), IPIV_p(cnt+1) )

		 		! Building B matrix and diisv vector

		 		diisv_p = 0.0_dp
		 		diisv_p(cnt+1) = 1.0_dp
		 		B_p = 1.0_dp
		 		do i = 1, cnt
			 		do j = 1, i
				 		B_p(i,j) = ddot((nocc*nvirt), incr_tens(:,:,i), 1, incr_tens(:,:,j), 1)
				 		B_p(j,i) = B_p(i,j)
			 		end do
		 		end do

		 		B_p(cnt+1,cnt+1) = 0.0_dp 

		
				! solving B*E = diisv system
		 		call dgesv( cnt+1, 1, B_p, cnt+1, IPIV_p, diisv_p, cnt+1, info)
		 		write(6,*) 'EXIT DGESV =', info	
 
		
		 		X=0.0_dp
		 		do i=1,cnt
	 	    		X = X+diisv_p(i)*X_tens(:,:,i)
		 		end do

		 		deallocate ( B_p, diisv_p, IPIV_p )

	     	end if

		end do
		
        ! building the polarizability tensor  
		do kk=1,3
			pol(kk,k) = 0.0_dp
			pol(kk,k) = 4.0_dp*ddot(nocc*nvirt, Q_ov(:,:,kk), 1, X, 1)
		end do

        ! output converged X
		write(6,*) 'converged matrix X relative to Q_', k
		do i=1,nocc
			write(6,*) X(i,:)
		end do

    end do

	! output polarizability
	write(6,*) '-----------------polarizability tensor-----------------'
    do i=1,3
    	write(6,*) pol (i,:)
	end do


    !free the memory(HF)
    deallocate( W, U, S_half, Shalf, temp1, temp2, h_til, C, C_k, C_k1, P_k, &
    &  P_k1, Fk_1, F_tens, E_tens)
    deallocate (h,s,f,p)

    !free the memory(polarizability)
    deallocate(sigma_mo_vec, C_vi, C_va, temp4)
    deallocate ( X_0, X, X_k, Y, Y_ao, G, G_mo, temp3, sigma_mo_mat, pol, incr, Q_mo, Q_ov, X_tens, incr_tens)
	
end program hartree_fock