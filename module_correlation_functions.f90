module module_correlation_functions
    
! =================================================================================================================
! GLOBAL VARIABLES
! =================================================================================================================

character(*),parameter  :: module_version = '1.16'
integer,allocatable     :: tstart(:)            ! variables for measuring the computation time
integer                 :: tindex,trate         ! variables for measuring the computation time
logical                 :: time_measurement     ! variables for measuring the computation time
logical                 :: show_progress

contains

! =================================================================================================================
! MAIN SUBROUTINES
! =================================================================================================================

subroutine set_opt_progress(opt_progress)
    implicit none
    logical,intent(in)  :: opt_progress
    show_progress = opt_progress
end subroutine set_opt_progress

subroutine compute_xi2(delta_r,L,scale,xi2,scalek,p,scale_min,opt_normalization,opt_verbose)

    ! BRIEF EXPLANATION
    ! Computes the isotropic 2-point correlation function of the 3D density field delta_r, defined on a
    ! cubic box of side-length L and volume V=L^3. Uses a deterministic algorithm that
    ! fully samples the 3D Fourier space.
    !
    ! INPUT VARIABLES
    ! delta_r(:,:,:)        3D density perturbation field of dimension (n,n,n)
    !                       lower and upper bounds are irrelevant
    ! L                     [length unit of density field] side-length of box
    ! opt_normalization     [1,2,3] normalization type of power spectrum (1: Obreschkow, 2: Wolstenhulme, 3: CAMB)
    ! opt_verbose           whether to write output
    !
    ! OUTPUT VARIABLES
    ! scale(:)              [unit of L] correlation scale r
    ! xi2(:)                isotropic two-point correlation function xi_2(r)
    ! scalek(:)             [1/unit of L] wave vector
    ! p(:)                  isotropic power spectrum p(k)

    implicit none
    
    ! variable declaration
    real,intent(in)                         :: delta_r(:,:,:)       ! Density field
    real,intent(in)                         :: L                    ! side-length of density field delta_r
    real,allocatable,intent(out)            :: scale(:)             ! [unit of L] list of scale lengths of correlation function
    real,allocatable,intent(out)            :: xi2(:)               ! values of correlation function
    real,allocatable,intent(out),optional   :: scalek(:)            ! [1/unit of L] wave vector
    real,allocatable,intent(out),optional   :: p(:)                 ! isotropic power spectrum
    real,intent(in)                         :: scale_min
    integer,intent(in),optional             :: opt_normalization    ! normalization type of power spectrum
    logical,intent(in),optional             :: opt_verbose
    real*4,allocatable                      :: ptotsqr(:)           ! total power over the modes of square length |k|^2
    real*4,allocatable                      :: ptot(:)              ! total power over the modes of length |k|
    complex,allocatable                     :: delta_k(:,:,:)       ! Fourier modes of density field delta_r
    integer,parameter                       :: dim = 3              ! number of dimensions
    integer                                 :: n                    ! n^dim = number of elements in delta_k
    integer                                 :: m                    ! number of elements in scale(:)
    real                                    :: dk                   ! linear spacing of Fourier modes
    integer                                 :: i                    ! general counter
    integer                                 :: ik                   ! index of norm of vector k
    integer                                 :: ik1,ik2,ik3          ! indices for the cartesian coordinates of the vector k
    integer                                 :: iksqr,iksqrmax
    integer                                 :: ik3sqr,ik23sqr
    integer                                 :: ikmin,ikmax          ! range for ik
    real                                    :: ikreal
    integer,allocatable                     :: nptotsqr(:)          ! number of modes per iksqr=ik**2
    integer,allocatable                     :: nptot(:)             ! number of modes per ik
    real*4                                  :: j0
    real                                    :: r                    ! correlation length (i.e. nearest-neighbor spacing)
    real                                    :: kr                   ! = |k|*|r|
    real*4,allocatable                      :: kptot(:)
    real,parameter                          :: pi = 3.14159265
    integer                                 :: normalization
    logical                                 :: verbose
    
    ! input check
    if ((size(delta_r,1).ne.size(delta_r,2)).or.(size(delta_r,1).ne.size(delta_r,3))) then
        print*,'ERROR: delta_r must have dimension (n,n,n).'
        stop
    end if
    
    ! handle options
    normalization = 3
    if (present(opt_normalization)) then
        if ((opt_normalization>=1).and.(opt_normalization<=3)) then
            normalization = opt_normalization
        end if
    end if
    if (present(opt_verbose)) then
        verbose = opt_verbose
    else
        verbose = .true.
    end if
    
    if (verbose) call tic('COMPUTE POWER SPECTRUM p(k) AND/OR TWO-POINT CORRELATION xi(r)')
    
    ! initialization
    n = size(delta_r,1)
    dk = 2*pi/L
    ikmax = (n-1)/2
    ikmin = -n/2
    
    ! make DFT
    call make_or_load_DFT(delta_r,delta_k,opt_verbose=verbose)
    
    ! calculate total power for every square mode |k|^2
    if (verbose) call tic('numerical integration')
    iksqrmax = 3*ikmin**2
    allocate(ptotsqr(0:iksqrmax),nptotsqr(0:iksqrmax))
    ptotsqr = 0
    nptotsqr = 0
    do ik3 = ikmin,ikmax
        ik3sqr = ik3**2
        do ik2 = ikmin,ikmax
            ik23sqr = ik2**2+ik3sqr
            do ik1 = ikmin,ikmax
                iksqr = ik1**2+ik23sqr
                ptotsqr(iksqr) = ptotsqr(iksqr)+cabs(delta_k(ik1,ik2,ik3))**2
                nptotsqr(iksqr) = nptotsqr(iksqr)+1
            end do
        end do
    end do
    deallocate(delta_k)
    
    ! smooth to one total power value one value per linear mode |k|
    allocate(ptot(0:ikmax),nptot(0:ikmax),kptot(0:ikmax))
    ptot = 0
    nptot = 0
    kptot = 0
    do iksqr = 1,floor((ikmax+0.4999)**2)
        ik = nint(sqrt(real(iksqr)))
        ptot(ik) = ptot(ik)+ptotsqr(iksqr)
        nptot(ik) = nptot(ik)+nptotsqr(iksqr)
        ikreal = (real(iksqr-1)**1.5+real(iksqr+1)**1.5)/real(2.0*iksqr) ! k^2-weighted mean of k
        kptot(ik) = kptot(ik)+ikreal*nptotsqr(iksqr)
    end do
    
    ! compute power spectra
    allocate(scalek(1:ikmax),p(1:ikmax))
    do ik = 1,ikmax
        if (nptot(ik)>0) then
            p(ik) = ptot(ik)/real(nptot(ik)) ! [-]
            scalek(ik) = kptot(ik)/real(nptot(ik))*dk ! [h/Mpc]
        else
            p(ik) = 0
            scalek(ik) = ik*dk ! [h/Mpc]
        end if
    end do

    ! make list of scale lengths
    call make_list_of_scale_lengths(n,L,scale_min,scale)
    m = size(scale,1)

    ! integration of 2-point correlation
    allocate(xi2(m))
    xi2 = 0
    do i = 1,m
        r = scale(i)
         do iksqr = 1,iksqrmax
            kr = sqrt(real(iksqr))*dk*r
            j0 = sin(kr+1e-20)/(kr+1e-20)
            xi2(i) = xi2(i)+j0*ptotsqr(iksqr)
         end do
    end do
    
    ! NB: the computation of xi2 above is mathematically *identical* (within hardly measurable numerical
    ! errors) to the following raw computation:
    !   real :: dkr
    !   do i = 1,m
    !       r = scale(i)
    !       dkr = dk*r
    !       do ik3 = ikmin,ikmax
    !           ik3sqr = ik3**2
    !           do ik2 = ikmin,ikmax
    !               ik23sqr = ik2**2+ik3sqr
    !               do ik1 = ikmin,ikmax
    !                   iksqr = ik1**2+ik23sqr
    !                   kr = sqrt(real(iksqr))*dkr
    !                   j0 = sin(kr+1e-20)/(kr+1e-20)
    !                   xi2(i) = xi2(i)+j0*cabs(delta_k(ik1,ik2,ik3))**2
    !               end do
    !           end do
    !       end do
    !   end do
    
    ! change normalization of power spectrum
    if (normalization==2) p = p*dk**dim   ! Wolstenhulme
    if (normalization==3) p = L**dim*p    ! CAMB convention
    
    if (verbose) call toc
    if (verbose) call toc

end subroutine compute_xi2

subroutine compute_xi3(delta_r,L,scale,xi3,scale_min,opt_accuracy)

    ! BRIEF EXPLANATION
    ! Computes the isotropic 3-point correlation function of the density field delta_r, defined on a
    ! cubic (square) box of side-length L and volume V=L^3 (L^2)
    !
    ! INPUT VARIABLES
    ! delta_r(:,:,:)        density perturbation field of dimension (n,n,1) in 2D or (n,n,n) in 3D
    !                       lower and upper bounds are irrelevant
    ! L                     side-length of box
    ! opt_accuracy          optional argument for functype=1, forcing an increased accuracy of computation (if>1)
    !
    ! OUTPUT VARIABLES
    ! scale(:)              correlation scale r
    ! xi2(:)                isotropic 3-point correlation function xi_3(r) for 3 points on a line
    !
    ! NB: This function is not fully optimized and stable. Not forseen to be released publically

    implicit none
    
    ! variable declaration
    real,intent(in)                 :: delta_r(:,:,:)   ! Density field
    real,intent(in)                 :: L                ! side-length of density field delta_r
    real,allocatable,intent(out)    :: scale(:)         ! list of scale lengths of correlation function
    real,allocatable,intent(out)    :: xi3(:)           ! values of correlation function
    real,intent(in)                 :: scale_min         ! minimum scale on which to evaluate the function
    real,optional,intent(in)        :: opt_accuracy     ! precision factor
    complex,allocatable             :: delta_k(:,:,:)   ! Phase-factors of Fourier modes of density field delta_r
    integer                         :: n                ! n^dim = number of elements in delta_k
    integer                         :: m                ! number of elements in scale(:)
    real                            :: dk               ! linear spacing of Fourier modes
    integer(kind=8)                 :: i                ! general counter
    integer                         :: ik,iq,itheta     ! index of vectors k and q and the angle spanned by them
    real                            :: k1,k2,k3         ! cartesian coordinates of vector k
    real                            :: q1,q2,q3         ! cartesian coordinates of vector q
    integer                         :: ik1,ik2,ik3      ! indices for the cartesian coordinates of the vector k
    integer                         :: iq1,iq2,iq3      ! indices for the cartesian coortinates of the vector q
    real*4,allocatable              :: b(:,:,:)         ! isotropic bi-spectrum
    real*4,allocatable              :: w(:,:,:)         ! weights for the integration of the bi-spectrum
    real                            :: utheta,uphi      ! spherical coordinates of random unit vector
    real                            :: ux,uy,uz         ! cartesian coordinates of random unit vector
    real                            :: alpha,sina,cosa  ! random rotation angle with sine and cosine
    real                            :: theta,cost,sint  ! angle between vectors k and q with sine and cosine
    real                            :: r                ! correlation length (i.e. nearest-neighbor spacing)
    real                            :: kqr              ! = |k-q|*|r|
    real                            :: rot(3,2)         ! some elements of a random rotation matrix
    integer(kind=8)                 :: niterations      ! nb of random rotations used to calculate the bi-spectrum
    integer                         :: ntheta           ! nb of angles angle(k,q) considered in the bi-spectrum
    integer                         :: ikmax            ! range of ik
    real                            :: accuracy
    real,parameter                  :: pi = 3.14159265

    call tic('COMPUTE CORRELATION FUNCTION xi_3(r)')
    
    ! input check
   if ((size(delta_r,1).ne.size(delta_r,2)).or.(size(delta_r,1).ne.size(delta_r,3))) then
      print*,'ERROR: delta_r must have dimension (n,n,n).'
      write(*,*) size(delta_r,1),size(delta_r,2),size(delta_r,3)
      stop
   end if
    
    ! handle options
    if (present(opt_accuracy)) then
        accuracy = opt_accuracy
    else
        accuracy = 1
    end if
    
    ! initialization
    niterations = nint(100*accuracy)
    ntheta        = nint(10*max(1.0,min(9.0,sqrt(accuracy))))
    write(*,*) ntheta,accuracy
    n = size(delta_r,1)
    dk = 2*pi/L
    ikmax = (n-1)/2
    
    ! make list of scale lengths
    call make_list_of_scale_lengths(n,L,max(L/real(n)*2.0,scale_min),scale)
    m = size(scale,1)
    allocate(xi3(m))
    
    ! make DFT
    call make_or_load_DFT(delta_r,delta_k)

    ! make isotropic bispectrum
    call tic('make isotropic bispectrum')
    call srand(1)
    allocate(b(0:ikmax/2,0:ikmax/2,0:ntheta-1))
    b = 0
    do i = 1,niterations
    
        ! make random rotation matrix
        utheta = acos(rand()*2.0-1.0)
        uphi = rand()*2*pi
        ux = sin(utheta)*cos(uphi)
        uy = sin(utheta)*sin(uphi)
        uz = cos(utheta)
        alpha = rand()*2*pi
        cosa = cos(alpha)
        sina = sin(alpha)
        rot(1,1) = cosa+ux**2*(1-cosa)
        rot(1,2) = ux*uy*(1-cosa)-uz*sina
        rot(2,1) = ux*uy*(1-cosa)+uz*sina
        rot(2,2) = cosa+uy**2*(1-cosa)
        rot(3,1) = ux*uz*(1-cosa)-uy*sina
        rot(3,2) = uy*uz*(1-cosa)+ux*sina
            
        k1 = rot(1,1)
        k2 = rot(2,1)
        k3 = rot(3,1)
        
        ! apply random rotation matrix to pairs of vectors {k,q}
        do itheta = 0,ntheta-1
            theta = (itheta+0.5)/real(ntheta)*pi
            cost = cos(theta)
            sint = sin(theta)
            q1 = rot(1,1)*cost+rot(1,2)*sint
            q2 = rot(2,1)*cost+rot(2,2)*sint
            q3 = rot(3,1)*cost+rot(3,2)*sint
            do ik = 0,ikmax/2
                ik1 = nint(k1*ik)
                ik2 = nint(k2*ik)
                ik3 = nint(k3*ik)
                do iq = 0,ikmax/2
                    iq1 = nint(q1*iq)
                    iq2 = nint(q2*iq)
                    iq3 = nint(q3*iq)
                    b(ik,iq,itheta) = b(ik,iq,itheta)+realpart(delta_k(ik1,ik2,ik3)*delta_k(iq1,iq2,iq3) &
                    & *delta_k(-ik1-iq1,-ik2-iq2,-ik3-iq3))
                end do
            end do
        end do
    end do
    b = b/niterations
    
    ! fast integration
    call tic('integrate isotropic bispectrum')
    allocate(w(0:ikmax/2,0:ikmax/2,0:ntheta-1))
    do i = 1,m
        r = scale(i)
        do ik = 0,ikmax/2
            do iq = 0,ikmax/2
                do itheta = 0,ntheta-1
                    theta = (itheta+0.5)/real(ntheta)*pi
                    kqr = sqrt(real(ik)**2+real(iq)**2-2*real(ik)*real(iq)*cos(theta)+1e-10)*dk*r
                    w(ik,iq,itheta) = sin(kqr+1e-20)/(kqr+1e-20)*sin(theta)/2
                end do
                w(ik,iq,:) = w(ik,iq,:)*(4*pi)**2*(ik*iq+0.0)**2 ! multiplication by S_dim
            end do
        end do
        w = w*pi/real(ntheta) ! multiplication by dtheta
        xi3(i) = sum(w*b)
    end do
    
    call toc
    call toc

end subroutine compute_xi3

subroutine compute_bispectrum(delta_r,L,bispectrum,angle_list,kleg_list,opt_accuracy,opt_normalization)

    ! BRIEF EXPLANATION
    ! Computes the equvilateral isotropic bispectrum B(k,q) of the 3D density field delta_r, defined on a
    ! cubic box of side-length L and volume V=L^3. The algorithm is uses a Monte Carlo
    ! integration that stochastically samples the triangles in 3D Fourier space. 
    !
    ! INPUT VARIABLES
    ! delta_r(:,:,:)        3D density perturbation field of dimension (n,n,n)
    !                       lower and upper bounds are irrelevant
    ! L                     side-length of box
    ! opt_accuracy          (default=100) this optional argument can be used to decrease (<100) or increase (>100) the
    !                       default number (100) of Monte Carlo iterations
    ! opt_normalization     [1,2,3] normalization type of power spectrum (1: Obreschkow, 2: Wolstenhulme, 3: CAMB)
    !
    ! OUTPUT VARIABLES
    ! b(:,:)                bispectrum
    ! angle_list            [deg]
    ! kleg_list             [units of 1/length]

    implicit none
    
    ! variable declaration
    real,parameter                  :: dangle_deg = 10
    integer,parameter               :: output_convention=3      ! 1: Obreschkow, 2: Wolstenhulme, 3: CAMB
    real,intent(in)                 :: delta_r(:,:,:)           ! Density field
    real,intent(in)                 :: L                        ! side-length of density field delta_r
    real,allocatable,intent(out)    :: bispectrum(:,:)          ! equilateral isotropic bispectrum
    real,allocatable,intent(out)    :: angle_list(:)            ! angle between vectors
    real,allocatable,intent(out)    :: kleg_list(:)             ! length of two equally long vectors
    real,optional,intent(in)        :: opt_accuracy             ! precision factor
    integer,intent(in),optional     :: opt_normalization        ! normalization type of power spectrum
    complex,allocatable             :: delta_k(:,:,:)           ! Phase-factors of Fourier modes of density field delta_r
    integer,parameter               :: dim = 3                  ! number of dimensions
    integer                         :: n                        ! n^dim = number of elements in delta_k
    real                            :: dk                       ! linear spacing of Fourier modes
    integer(kind=8)                 :: i                        ! general counter
    real*4,allocatable              :: rnd(:,:),k(:,:),q(:,:)
    real*8,allocatable              :: w(:,:,:)                 ! counter (needs to be kind=8 in case of high number of interations)
    real*8,allocatable,target       :: b(:,:,:)                 ! bispectrum (needs to be kind=8 in case of high number of interations)
    real*8,pointer                  :: bp(:)
    integer(kind=8)                 :: n_iterations
    integer*4                       :: ikmax,ikleg,ianglemax
    integer                         :: super_iteration
    integer                         :: n_super_iterations
    real                            :: accuracy
    real,parameter                  :: pi = 3.14159265
    integer*4                       :: ik(3),iq(3),ip(3)
    real*4                          :: dangle,angle,base_length
    integer                         :: normalization
    integer                         :: rnd_size
    integer*4,allocatable           :: seed_array(:),iangle(:),ikleg_max(:)
    integer                         :: tstop
    real                            :: dt,remaining,progress
    
    call tic('COMPUTE BISPECTRUM [Stochastic algorithm]')
    
    ! input check
    if ((size(delta_r,1).ne.size(delta_r,2)).or.(size(delta_r,1).ne.size(delta_r,3))) then
        print*,'ERROR: delta_r must have dimension (n,n,n).'
        stop
    end if
    
    ! handle options
    if (present(opt_accuracy)) then
        if (opt_accuracy<=0) then
            print*,'ERROR: opt_accuracy must be larger than 0.'
            stop
        else
            accuracy = opt_accuracy
        end if
        write(*,'(A,F8.1,A)') '+ Accuracy of Monte Carlo sampling = ',accuracy,' (set by user)'
    else
        accuracy = 1
        write(*,'(A,F8.1,A)') '+ Accuracy of Monte Carlo sampling = ',accuracy,' (set by automatically)'
    end if
    normalization = 3
    if (present(opt_normalization)) then
        if ((opt_normalization>=1).and.(opt_normalization<=3)) then
            normalization = opt_normalization
        end if
    end if
    
    ! initialization
    n_iterations = 100000
    n_super_iterations = nint(50*accuracy)
    n = size(delta_r,1)
    dk = 2*pi/L
    ikmax = (n-1)/2
    
    ! define angular discretization of bispectrum
    dangle = dangle_deg/180.0*pi
    ianglemax = nint(pi/dangle)
    allocate(angle_list(0:ianglemax))
    angle_list = real((/(i,i=0,ianglemax)/))*dangle*180.0/pi
    allocate(kleg_list(1:ikmax))
    kleg_list = real((/(i,i=1,ikmax)/))*dk
    
    ! make DFT
    call make_or_load_DFT(delta_r,delta_k)
    
    ! memory allocations
    allocate(b(1:ikmax,0:ianglemax,n_super_iterations))
    allocate(w(1:ikmax,0:ianglemax,n_super_iterations))
    b = 0
    w = 0
    call random_seed(size=rnd_size)
    allocate(seed_array(rnd_size))
    seed_array = 1
    allocate(rnd(n_iterations,2))
    
    ! super iteration
    call tic('numerical integration')
    
    !$OMP PARALLEL PRIVATE(rnd,k,q,angle,base_length,i,ikleg_max,ikleg,ik,iq,ip,bp,iangle)
    !$OMP DO SCHEDULE(DYNAMIC)

    do super_iteration = 1,n_super_iterations
    
        ! print progress report
        if (show_progress) then
            progress = real(super_iteration)/n_super_iterations
            call system_clock(tstop)
            dt = real(tstop-tstart(tindex))/trate
            remaining = max(0.0,dt/progress*(1-progress))
            if (remaining<3600) then
                write(*,'(A,f6.2,A,f6.1,A)',advance='no') char(13)//'+ numerical integration (progress: ', &
                & progress*100.0,'%, remaining time: ',remaining,'s)'
            else
                write(*,'(A,f6.2,A,f6.1,A)',advance='no') char(13)//'+ numerical integration (progress: ', &
                & progress*100.0,'%, remaining time: ',remaining/3600,'h)'
            end if
        end if
                         
        ! chose random unit vectors (k,q) for fast integration of bispectrum
        !$OMP CRITICAL
        call random_seed(put=seed_array*super_iteration+1000)
        call random_number(rnd)
        !$OMP END CRITICAL
        rnd(:,1) = rnd(:,1)*2.0*pi     ! longitude
        rnd(:,2) = rnd(:,2)*2.0-1.0    ! cos(latitude)
        allocate(k(n_iterations,dim))
        k(:,1) = cos(rnd(:,1))*sqrt(1-rnd(:,2)**2)
        k(:,2) = sin(rnd(:,1))*sqrt(1-rnd(:,2)**2)
        k(:,3) = rnd(:,2)
        
        allocate(q(n_iterations,dim))
        q = cshift(k,1)
        
        ! calculate angles between k and q
        allocate(iangle(n_iterations))
        allocate(ikleg_max(n_iterations))
        do i = 1,n_iterations
            angle = acos(sum(k(i,:)*q(i,:)))
            iangle(i) = nint(angle/dangle)
            base_length = 2*sin((pi-angle)/2.0)
            ikleg_max(i) = min(ikmax,floor(ikmax/max(base_length,1e-10)-0.7)) ! 0.7 not really necessary, but cuts 1-2 of the uncertain high-k values
        end do
        
        ! fast summation
        do i = 1,n_iterations
            bp => b(1:ikleg_max(i),iangle(i),super_iteration)
            do ikleg = 1,ikleg_max(i)
                ik = nint(ikleg*k(i,:))
                iq = nint(ikleg*q(i,:))
                ip = -ik-iq
                !if (maxval(abs(ip))>ikmax) stop
                bp(ikleg) = bp(ikleg)+&
                    &real(delta_k(ik(1),ik(2),ik(3))*delta_k(iq(1),iq(2),iq(3))*delta_k(ip(1),ip(2),ip(3)))
            end do
            w(1:ikleg_max(i),iangle(i),super_iteration) = w(1:ikleg_max(i),iangle(i),super_iteration)+1
        end do
        
        deallocate(k,q,iangle,ikleg_max)
        
    end do
    
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL
    
    if (show_progress) then
        write(*,'(A,f6.2,A)',advance='no') char(13)//'+ numerical integration (progress: ',100.0,'%, remaining time: 0.0)      '
    end if
    
    ! finalize bispectrum
    b(:,:,1) = sum(b,3)
    w(:,:,1) = sum(w,3)
    allocate(bispectrum(1:ikmax,0:ianglemax))
    where(w(:,:,1)>0) bispectrum = real(b(:,:,1))/real(w(:,:,1))
    deallocate(b,w)
    if (normalization==2) bispectrum = bispectrum*dk**dim    ! Wolstenhulme convention
    if (normalization==3) bispectrum = L**(3*dim)/(2*pi)**(3*dim/2)*(bispectrum*dk**dim) ! CAMB convention (Bcamb = V^3 / (2pi)^9/2 Bwolst)

    call toc
    call toc

end subroutine compute_bispectrum

subroutine compute_lic(delta_r,L,scale,value,error,scale_min,opt_accuracy,opt_edgeworth)

    ! BRIEF EXPLANATION
    ! Computes the isotropic line-correlation function l(r) of the 3D density field delta_r, defined on a
    ! cubic box of side-length L and volume V=L^3. The algorithm is uses a Monte Carlo
    ! integration that stochastically samples the phase-triangles in 3D Fourier space. 
    ! This function uses the new definition of l(r) with |k|,|q|,|k+q|<2*pi/r instead of |k|,|q|<pi/r.
    !
    ! INPUT VARIABLES
    ! delta_r(:,:,:)    density perturbation field of dimension (n,n,n) in 3D
    !                   lower and upper bounds are irrelevant
    ! L                 side-length of box
    ! opt_accuracy      (default=100) this optional argument can be used to decrease (<100) or increase (>100) the
    !                   default number (100) of Monte Carlo iterations
    ! opt_edgeworth     optional argument; if set to .true., the isotropic line-correlation function l(r), is
    !                   computed using the so-called Edgeworth exapansion without (N>3)-point correlations in
    !                   the loop corrections
    !
    ! OUTPUT VARIABLES
    ! scale(:)          [unit of L] correlation scale r
    ! value(:)          value of correlation function l(r)
    ! error(:)          statistical uncertainty of correlation function l(r)

    implicit none
    
    ! variable declaration
    real,intent(in)                 :: delta_r(:,:,:)       ! Density field
    real,intent(in)                 :: L                    ! side-length of density field delta_r
    real,allocatable,intent(out)    :: scale(:)             ! list of scale lengths of correlation function
    real,allocatable,intent(out)    :: value(:)             ! values of correlation function
    real,allocatable,intent(out)    :: error(:)             ! statistical uncertainties of correlation function
    real,intent(in)                 :: scale_min             ! minimum scale
    real,optional,intent(in)        :: opt_accuracy         ! precision factor
    logical,optional,intent(in)     :: opt_edgeworth        ! decides whether l(r) is calculated in the Edgeworth approximation
    complex,allocatable             :: delta_k(:,:,:)       ! Fourier modes of density field delta_r
    integer,parameter               :: dim = 3              ! number of dimensions
    integer                         :: n                    ! n^dim = number of elements in delta_k
    integer                         :: m                    ! number of elements in scale
    real                            :: dk                   ! linear spacing of Fourier modes
    integer(kind=8)                 :: i                    ! general counter
    complex,allocatable             :: epsil_k(:,:,:)       ! Phase-factors of Fourier modes of density field delta_r
    real*4,allocatable              :: rnd(:,:)
    real*4,allocatable              :: k(:,:)
    real*4,allocatable              :: q(:,:)
    real*4,allocatable              :: w(:)
    real*4,allocatable              :: x(:)
    integer                         :: ir
    real                            :: r                    ! correlation length (i.e. nearest-neighbor spacing)
    real(kind=4)                    :: n_modes              ! number of modes
    integer(kind=8)                 :: n_iterations
    real                            :: ikmax
    integer(kind=8)                 :: n_iterations_max
    integer                         :: super_iteration
    integer                         :: n_super_iterations
    real,allocatable                :: svalue(:,:)
    real                            :: accuracy
    real                            :: standard_deviation
    real,parameter                  :: pi = 3.14159265
    integer                         :: tstop
    real                            :: dt,remaining,progress
    logical                         :: edgeworth
    real,allocatable                :: scaler(:)            ! [unit of L] list of scale lengths of correlation function
    real,allocatable                :: xi2(:)               ! values of correlation function
    real,allocatable                :: scalek(:)            ! [1/unit of L] wave vector
    real,allocatable                :: p(:)
    real*4,allocatable              :: sqrtp(:)
    real*4                          :: kvalue,f
    integer                         :: iktop,index
    
    ! handle options
    if (present(opt_edgeworth)) then
        edgeworth = opt_edgeworth
    else
        edgeworth = .false.
    end if
    
    if (edgeworth) then
        call tic('COMPUTE LINE CORRELATION l(r) IN EDGEWORTH APPROXIMATION [Stochastic algorithm]')
    else
        call tic('COMPUTE LINE CORRELATION FUNCTION l(r) [Stochastic algorithm]')
    end if
    
    ! input check
    if ((size(delta_r,1).ne.size(delta_r,2)).or.(size(delta_r,1).ne.size(delta_r,3))) then
        print*,'ERROR: delta_r must have dimension (n,n,n).'
        stop
    end if
    
    ! handle options
    if (present(opt_accuracy)) then
        if (opt_accuracy<=0) then
            print*,'ERROR: opt_accuracy must be larger than 0.'
            stop
        else
            accuracy = opt_accuracy
        end if
        write(*,'(A,F8.1,A)') '+ Accuracy of Monte Carlo sampling = ',accuracy,' (set by user)'
    else
        accuracy = 1
        write(*,'(A,F8.1,A)') '+ Accuracy of Monte Carlo sampling = ',accuracy,' (set by automatically)'
    end if
    
    ! initialization
    n_iterations_max = 100000
    n_super_iterations = nint(500*accuracy)
    n = size(delta_r,1)
    dk = 2*pi/L
    
    ! make list of scale lengths
    call make_list_of_scale_lengths(n,L,max(L/real(n)*2.0,scale_min),scale)
    m = size(scale,1)
    allocate(svalue(n_super_iterations,m),value(m),error(m))
        
    ! make DFT
    call make_or_load_DFT(delta_r,delta_k)
    
    ! if needed, evaluate power spectrum and make lookup table
    if (edgeworth) then
        call tic('compute powerspectrum')
        iktop = (n-1)/2
        call compute_xi2(delta_r,L,scaler,xi2,scalek,p,0.0,opt_normalization=1,opt_verbose=.false.)
        allocate(sqrtp(0:(iktop+2)**2))
        index = 1
        do i = 0,ubound(sqrtp,1)
            kvalue = sqrt(real(i))*dk
            if ((kvalue>scalek(index)).and.(index<ubound(p,1))) index=index+1
            if (index>1) then
                f = (scalek(index)-kvalue)/(scalek(index)-scalek(index-1))
                sqrtp(i) = sqrt(p(index-1)*f+p(index)*(1-f))
            else
                sqrtp(i) = sqrt(p(index))
            end if
            if (index>1) then
                f = (scalek(index)-kvalue)/(scalek(index)-scalek(index-1))
                sqrtp(i) = sqrt(p(index-1)*f+p(index)*(1-f))
            else
                sqrtp(i) = sqrt(p(index))
            end if
        end do
    end if

    ! make phase factors
    if (.not.edgeworth) then
        call tic('compute phase factors')
        call compute_phase_factors(delta_k,epsil_k)
        deallocate(delta_k)
    end if

    ! super iteration
    call tic('numerical integration')
    
    !$OMP PARALLEL PRIVATE(r,ikmax,x,n_modes,n_iterations,i,rnd,k,q,w,ir)
    !$OMP DO SCHEDULE(DYNAMIC)
    
    do super_iteration = 1,n_super_iterations
    
        ! print progress report
        if (show_progress) then
            progress = real(super_iteration)/n_super_iterations
            call system_clock(tstop)
            dt = real(tstop-tstart(tindex))/trate
            remaining = max(0.0,dt/progress*(1-progress))
            if (remaining<3600) then
                write(*,'(A,f6.2,A,f6.1,A)',advance='no') char(13)//'+ numerical integration (progress: ', &
                & progress*100.0,'%, remaining time: ',remaining,'s)'
            else
                write(*,'(A,f6.2,A,f6.1,A)',advance='no') char(13)//'+ numerical integration (progress: ', &
                & progress*100.0,'%, remaining time: ',remaining/3600,'h)'
            end if
        end if
        
        ! make array of mode-couples (k,q) for fast integration of normalized bispectrum
        allocate(rnd(n_iterations_max,dim))
        call random_number(rnd)
        rnd(:,1) = rnd(:,1)**(1.0/real(dim))     ! radial coordinate
        rnd(:,2) = rnd(:,2)*2.0*pi                  ! longitude
        rnd(:,3) = rnd(:,3)*2.0-1.0                ! cos(latitude)
        allocate(k(n_iterations_max,dim))
        k(:,1) = rnd(:,1)*cos(rnd(:,2))*sqrt(1-rnd(:,3)**2)
        k(:,2) = rnd(:,1)*sin(rnd(:,2))*sqrt(1-rnd(:,3)**2)
        k(:,3) = rnd(:,1)*rnd(:,3)
        deallocate(rnd)
        allocate(q(n_iterations_max,dim))
        q = cshift(k,1)
        
        ! calculation of weighting function w_D
        allocate(w(n_iterations_max))
        where(sum((k+q)**2,2)<=1)
            w = sqrt(sum((k-q)**2,2))*pi*2
        elsewhere
            w = 0
        end where
        where(w.ne.0) w = sin(w+1e-20)/(w+1e-20)
        
        ! fast integration of normalized bispectrum
        do ir = 1,m
    
            r = scale(ir)
            ikmax = min(real((n-1)/2),2*pi/r/dk)
            n_modes = 4.0*pi/3.0*ikmax**3
            n_iterations = nint(min(real(n_iterations_max),n_modes**2))
            
            allocate(x(n_iterations))
            x = 0
            
            if (edgeworth) then
                forall (i=1:n_iterations,w(i).ne.0)
                    x(i) = w(i)*real(delta_k(nint(k(i,1)*ikmax),nint(k(i,2)*ikmax),nint(k(i,3)*ikmax))* &
                    & delta_k(nint(q(i,1)*ikmax),nint(q(i,2)*ikmax),nint(q(i,3)*ikmax))* &
                    & delta_k(-nint(k(i,1)*ikmax)-nint(q(i,1)*ikmax), &
                    & -nint(k(i,2)*ikmax)-nint(q(i,2)*ikmax),-nint(k(i,3)*ikmax)-nint(q(i,3)*ikmax)))/ &
                    & sqrtp(sum(nint(k(i,:)*ikmax)**2))/ &
                    & sqrtp(sum(nint(q(i,:)*ikmax)**2))/ &
                    & sqrtp(sum((nint(k(i,:)*ikmax)+nint(q(i,:)*ikmax))**2))
                end forall
            else
                forall (i=1:n_iterations,w(i).ne.0)
                    x(i) = w(i)*real(epsil_k(nint(k(i,1)*ikmax),nint(k(i,2)*ikmax),nint(k(i,3)*ikmax))* &
                    & epsil_k(nint(q(i,1)*ikmax),nint(q(i,2)*ikmax),nint(q(i,3)*ikmax))* &
                    & epsil_k(-nint(k(i,1)*ikmax)-nint(q(i,1)*ikmax), &
                    & -nint(k(i,2)*ikmax)-nint(q(i,2)*ikmax),-nint(k(i,3)*ikmax)-nint(q(i,3)*ikmax)))
                end forall
            end if
            
            svalue(super_iteration,ir) = sum(x)/real(n_iterations)*n_modes**2*(r/L)**(3.0*dim/2.0)
            deallocate(x)
            
        end do
        
        deallocate(k,q,w)
        
    end do
    
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL
    
    if (show_progress) then
        write(*,'(A,f6.2,A)',advance='no') char(13)//'+ numerical integration (progress: ',100.0,'%, remaining time: 0.0)      '
    end if
    
    ! pre-factor (sqrt(pi)/2)^3 of Edgeworth approx
    if (edgeworth) then
        svalue = (sqrt(pi)/2.0)**3*svalue
    end if
    
    ! calculate average correlation values
    do ir = 1,m
        value(ir) = sum(svalue(:,ir))/real(n_super_iterations)
    end do
    
    ! standard error of the calculation (vanishes as number of iteration becomes infinite)
    do ir = 1,m
        standard_deviation = sqrt(sum((svalue(:,ir)-value(ir))**2)/real(n_super_iterations))
        error(ir) = standard_deviation/sqrt(real(n_super_iterations))
    end do
    
    call toc
    call toc

end subroutine compute_lic

subroutine compute_lic_deterministic(delta_r,L,scale,value,scale_min,opt_accuracy,opt_performance)

    ! BRIEF EXPLANATION
    ! Computes the isotropic line correlation function l(r) of the 3D density field delta_r, defined on a
    ! cubic box of side-length L and volume V=L^3. This function uses a deterministic algorithm that samples
    ! *ALL* triagnles in Fourier space (unless opt_completeness is set <1), thus producing the exact line correlation
    ! of the given density field. This function uses the new definition of l(r) with |k|,|q|,|k+q|<2*pi/r instead
    ! of |k|,|q|<pi/r.
    !
    ! INPUT VARIABLES
    ! delta_r(:,:,:)    density perturbation field of dimension (n,n,n) in 3D
    !                   lower and upper bounds are irrelevant
    ! L                 side-length of box
    ! opt_accuracy      optional argument completeness of Fourier space sampling 0..1
    !
    ! OUTPUT VARIABLES
    ! scale(:)          [unit of L] correlation scale r
    ! value(:)          value of correlation function l(r)
    ! opt_performance   vector: 1: total nb of terms in the sum_{k,q}, 2: terms/sec

    implicit none
    
    ! variable declaration
    real,intent(in)                 :: delta_r(:,:,:)       ! Density field
    real,intent(in)                 :: L                    ! side-length of density field delta_r
    real,allocatable,intent(out)    :: scale(:)             ! list of scale lengths of correlation function
    real,allocatable,intent(out)    :: value(:)             ! values of correlation function
    real,intent(in)                 :: scale_min             ! minimum scale
    real,optional,intent(in)        :: opt_accuracy         ! completeness of Fourier space sampling 0..1
    integer*8,optional,intent(out)  :: opt_performance(2)   ! optional speed performance measurement
    complex,allocatable             :: delta_k(:,:,:)       ! Fourier modes of density field delta_r
    integer,parameter               :: dim = 3              ! number of dimensions
    integer                         :: n                    ! n^dim = number of elements in delta_k
    integer                         :: m                    ! number of elements in scale
    real                            :: dk                   ! linear spacing of Fourier modes
    complex,allocatable,target      :: epsil_k(:,:,:)       ! Phase-factors of Fourier modes of density field delta_r
    integer                         :: ir
    integer*4                       :: ikmin,ikmax
    integer*4                       :: ik1,ik2,ik3
    integer*4                       :: iq1,iq2,iq3
    integer*4                       :: ip1,ip2,ip3
    integer*4                       :: sk
    integer*4                       :: sq1,sp1
    integer*4                       :: sp12,sq12
    integer*4                       :: sdkq1,sdqp1,sdpk1
    integer*4                       :: sdkq12,sdqp12,sdpk12,sq12msp12msk3
    integer*4                       :: iksqrmax,iksqr,ik1sqr,ik12sqr,sq12msp12,iq3maxx,skmsq12,wmax
    real                            :: b
    real*8,allocatable,target       :: x(:,:)
    real*8,pointer                  :: px(:)
    complex,pointer                 :: pepsil_q_pp1(:),pepsil_q_pn1(:),pepsil_p_pp1(:),pepsil_p_pn1(:)
    complex,pointer                 :: pepsil_q_pp2(:),pepsil_q_pn2(:),pepsil_p_pp2(:),pepsil_p_pn2(:)
    real,parameter                  :: pi = 3.14159265
    integer*4,allocatable           :: ik3max(:,:),ik2max(:)
    real*4                          :: f
    real,allocatable                :: rnd(:,:,:),rsqr(:,:,:)
    integer*4                       :: i,iq3min,iq1max,iq2max,iq3max,twoik3
    real*4,allocatable,target       :: w(:,:)
    real*4,pointer                  :: pw(:,:)
    complex                         :: ep3p2p1,ep3p2m1,ep3m2p1,ep3m2m1,ep3p1p2,em3p1p2,ep3p1m2,em3p1m2
    integer*4,allocatable           :: ir_max(:)
    logical,allocatable             :: selected(:,:,:)
    real*4                          :: completeness,threshold,threshold_max,threshold_min,fraction,ff,weight
    real,parameter                  :: permutations(0:1,0:1) = reshape((/6,3,3,1/)/6.0,(/2,2/))
    real*8,allocatable              :: n_terms(:),n_terms_all(:)
    real*8                          :: n_terms_now
    integer                         :: tstop,i0
    real                            :: progress,dt,remaining
    
    call tic('COMPUTE LINE CORRELATION FUNCTION l(r) [Deterministic algorithm]')
    
    ! input check
    if ((size(delta_r,1).ne.size(delta_r,2)).or.(size(delta_r,1).ne.size(delta_r,3))) then
        print*,'ERROR: delta_r must have dimension (n,n,n).'
        stop
    end if
    
    ! handle options
    if (present(opt_accuracy)) then
        if (opt_accuracy<=0) then
            print*,'ERROR: opt_accuracy in must be larger than 0 and smaller or equal than 1.'
            stop
        else if (opt_accuracy>1) then
            completeness = 1
            write(*,'(A,F8.6,A)') '+ Fourier space completeness (accuracy) = ',completeness,' (set automatically)'
            write(*,'(A)') '  Note: accuracies>1 but are ignored in the determinstic algorithm.'
        else
            completeness = opt_accuracy
            write(*,'(A,F8.6,A)') '+ Fourier space completeness (accuracy) = ',completeness,' (set by user)'
        end if
    else
        completeness = 1
        write(*,'(A,F8.6,A)') '+ Fourier space completeness (accuracy) = ',completeness,' (set automatically)'
    end if
    
    ! initialization
    n = size(delta_r,1)
    dk = 2*pi/L
    ikmax = (n-1)/2
    ikmin = -n/2
        
    ! make DFT
    call make_or_load_DFT(delta_r,delta_k)
    
    ! make phase factors
    call tic('compute phase factors')
    call compute_phase_factors(delta_k,epsil_k)
    deallocate(delta_k)
    
    ! make list of scale lengths
    call make_list_of_scale_lengths(n,L,max(L/real(n)*2.0,scale_min),scale)
    m = size(scale,1)
    if (scale(1)<2*L/real(n)) then
        write(*,*)
        print*,'ERROR: the scale contains radii that are too small for the given discretization.'
        print*,'       the minimal r in l(r) is 2*dr, where dr=L/n is the cell separation.'
        stop
    end if
    
    ! find the maximal value of ir, for any given max(|ik|,|iq|,|ip|)^2
    allocate(ir_max(0:3*ikmax**2))
    do i = 0,3*ikmax**2
        do ir = 1,m
            if ((2.0*pi/scale(ir)/dk)**2.0<=i) exit
        end do
        ir_max(i) = ir-1
    end do
    
    ! find maximum |ik+-iq|^2
    do i = 0,3*ikmax**2
        if (ir_max(i)==0) exit
    end do
    iksqrmax = min(ikmax**2,i-1)

    ! lookup table for Bessel function
    wmax = 4*iksqrmax                                
    allocate(w(1:wmax+1,m))
    do ir = 1,m
        w(1,ir) = 1.0    
        do i = 1,wmax
            f = sqrt(real(i))*dk*scale(ir)
            w(i+1,ir) = sin(f)/f
        end do
    end do

    ! define relevant sphere/half-sphere in 3D Fourier space
    allocate(ik2max(0:ikmax),ik3max(0:ikmax,0:ikmax))
    do ik1 = 0,ikmax
        ik1sqr = ik1**2
        if (iksqrmax-ik1sqr>=0) then
            ik2max(ik1) = min(ik1,floor(sqrt(real(iksqrmax-ik1sqr))))
        else
            ik2max(ik1) = -1
        end if
        do ik2 = 0,ikmax
            if (iksqrmax-ik1sqr-ik2**2>=0) then
                ik3max(ik1,ik2) = floor(sqrt(real(iksqrmax-ik1sqr-ik2**2)))
            else
                ik3max(ik1,ik2) = -1
            end if
        end do
    end do

    ! make selection if incomplete Fourier space sampling is desired
    threshold = 1.0
    allocate(selected(0:ikmax,0:ikmax,0:ikmax))
    if (completeness<1.0) then
        allocate(rnd(0:ikmax,0:ikmax,0:ikmax))
        allocate(rsqr(0:ikmax,0:ikmax,0:ikmax))
        call random_number(rnd)
        do ik1 = 0,ikmax
            ik1sqr = ik1**2
            do ik2 = 0,ikmax
                ik12sqr = ik1sqr+ik2**2
                do ik3 = 0,ikmax
                    iksqr = ik12sqr+ik3**2
                    if (iksqr>iksqrmax) then
                        rsqr(ik1,ik2,ik3) = -1.0
                    else
                        rsqr(ik1,ik2,ik3) = real(iksqrmax)/max(1,iksqr)
                    end if
                end do
            end do
        end do
        threshold_max = 1.0
        threshold_min = 1e-8
        do while ((abs(threshold_max/threshold_min-1)>0.01).and.(threshold_max>0.0))
            threshold = (threshold_max+threshold_min)/2.0
            fraction = count(threshold*rsqr>rnd)/real(pi/6.0*(sqrt(real(iksqrmax))+0.707)**3)
            if (fraction>completeness) then
                threshold_max = threshold
            else
                threshold_min = threshold
            end if
        end do
        threshold = (threshold_max+threshold_min)/2.0
        selected = threshold*rsqr>rnd
    else
        selected = .true.
    end if

    call tic('numerical integration')

    allocate(x(m,0:ikmax),n_terms(0:ikmax),n_terms_all(0:ikmax))
    x = 0
    n_terms = 0
    n_terms_all = 0

    ! EXPLOITED SYMMETRIES
    ! + 6: compute only one of the 6 permultations of the 3 vectors |k+q|<=|q|<=|k|
    ! + 2: coord1 <-> -coord1 same triagle shapes
    ! + 2: coord2 <-> -coord2 same triagle shapes
    ! + 2: coord3 <-> -coord3 same triagle shapes, complex conjugate for epsilon => 2*real part of ik3>=0
    ! + 6: could compute w() for only one of the 6 permultations 3 coordinates, but only the permutation of
    !      coord1 and coord2 can be used to reduce time, because it is a great advantate to only loop over the
    !      first index of epsil, here defined as the 3rd coordinate
    
    !$OMP PARALLEL PRIVATE(ik1,ik2,ik3,twoik3,sk,ir,px,pw,weight, &
    !$OMP ep3p2p1,ep3p2m1,ep3m2p1,ep3m2m1,ep3p1p2,em3p1p2,ep3p1m2,em3p1m2, &
    !$OMP pepsil_q_pp1,pepsil_q_pn1,pepsil_p_pp1,pepsil_p_pn1, &
    !$OMP pepsil_q_pp2,pepsil_q_pn2,pepsil_p_pp2,pepsil_p_pn2, &
    !$OMP iq1max,iq1,ip1,sq1,sp1,sdkq1,sdqp1,sdpk1, &
    !$OMP iq2max,iq2,ip2,sq12,sp12,sq12msp12,sq12msp12msk3,skmsq12, &
    !$OMP iq3min,iq3max,iq3maxx,iq3,i0,sdkq12,sdqp12,sdpk12,ip3,i,b,ff,n_terms_now)
    !$OMP DO SCHEDULE(DYNAMIC)

    do ik1 = 0,ikmax
        do ik2 = 0,ik2max(ik1)
    
        ! print progress report
        if (show_progress) then
            ! NB: the l(r)=sum_{k,q} for r=2*dr (thus kmax=2pi/(2dr)=dk*L/(2dr)=dk*n/2=dk*ikmax) has, a priori,
            !     (4*pi/3*ikmax^3)**2 [~0.274 n^6] terms. However, the truncation condition |k+q|<kmax reduced the number
            !     of couples {k,q} by a factor ~0.469. Thus the total number of terms to evaluate for l(rmin=2dr) is
            !     ~0.469*(4*pi/3*ikmax^3)**2 [~0.1285 n^6]. We here replace ikmax for sqrt(iksqrmax) to get the best estimate
            !     given the actual value of rmin used here.
            progress = min(1.0,max(1.0,real(sum(n_terms)))/(0.469*real(4*pi/3*iksqrmax**1.5)**2.0))
            call system_clock(tstop)
            dt = real(tstop-tstart(tindex))/trate
            remaining = max(0.0,dt/progress*(1-progress))
            if (remaining<3600) then
                write(*,'(A,f6.2,A,f6.1,A)',advance='no') char(13)//'+ numerical integration (progress: ', &
                & progress*100.0,'%, remaining time: ',remaining,'s)'
            else
                write(*,'(A,f6.2,A,f6.1,A)',advance='no') char(13)//'+ numerical integration (progress: ', &
                & progress*100.0,'%, remaining time: ',remaining/3600,'h)'
            end if
        end if

        do ik3 = 0,ik3max(ik1,ik2)
            if (selected(ik1,ik2,ik3)) then

                twoik3 = 2*ik3
                sk = ik1**2+ik2**2+ik3**2
                ir = ir_max(sk) ! index of r=2pi/max(|k|,|q|,|k+q|)
                px => x(1:ir,ik1)
                pw => w(:,1:ir)
                weight = (1+transfer(ik1>0,i))*(1+transfer(ik2>0,i))*(1+transfer(ik3>0,i))
                weight = weight*(1+transfer(ik1>ik2,i))*max(1.0,sk/(threshold*real(iksqrmax)))
                
                ep3p2p1 = weight*epsil_k(+ik3,+ik2,+ik1)
                ep3p2m1 = weight*epsil_k(+ik3,+ik2,-ik1)
                ep3m2p1 = weight*epsil_k(-ik3,+ik2,-ik1) ! conjg
                ep3m2m1 = weight*epsil_k(-ik3,+ik2,+ik1) ! conjg
                ep3p1p2 = weight*epsil_k(+ik3,+ik1,+ik2)
                em3p1p2 = weight*epsil_k(-ik3,+ik1,+ik2)
                ep3p1m2 = weight*epsil_k(+ik3,+ik1,-ik2)
                em3p1m2 = weight*epsil_k(-ik3,+ik1,-ik2)
                    
                iq1max = floor(sqrt(real(sk))) ! condition |q|<=|k|
                
                do iq1 = max(-iq1max,-ikmax-ik1),min(iq1max,ikmax-ik1)
                    ip1 = -ik1-iq1
                    sq1 = iq1**2
                    sp1 = ip1**2
                    sdkq1 = (ik1-iq1)**2
                    sdqp1 = (iq1-ip1)**2
                    sdpk1 = (ip1-ik1)**2
                    iq2max = floor(sqrt(real(sk-sq1))) ! condition |q|<=|k|
                    
                    do iq2 = max(-iq2max,-ikmax-ik2),min(iq2max,ikmax-ik2)
                        ip2 = -ik2-iq2
                        sq12 = sq1+iq2**2
                        sp12 = sp1+ip2**2
                        sq12msp12 = sq12-sp12
                        sq12msp12msk3 = sq12msp12-ik3**2
                        iq3max = floor(sqrt(real(sk-sq12))) ! condition |q|<=|k|
                        
                        if (ik3>0) then
                            iq3maxx = floor(sq12msp12msk3/(2.0*ik3)) ! condition |p|=|k+q|<=|q|
                        else
                            iq3maxx = (transfer(sq12msp12>=0,i)*3-2)*ikmax ! condition |p|=|k+q|<=|q|
                            ! this is an optimized way of ensuring that there is:
                            ! + no condition additional on iq3 (iq3maxx=ikmax), if sq12msp12>=0;
                            ! + no integration over iq3 (iq3maxx=-2ikmax), if sq12msp12<0, which would imply |p|>|q|
                        end if
                        iq3min = max(-iq3max,-ikmax-ik3)
                        iq3max = min(+iq3max,+ikmax-ik3,+iq3maxx)
                        
                        if (iq3max>=iq3min) then
                        
                            skmsq12 = sk-sq12
                            sdkq12 = sdkq1+(ik2-iq2)**2+1
                            sdqp12 = sdqp1+(iq2-ip2)**2+1
                            sdpk12 = sdpk1+(ip2-ik2)**2+1
                        
                            pepsil_q_pp1 => epsil_k(:,+iq2,+iq1)
                            pepsil_q_pn1 => epsil_k(:,+iq2,-iq1)
                            pepsil_p_pp1 => epsil_k(:,+ip2,+ip1)
                            pepsil_p_pn1 => epsil_k(:,+ip2,-ip1)
                            pepsil_q_pp2 => epsil_k(:,+iq1,+iq2)
                            pepsil_q_pn2 => epsil_k(:,+iq1,-iq2)
                            pepsil_p_pp2 => epsil_k(:,+ip1,+ip2)
                            pepsil_p_pn2 => epsil_k(:,+ip1,-ip2)
                            i0 = 1-ikmin ! needed because the custom indexation of pointers, pepsil_p_pn2(ikmin:ikmax) => epsil_k(:,+ip1,-ip2) does not work with all compilers
                            
                            ! NB: it's worth selecting those loops where all ff=6, i.e. ~80% of all loops, via
                            if ((iq3min**2<skmsq12).and.(iq3max**2<skmsq12).and.(twoik3*iq3max<sq12msp12msk3)) then
                                do iq3 = iq3min,iq3max
                                    ip3 = -ik3-iq3
                                    b = real(ep3p2p1*pepsil_q_pp1(i0+iq3)*pepsil_p_pp1(i0+ip3)) &
                                      + real(ep3p2m1*pepsil_q_pn1(i0+iq3)*pepsil_p_pn1(i0+ip3)) &
                                      + real(ep3m2p1*pepsil_q_pn1(i0-iq3)*pepsil_p_pn1(i0-ip3)) & ! conjg
                                      + real(ep3m2m1*pepsil_q_pp1(i0-iq3)*pepsil_p_pp1(i0-ip3)) & ! conjg
                                      + real(ep3p1p2*pepsil_q_pp2(i0+iq3)*pepsil_p_pp2(i0+ip3)) &
                                      + real(em3p1p2*pepsil_q_pp2(i0-iq3)*pepsil_p_pp2(i0-ip3)) &
                                      + real(ep3p1m2*pepsil_q_pn2(i0+iq3)*pepsil_p_pn2(i0+ip3)) &
                                      + real(em3p1m2*pepsil_q_pn2(i0-iq3)*pepsil_p_pn2(i0-ip3))
                                    px = px+b*(pw(sdkq12+(ik3-iq3)**2,:)+pw(sdqp12+(iq3-ip3)**2,:)+pw(sdpk12+(ip3-ik3)**2,:))
                                end do
                            else
                                do iq3 = iq3min,iq3max
                                    ip3 = -ik3-iq3
                                    ff = permutations(transfer(iq3**2==skmsq12,1),transfer(twoik3*iq3==sq12msp12msk3,1))
                                    b = real(ep3p2p1*pepsil_q_pp1(i0+iq3)*pepsil_p_pp1(i0+ip3)) &
                                      + real(ep3p2m1*pepsil_q_pn1(i0+iq3)*pepsil_p_pn1(i0+ip3)) &
                                      + real(ep3m2p1*pepsil_q_pn1(i0-iq3)*pepsil_p_pn1(i0-ip3)) & ! conjg
                                      + real(ep3m2m1*pepsil_q_pp1(i0-iq3)*pepsil_p_pp1(i0-ip3)) & ! conjg
                                      + real(ep3p1p2*pepsil_q_pp2(i0+iq3)*pepsil_p_pp2(i0+ip3)) &
                                      + real(em3p1p2*pepsil_q_pp2(i0-iq3)*pepsil_p_pp2(i0-ip3)) &
                                      + real(ep3p1m2*pepsil_q_pn2(i0+iq3)*pepsil_p_pn2(i0+ip3)) &
                                      + real(em3p1m2*pepsil_q_pn2(i0-iq3)*pepsil_p_pn2(i0-ip3))
                                    px = px+ff*b*(pw(sdkq12+(ik3-iq3)**2,:)+pw(sdqp12+(iq3-ip3)**2,:)+pw(sdpk12+(ip3-ik3)**2,:))
                                end do
                            end if
                            
                            n_terms_now = nint((iq3max-iq3min+1)*6*8*2 &
                             & *(1-0.5*transfer(ik1==0,1))*(1-0.5*transfer(ik2==0,1))*(1-0.5*transfer(ik3==0,1)) &
                             & *(1-0.5*transfer(iq1==0,1))*(1-0.5*transfer(iq2==0,1)) &
                             & *max(1.0,sk/(threshold*real(iksqrmax))))
                             
                            ! 6 permutations of k,q,p=-k-q, which are evaluated by the one particular permutation
                            !   that satisfies |k|>=|q|>=|p|. To apply this condition on |p|, the terms
                            !   w(|k-q|r) had to be symmetried to w(|k-q|r)w(|q-p|r)w(|p-k|r)
                            ! 8 terms im the innner-most loop, which exploit some symmetries of
                            !   w(|k-q|r)w(|q-p|r)w(|p-k|r):
                            !   + coord1<=>-coord1 => ik1>=0
                            !   + coord2<=>-coord2 => ik2>=0
                            !   + and *one* of sevearal possible permutations coord2<=>coord1 => ik2<=ik1
                            !     the other permutations are not exploited, to maintain the summation over the first
                            !     index of epsil_k(:,:,:) in the inner-most loop
                            ! 2 from the symmetry coord3<=>-coord3, encoded in taking the real parts of the 8 terms
                            
                            n_terms(ik1) = n_terms(ik1)+n_terms_now
                            n_terms_all(ik1) = n_terms_all(ik1)+n_terms_now*ir
                            
                        end if
                    end do
                    
                end do
            end if
        end do
        end do
    end do
        
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL
    
    if (show_progress) then
        write(*,'(A,f6.2,A)',advance='no') char(13)//'+ numerical integration (progress: ',100.0,'%, remaining time: 0.0)      '
    end if
    
    allocate(value(m))
    value = real((scale/L)**(3.0*dim/2.0)*sum(x,2)/4) ! factor 4 to balance some symmetries
    
    if (present(opt_performance)) then
        call system_clock(tstop)
        dt = real(tstop-tstart(tindex))/trate
        opt_performance = (/int(sum(n_terms_all),8),int(sum(n_terms_all)/dt,8)/)
    end if

    call toc
    call toc

end subroutine compute_lic_deterministic

subroutine load_density_perturbation_from_file(filename,delta_r)

    ! variable declaration
    implicit none
    character(len=*),intent(in)     :: filename
    real,allocatable,intent(out)    :: delta_r(:,:,:)
    integer*8                       :: file_size
    real*4,allocatable              :: rho4(:,:,:)
    real*8,allocatable              :: rho8(:,:,:)
    integer                         :: n
    character(len=100)              :: string
    integer                         :: bytes
    
    ! initialization
    call tic('LOAD DENSITY PERTURBATION FIELD')
    
    ! determine size of grid
    inquire(file=trim(filename), size=file_size)
    n = nint(real(file_size/4)**(1.0/3.0))
    if (int(n,8)**3*4_8==file_size) then
        bytes = 4
    else
        n = nint(real(file_size/8)**(1.0/3.0))
        if (int(n,8)**3*8_8==file_size) then
            bytes = 8
        else
            write(*,'(A)')
            write(*,'(A)') 'Format of input file not recognized. Consider specifying a different format using -input.'
            stop
        end if
    end if
    write(string,'(I4)') n
    write(*,'(3A)') '  number of cells = ',trim(adjustl(string)),'^3'
    
    ! load file and normalize density field
    allocate(delta_r(0:n-1,0:n-1,0:n-1))
    open(1,file=trim(filename),action='read',form='unformatted',access='stream')
    if (bytes==4) then
        allocate(rho4(n,n,n))
        read(1) rho4
        delta_r = real(rho4)
        deallocate(rho4)
    else if (bytes==8) then
        allocate(rho8(n,n,n))
        read(1) rho8
        delta_r = real(rho8)
        deallocate(rho8)
    end if
    close(1)
    
    call toc

end subroutine load_density_perturbation_from_file

subroutine make_density_perturbation_field_from_particle_file(filename,filetype,species,cellsize,&
    interpolation_method,delta_r,L,subsample_probability,subsample_seed,nparticles)
    
    ! BRIEF EXPLANATION
    ! Computes a 3D density field delta_r from a set of points with equal or different mass in binary, ascii or gadget format.
    ! The density field must be defined on a cubic volume with periodic boundary conditions.

    ! variable declaration
    implicit none
    character(len=*),intent(in)     :: filename
    integer,intent(in)              :: filetype
    integer,intent(in)              :: species
    real,intent(in)                 :: cellsize ! [Mpc/h]
    integer*4,intent(in)            :: interpolation_method ! 1: nearest neighbor, 2: top-hat (has disadvangate of loosing small power, but advantage of better dealing with unsufficiently sampled fields)
    real,allocatable,intent(out)    :: delta_r(:,:,:)
    real,intent(inout)              :: L ! if 0 at input then, L will be determined automatically
    real,intent(in)                 :: subsample_probability
    integer*4,intent(in)            :: subsample_seed
    integer*8,intent(out),optional  :: nparticles ! number of sph particles a side
    real*4,allocatable              :: rho(:,:,:)
    real*4,allocatable              :: x(:),y(:),z(:),mass(:)
    real*4                          :: xempty,yempty,zempty,Lauto
    integer*4                       :: np(6)
    integer*4                       :: npart,npartnew
    integer*8                       :: npart_tot,npartnew_tot
    integer*4                       :: ipart,i
    integer*4                       :: n,m,j,k
    real*4                          :: meanrho
    real*4                          :: h,smoothing,fi,fj,fk
    character(len=100)              :: string,string2
    real*4                          :: fiavg,fjavg,fkavg
    real,parameter                  :: pi = 3.141592653589793
    character(len=255)              :: fn
    integer*4                       :: ifile,nfiles,stat
    complex*16                      :: ci,cj,ck,cf
    real*4,parameter                :: target_offset = 0.25 ! target_offset from cell centre [-0.5..0.5]
    real*4,allocatable              :: rnd(:)
    
    
    call tic('MAKE DENSITY PERTURBATION FIELD FROM PARTICLE FILE')
    
    ! *** initialization 
    call tic('initialize data')
    call get_number_of_files
    
    ! determine total number of particles and box length
    npart_tot = 0
    Lauto = 0
    do ifile = 1,nfiles
        if (nfiles==1) then
            fn = trim(filename)
        else    
            write(fn,'(A,A,I0)') trim(filename),'.',ifile-1
        end if
        call load_single_file(trim(fn),filetype)
        Lauto = max(Lauto,real(ceiling(maxval((/maxval(x),maxval(y),maxval(z)/)))))
        npart_tot = npart_tot+npart
    end do
    
    m = nint(real(npart_tot,8)**(1.0/3.0))
    if (present(nparticles)) nparticles = m
    call toc
    write(*,'(A,I0)') '  number of files = ',nfiles
    write(string,'(I10)') npart_tot
    write(string2,'(I4)') m
    write(*,'(5A)') '  number of particles = ',trim(adjustl(string)),' (= ',trim(adjustl(string2)),'^3)'
    if (L>0) then
        write(string,'(F11.2)') L
        write(*,'(3A)') '  simulation box side length = ',trim(adjustl(string)),' (set by user)'
    else
        L = Lauto
        write(string,'(F11.2)') L
        write(*,'(3A)') '  simulation box side length = ',trim(adjustl(string)),' (determined automatically)'
    end if    
    
    ! determine grid size
    if (cellsize>0.0) then
        n = nint(L/cellsize)
    else if (cellsize<0.0) then
        n = -nint(cellsize)
    else
        n = min(512,nint(real(npart_tot)**(1.0/3.0)))
    end if
    write(string,'(I5)') n
    write(*,'(A)') '  grid size = '//trim(adjustl(string))//'^3'
    
    ! *** pre-analysis (find mean distance from cell center [-0.5..0.5])
    if (interpolation_method>0) then
       call tic('adjust particle positions')
    
       ci = 0
       cj = 0
       ck = 0
       cf = (0.0,6.283185307179586232)/L*n
       do ifile = 1,nfiles
    
           ! load particles
           if (nfiles>1) then
               write(fn,'(A,A,I0)') trim(filename),'.',ifile-1
               call load_single_file(trim(fn),filetype)
           end if
        
           ! determine mean phase
           ci = ci+sum(exp(cf*x))*npart
           cj = cj+sum(exp(cf*y))*npart
           ck = ck+sum(exp(cf*z))*npart
        
       end do
       ci = ci/npart_tot*exp((0.0,pi))
       cj = cj/npart_tot*exp((0.0,pi))
       ck = ck/npart_tot*exp((0.0,pi))
       fiavg = real(atan2(aimag(ci),real(ci))/2.0/pi,4) ! mean offset from the bin centre
       fjavg = real(atan2(aimag(cj),real(cj))/2.0/pi,4)
       fkavg = real(atan2(aimag(ck),real(ck))/2.0/pi,4)
       call toc()
       write(*,'(A,3F8.4)') '  mean offsets from cell centers',fiavg,fjavg,fkavg
       write(*,'(A,3F8.4)') '  shifted to',target_offset,target_offset,target_offset
    end if

    ! *** interpolate on a grid
    call tic('interpolate particles on grid')
    
    allocate(rho(0:n-1,0:n-1,0:n-1))
    rho = 0
    npartnew_tot = 0
    
    ! set random seed for subsampling
    if (subsample_probability<1.0) call set_seed(subsample_seed)
    
    do ifile = 1,nfiles
    
        ! load particles
        if (nfiles>1) then
            write(fn,'(A,A,I0)') trim(filename),'.',ifile-1
            call load_single_file(trim(fn),filetype)
        end if
        
        ! renormalize (x,y,z) coordinates to interval [0,n]
        x = x/L*n
        y = y/L*n
        z = z/L*n
        
        ! subsampling
        if (subsample_probability<1.0) then
            if (allocated(rnd)) deallocate(rnd)
            allocate(rnd(npart))
            call random_number(rnd)
            npartnew = 0
            do ipart = 1,npart
               if (rnd(ipart)<subsample_probability) then
                  npartnew = npartnew+1
                  x(npartnew) = x(ipart)
                  y(npartnew) = y(ipart)
                  z(npartnew) = z(ipart)
               end if
            end do
            deallocate(rnd)
            npart = npartnew
            npartnew_tot = npartnew_tot+npart
            x = x(1:npart)
            y = y(1:npart)
            z = z(1:npart)
        end if
        
        if (interpolation_method>0) then
            ! translate particles, such that, on average, they sit half-way between the centres and the edges of the cells
            ! this translation tends to keep the power spectrun and bispectrum correct. it is particularly important for
            ! the nearest neighbor algorithm, but also for high-z applications of the smoothing algorithm:
            ! e.g. line-correlation of L1000_N128_CDM1/snapshot_001
            x = modulo(x-fiavg+target_offset,real(n))
            y = modulo(y-fjavg+target_offset,real(n))
            z = modulo(z-fkavg+target_offset,real(n))
        end if
    
        ! nearest neighbor interpolation (no smoothing)
        if (abs(interpolation_method)==1) then
            do ipart = 1,npart
                i = min(n-1,floor(x(ipart))) ! index of bin
                j = min(n-1,floor(y(ipart)))
                k = min(n-1,floor(z(ipart)))
                rho(i,j,k) = rho(i,j,k)+mass(ipart)
            end do
        end if
    
        ! top-hat smoothing
        if (abs(interpolation_method)==2) then
            smoothing = 1.0/max(n,m) ! [0..1] smoothing scale in box lengths
            h = smoothing*n/2.0 ! half smoothing scale in bin lengths
            do ipart = 1,npart
                do k = floor(z(ipart)-h),floor(z(ipart)+h)
                    fk = (min(real(k+1),z(ipart)+h)-max(real(k),z(ipart)-h))/2/h
                    fk = fk*mass(ipart)
                    do j = floor(y(ipart)-h),floor(y(ipart)+h)
                        fj = (min(real(j+1),y(ipart)+h)-max(real(j),y(ipart)-h))/2/h
                        do i = floor(x(ipart)-h),floor(x(ipart)+h)
                            fi = (min(real(i+1),x(ipart)+h)-max(real(i),x(ipart)-h))/2/h
                            rho(modulo(i,n),modulo(j,n),modulo(k,n)) = rho(modulo(i,n),modulo(j,n),modulo(k,n))&
                            &+(fi*fj*fk)
                        end do
                    end do
                end do
            end do
        end if
    
    end do

    ! compute density perturbation field
    deallocate(x,y,z,mass)
    allocate(delta_r(0:n-1,0:n-1,0:n-1))
    meanrho = real(sum(real(rho,8))/real(n,8)**3,4)
    delta_r = real(rho/meanrho-1.0)
    deallocate(rho)
    
    call toc
    
    if (interpolation_method==1) then
        write(*,'(A)') '  method: nearest neighbor interpolation; conserves small-scale power,'
        write(*,'(A)') '          but can be subject to beating artefacts in all correlations'
    else if (interpolation_method==2) then
        write(*,'(A)') '  method: robust interpolation with top-hat smoothing'
        write(*,'(A,Es9.3)') '          => some small-scale power lost for k>',min(n,m)/L
    end if
    
    if (subsample_probability<1.0) then
        write(*,'(A,I0)') '  subsampled number of particles = ',npartnew_tot
    end if
    
    call toc
    
    contains
    
    subroutine load_single_file(filename,filetype)
    
        ! loads particles into
        ! x(:),y(:),z(:),mass(:)
        ! and updates the variables
        ! npart = number of particles
    
        implicit none
        character(*),intent(in) :: filename
        integer*4,intent(in)    :: filetype
        logical*4               :: loadmasses
        integer*8               :: file_size,i
        integer*4               :: bytes_per_particle
        integer*8               :: empty8
        real*4                  :: empty4
        real*8                  :: x8,y8,z8 
        logical                 :: option_8byte = .false.
        
        if (allocated(x)) deallocate(x)
        if (allocated(y)) deallocate(y)
        if (allocated(z)) deallocate(z)
        if (allocated(mass)) deallocate(mass)
        
        loadmasses = filetype<0
        
        if (abs(filetype) == 2) then ! Simple binary file
            
            ! determine number of particles
            if (loadmasses) then
                bytes_per_particle = 16
            else
                bytes_per_particle = 12
            end if
            if (option_8byte) bytes_per_particle = bytes_per_particle*2
            inquire(file=trim(filename), size=file_size)
            call check_npart(file_size/int(bytes_per_particle,8))
            npart = int(file_size/int(bytes_per_particle,8),4)
            if (int(npart,8)*int(bytes_per_particle,8).ne.file_size) then
                write(*,'(A)')
                write(*,'(A)') 'Format of input file not recognized. Consider specifying a different format using -input.'
                stop
            end if
            
            ! load particles
            allocate(x(npart),y(npart),z(npart),mass(npart))
            open(1,file=trim(filename),action='read',form='unformatted',status='old',access='stream')
            if (loadmasses) then
                read(1) (x(i),y(i),z(i),mass(i),i=1,npart)
            else
                if (option_8byte) then
                   do i=1,npart
                      read(1) x8,y8,z8 
                      x(i) = real(x8,4)
                      y(i) = real(y8,4)
                      z(i) = real(z8,4)
                   end do
                else
                   read(1) (x(i),y(i),z(i),i=1,npart)
                end if
                mass = 1.0
            end if
            close(1)
            
        else if (abs(filetype) == 3) then ! Simple ascii file
        
            ! determine number of particles
            npart = 0
            open(1,file=trim(filename),action='read',form='formatted',status='old')
            stat = 0
            do while (stat==0)
                read(1,*,IOSTAT=stat) xempty
                if (stat.ne.0) exit
                npart = npart+1
            end do
            close(1)
            call check_npart(npart*1_8)
    
            allocate(x(npart),y(npart),z(npart),mass(npart))
            open(1,file=trim(filename),action='read',form='formatted',status='old')
            if (loadmasses) then
                do i = 1,npart
                    read(1,*) x(i),y(i),z(i),mass(i)
                end do
            else
                do i = 1,npart
                    read(1,*) x(i),y(i),z(i)
                end do
                mass = 1.0
            end if
            close(1)
            
        else if (abs(filetype) == 4) then ! Gadget file binary
    
            ! determine number of particles
            open(1,file=trim(filename),action='read',form='unformatted',status='old')
            read(1) np    ! np is a 6 element array of integers*4
            npart = np(species)
            call check_npart(npart*1_8)
            
            ! load particles
            if (sum(np(1:species-1))>0) then
                read(1) (xempty,yempty,zempty,i=1,sum(np(1:species-1)))
            end if
            allocate(x(npart),y(npart),z(npart),mass(npart))
            read(1) (x(i),y(i),z(i),i=1,npart) ! x,y,z are arrays of positions, real*4
            if (loadmasses) then
                read(1) (mass(i),i=1,npart)
            else
                mass = 1.0
            end if
            close(1)
            
        else if (filetype == 5) then ! Simple binary file
    
            ! determine number of particles
            bytes_per_particle = 36
            inquire(file=trim(filename),size=file_size)
            call check_npart(file_size/int(bytes_per_particle,8))
            npart = int(file_size/int(bytes_per_particle,8),4)
            if (int(npart,8)*int(bytes_per_particle,8).ne.file_size) then
                write(*,'(A)')
                write(*,'(A)') 'Format of input file not recognized. Consider specifying a different format using -input.'
                stop
            end if
            
            ! load particles
            allocate(x(npart),y(npart),z(npart),mass(npart))
            open(1,file=trim(filename),action='read',form='unformatted',status='old',access='stream')
            read(1) (empty8,empty4,x(i),y(i),z(i),empty4,empty4,empty4,i=1,npart)
            mass = 1.0
            close(1)
            
        end if
      
    end subroutine load_single_file
    
    subroutine check_npart(guess)
    
        implicit none
        integer*8,intent(in)    :: guess
        if (guess>=int(huge(npart),8)) then
            write(*,'(A)')
            write(*,'(A)') 'No single file can contain more than 2^31 particles.'
            stop   
        end if
    
    end subroutine check_npart
    
    subroutine get_number_of_files
        implicit none
        character(len=255)  :: fn
        logical             :: single_file_exists,multiple_file_exists
        
        inquire(file=trim(filename),exist=single_file_exists)
        inquire(file=trim(filename)//'.0',exist=multiple_file_exists)
        
        if (single_file_exists.and.multiple_file_exists) then
            write(*,'(A)')
            write(*,'(A)') 'Filename ambiguous. Single and multiple files exist for'
            write(*,'(A)') trim(filename)
            stop
        end if
        
        if ((.not.single_file_exists).and.(.not.multiple_file_exists)) then
            write(*,'(A)')
            write(*,'(A)') 'File does not exist.'
            write(*,'(A)') trim(filename)
            stop
        end if
        
        if (single_file_exists) then
            nfiles = 1
        end if
        
        if (multiple_file_exists) then
            nfiles = 0
            do
                write(fn,'(A,A,I0)') trim(filename),'.',nfiles
                inquire(file=trim(fn),exist=multiple_file_exists)
                if (.not.multiple_file_exists) exit
                nfiles = nfiles+1
            end do
        end if
    end subroutine get_number_of_files
    
end subroutine make_density_perturbation_field_from_particle_file

! =================================================================================================================
! AUXILIARY SUBROUTINES
! =================================================================================================================

subroutine make_list_of_scale_lengths(n,L,rmin,scale)

    ! produces correlation lengths which are all multiples of dr=L/n

    implicit none
    
    ! variable declaration
    integer,intent(in)              :: n                ! n^dim = number of elements in delta_k
    real,intent(in)                 :: L                ! side-length of density field delta_r
    real,intent(in)                 :: rmin             ! minimal scale requested
    real,allocatable,intent(out)    :: scale(:)         ! array with different correlation functions
    real,allocatable                :: tmp(:)
    integer                         :: i,j,m
    real                            :: r,dr

    m = floor(real(n)/2)
    allocate(tmp(m))
    dr = L/real(n)
    
    j = 0
    do i = 1,m
        r = i*dr
        if (r>=rmin) then
            j = j+1
            tmp(j) = r
        end if
    end do
    
    allocate(scale(j))
    scale(:) = tmp(1:j)
    
end subroutine make_list_of_scale_lengths

subroutine make_list_of_scale_lengths_old(n,L,scale,to_cell_spacing)

    ! used in versions up to 1.15
    ! produces correlation lengths which are all multiples of dr=L/n

    implicit none
    
    ! variable declaration
    integer,intent(in)              :: n                ! n^dim = number of elements in delta_k
    real,intent(in)                 :: L                ! side-length of density field delta_r
    real,allocatable,intent(out)    :: scale(:)         ! array with different correlation functions
    logical,intent(in),optional     :: to_cell_spacing  ! if true the smallest scale is set for 2-point correlation
    real                            :: dr               ! side-length of a grid cell of delta_r
    real,allocatable                :: rlist_tmp(:)     ! temporary list of correlation lengths
    integer                         :: m                ! number of correlation lengths
    integer                         :: ikmax
    real                            :: r,ikreal
    real,parameter                  :: pi = 3.14159265

    dr = L/real(n)
    allocate(rlist_tmp(0:n))
    rlist_tmp(0) = 1e10
    m = 0
    ikmax = n/2
    if (present(to_cell_spacing)) then
        if (to_cell_spacing) then
            ikmax = n
        end if
    end if
    
    ikreal = 4
    do while (ikreal<=ikmax)
        r = nint(n/ikreal)*dr
        if ((r<rlist_tmp(m)).and.(r>=dr).and.(r>=10)) then
            m = m+1
            rlist_tmp(m) = r
        end if
        ikreal = ikreal+1
    end do
    
    allocate(scale(m))
    scale(:) = rlist_tmp(m:1:-1)
    
end subroutine make_list_of_scale_lengths_old

subroutine compute_phase_factors(delta_k,epsil_k)

    ! variable declaration
    implicit none
    complex,intent(in)              :: delta_k(:,:,:)
    complex,allocatable,intent(out) :: epsil_k(:,:,:)
    integer                         :: n,bmin,bmax
    integer                         :: i,j,k
    real                            :: limit
    
    ! computation
    n = size(delta_k,1)
    bmin = -n/2
    bmax = (n-1)/2
    
    allocate(epsil_k(bmin:bmax,bmin:bmax,bmin:bmax))
    epsil_k(bmin:bmax,bmin:bmax,bmin:bmax) = delta_k
    limit = maxval(cabs(epsil_k))*epsilon(abs(epsil_k(0,0,0)))*10
    do k = bmin,bmax
        do j = bmin,bmax
            do i = bmin,bmax
                if (cabs(epsil_k(i,j,k))>limit) then
                    epsil_k(i,j,k) = epsil_k(i,j,k)/cabs(epsil_k(i,j,k))
                else
                    epsil_k(i,j,k) = 0.0
                end if
            end do
        end do
    end do
    epsil_k(0,0,0) = 0 ! Recommended value for best results at largest correlation lengths
    !if (epsil_k(bmin,bmin,bkmin)==0.0) then
    !    write(*,'(A)') '+ WARNING: numerical accuracy below the dynamic range of the spectrum.'
    !    write(*,'(A)') '           Try to add more noise to the input density field.'
    !end if

end subroutine compute_phase_factors

subroutine set_seed(seed)
   integer,intent(in)      :: seed
   integer                 :: rnd_size
   integer*4,allocatable   :: seed_array(:)
   call random_seed(size=rnd_size)
   allocate(seed_array(rnd_size))
   seed_array = seed
   call random_seed(put=seed_array)
end subroutine set_seed

function nodot(strin) result(strout)
   ! removed dot from on the RHS of string
   implicit none
   character(*),intent(in) :: strin
   character(len=255)      :: strout
   integer*4               :: l
   l = len(trim(strin))
   if (strin(l:l)=='.') then
      strout = strin(1:l-1)
   else
      strout = strin
   end if
end function nodot

subroutine make_or_load_DFT(x,y,opt_verbose)

    implicit none
    
    character(len=*),parameter      :: filename = '.last_fourier_transform'
    real,intent(in)                 :: x(:,:,:)
    complex,allocatable,intent(out) :: y(:,:,:)
    logical,intent(in),optional     :: opt_verbose
    logical                         :: exists
    real*8                          :: checksum(3),ref(3)
    integer                         :: i,j,k
    integer                         :: bmin,bmax,bkmin,bkmax
    integer,parameter               :: fileid = 90
    logical                         :: verbose
    
    ! handle input
    if (present(opt_verbose)) then
        verbose = opt_verbose
    else
        verbose = .true.
    end if

    ! check if Fourier transform was performed before
    bmin = -size(x,1)/2
    bmax = (size(x,1)-1)/2
    bkmin = -size(x,3)/2
    bkmax = (size(x,3)-1)/2
    
    checksum(1) = sum(x**2)*log10(epsilon(abs(x(lbound(x,1),lbound(x,2),lbound(x,3)))))
    checksum(2) = sum(x(lbound(x,1):lbound(x,1)+bmax,:,:))
    checksum(3) = sum(x(:,lbound(x,2):lbound(x,2)+bmax,:))

    inquire(file=filename,exist=exists)
    if (exists) then
        open(fileid,file=filename,action='read',form='unformatted')
        read(fileid) ref
        close(fileid)
        exists = (all(checksum.eq.ref))
    end if

    if (exists) then
    
        if (verbose) call tic('load pre-computed Fourier transform') ! Loads FT from HD
        allocate(y(bmin:bmax,bmin:bmax,bkmin:bkmax))
        open(fileid,file=filename,action='read',form='unformatted')
        read(fileid) checksum
        read(fileid) (((y(i,j,k),k=bkmin,bkmax),j=bmin,bmax),i=bmin,bmax)
        close(fileid)
        
    else
    
        if (verbose) call tic('compute Fourier transform') ! Computes new FT
        call DFT(x,y)
        open(fileid,file=filename,action='write',form='unformatted')
        write(fileid) checksum
        write(fileid) (((y(i,j,k),k=bkmin,bkmax),j=bmin,bmax),i=bmin,bmax)
        close(fileid)
        
    end if

    if (verbose) call toc
    
end subroutine make_or_load_DFT

subroutine DFT(x,y)

    ! computes the centered discrete fourier transform (DFT), if the number of elements per
    ! dimension of x is odd. Otherwise, a 'pseudo-DFT' is computed where the smallest mode
    ! remains the trivial mode (q=0) rather than becoming the half mode (q=1/2).

    ! variable declaration
    implicit none
    real,intent(in)                 :: x(:,:,:)
    complex,allocatable,intent(out) :: y(:,:,:)
    real,allocatable                :: z_real(:,:,:)
    real,allocatable                :: z_imag(:,:,:)
    integer                         :: n,bmin,bmax
    integer                         :: nk,bkmin,bkmax
    
    ! initiate
    n = size(x,1)
    bmin = -n/2
    bmax = (n-1)/2
    nk = size(x,3)
    bkmin = -nk/2
    bkmax = (nk-1)/2

    ! computation
    allocate(z_real(n,n,nk)); z_real(:,:,:) = x(:,:,:)
    allocate(z_imag(n,n,nk)); z_imag(:,:,:) = 0
    call fft(z_real,z_imag,n*n*nk,n,n,-1)
    call fft(z_real,z_imag,n*n*nk,n,n*n,-1)
    call fft(z_real,z_imag,n*n*nk,nk,n*n*nk,-1)
    z_real = cshift(z_real,bmin,1)
    z_real = cshift(z_real,bmin,2)
    z_real(:,:,:) = cshift(z_real,bkmin,3)
    z_imag = cshift(z_imag,bmin,1)
    z_imag = cshift(z_imag,bmin,2)
    z_imag(:,:,:) = cshift(z_imag,bkmin,3)
    allocate(y(bmin:bmax,bmin:bmax,bkmin:bkmax))
    y(:,:,:) = cmplx(z_real,z_imag)/(real(n)**2*real(nk))
    deallocate(z_real,z_imag)
    
    ! NOTE on the use of FFTW3: FFTW3, which requries the flag -lfftw3, is not optimal here for two reasons
    ! 1) Although the execution of a Fourier transform is very fast using FFTW3, the planing
    !    itself requires a lot of time. FFTW3 is very efficient in applications where the planing
    !    is not the time critical part. In the case of the line-correlation function, where only
    !    one transform of a certain size is required, planing is becomes crucial, and FFTW3 becomes slow.
    ! 2) The planer routines in FFTW3 cannot be called from several threads at the same time.
    !    Therefore, the individual threads need to wait.
    
end subroutine DFT

subroutine IDFT(x,z_real)

    ! computes the real inverse centered discrete fourier transform (DFT) of x

    ! variable declaration
    implicit none
    complex,intent(in)              :: x(:,:,:)
    real,allocatable,intent(out)    :: z_real(:,:,:)
    real,allocatable                :: z_imag(:,:,:)
    integer                         :: n,bmin,bmax
    integer                         :: nk,bkmin,bkmax
    
    ! initiate
    n = size(x,1)
    bmin = -n/2
    bmax = (n-1)/2
    nk = size(x,3)
    bkmin = -nk/2
    bkmax = (nk-1)/2
        
    ! computation
    allocate(z_real(0:n-1,0:n-1,0:nk-1))
    allocate(z_imag(0:n-1,0:n-1,0:nk-1))
    z_real(:,:,:) = real(x(:,:,:))
    z_imag(:,:,:) = aimag(x(:,:,:))
    z_real(:,:,:) = cshift(cshift(cshift(z_real,-bmin,1),-bmin,2),-bkmin,3)
    z_imag(:,:,:) = cshift(cshift(cshift(z_imag,-bmin,1),-bmin,2),-bkmin,3)
    call fft(z_real,z_imag,n*n*nk,n,n,1)
    call fft(z_real,z_imag,n*n*nk,n,n*n,1)
    call fft(z_real,z_imag,n*n*nk,nk,n*n*nk,1)
    if (sum(abs(z_imag))/sum(abs(z_real))>0.001) print*,'WARINING in IDFT, try to use compiler option -fdefault-real-8'
    deallocate(z_imag)

end subroutine IDFT

subroutine fft(a,b,ntot,n,nspan,isn)

    ! ------------------------------------------------------------------------------------------
    ! This routine was downloaded from http://www.netlib.org/go/fft.f
    ! and dimensions have been changed to compute FFTs of 3D arrays up to a size of 4096^3 cells
    ! It should be referenced as: "An N-dimensional Fortran mixed-radix FFT from the 'Golden Oldies' library,
    ! by R. C. Singleton, Stanford Research Institute, Sept. 1968, http://www.netlib.org/go/index.html'
    ! ------------------------------------------------------------------------------------------
    ! We use this FFT for the following reasons
    ! + no initialization required, and therefore very efficient in the case where many arrays of different
    !   size need to be transformed. (on a single call faster than FFTW3 for n^3 < 300^3 ~ 10^7)
    ! + fully compatible with openmp in the sense that it can be called from multiple threads
    !   simultaneously (not the case for FFTW3)
    ! + no additional library required (unlike FFTW3, FFTPACK5, ...)
    ! + compatible with any real kind (namely single and double) without changing the code
    ! + free license
    ! ------------------------------------------------------------------------------------------
    ! Date: Wed, 9 Aug 1995 09:38:49 -0400
    ! From: ldm@apollo.numis.nwu.edu
    !  multivariate complex fourier transform, computed in place
    !    using mixed-radix fast fourier transform algorithm.
    !  by r. c. singleton, stanford research institute, sept. 1968
    !  arrays a and b originally hold the real and imaginary
    !    components of the data, and return the real and
    !    imaginary components of the resulting fourier coefficients.
    !  multivariate data is indexed according to the fortran
    !    array element successor function, without limit
    !    on the number of implied multiple subscripts.
    !    the subroutine is called once for each variate.
    !    the calls for a multivariate transform may be in any order.
    !  ntot is the total number of complex data values.
    !  n is the dimension of the current variable.
    !  nspan/n is the spacing of consecutive data values
    !    while indexing the current variable.
    !  the sign of isn determines the sign of the complex
    !    exponential, and the magnitude of isn is normally one.
    !  a tri-variate transform with a(n1,n2,n3), b(n1,n2,n3)
    !    is computed by
    !      call fft(a,b,n1*n2*n3,n1,n1,1)
    !      call fft(a,b,n1*n2*n3,n2,n1*n2,1)
    !      call fft(a,b,n1*n2*n3,n3,n1*n2*n3,1)
    !  for a single-variate transform,
    !    ntot = n = nspan = (number of complex data values), e.g.
    !      call fft(a,b,n,n,n,1)
    !  the data can alternatively be stored in a single complex array c
    !    in standard fortran fashion, i.e. alternating real and imaginary
    !    parts. then with most fortran compilers, the complex array c can
    !    be equivalenced to a real array a, the magnitude of isn changed
    !    to two to give correct indexing increment, and a(1) and a(2) used
    !    to pass the initial addresses for the sequences of real and
    !    imaginary values, e.g.
    !       complex c(ntot)
    !       real    a(2*ntot)
    !       equivalence (c(1),a(1))
    !       call fft(a(1),a(2),ntot,n,nspan,2)
    !  arrays at(maxf), ck(maxf), bt(maxf), sk(maxf), and np(maxp)
    !    are used for temporary storage.  if the available storage
    !    is insufficient, the program is terminated by a stop.
    ! ------------------------------------------------------------------------------------------
    dimension a(*),b(*) ! pointers to input/output arrays
    !    maxf must be .ge. the maximum prime factor of n.
    !    maxp must be .gt. the number of prime factors of n.
    !    in addition, if the square-free portion k of n has two or
    !    more prime factors, then maxp must be .ge. k-1.
    !  array storage in nfac for a maximum of 15 prime factors of n.
    !  if n has more than one square-free factor, the product of the
    !    square-free factors must be .le. 4096
    dimension nfac(30),np(4095)
    !  array storage for maximum prime factor of 4099
    dimension at(4099),ck(4099),bt(4099),sk(4099)
    equivalence (i,ii)
    !  the following two constants should agree with the array dimensions.
    maxp=4095
    maxf=4093
    ! ------------------------------------------------------------------------------------------
    
    if(n .lt. 2) return
    inc=isn
    c72=0.30901699437494742
    s72=0.95105651629515357
    s120=0.86602540378443865
    rad=6.2831853071796
    if(isn .ge. 0) go to 10
    s72=-s72
    s120=-s120
    rad=-rad
    inc=-inc
    10 nt=inc*ntot
    ks=inc*nspan
    kspan=ks
    nn=nt-inc
    jc=ks/n
    radf=rad*float(jc)*0.5
    i=0
    jf=0
    !  determine the factors of n
    m=0
    k=n
    go to 20
    15 m=m+1
    nfac(m)=4
    k=k/16
    20 if(k-(k/16)*16 .eq. 0) go to 15
    j=3
    jj=9
    go to 30
    25 m=m+1
    nfac(m)=j
    k=k/jj
    30 if(mod(k,jj) .eq. 0) go to 25
    j=j+2
    jj=j**2
    if(jj .le. k) go to 30
    if(k .gt. 4) go to 40
    kt=m
    nfac(m+1)=k
    if(k .ne. 1) m=m+1
    go to 80
    40 if(k-(k/4)*4 .ne. 0) go to 50
    m=m+1
    nfac(m)=2
    k=k/4
    50 kt=m
    j=2
    60 if(mod(k,j) .ne. 0) go to 70
    m=m+1
    nfac(m)=j
    k=k/j
    70 j=((j+1)/2)*2+1
    if(j .le. k) go to 60
    80 if(kt .eq. 0) go to 100
    j=kt
    90 m=m+1
    nfac(m)=nfac(j)
    j=j-1
    if(j .ne. 0) go to 90
    !  compute fourier transform
    100 sd=radf/float(kspan)
    cd=2.0*sin(sd)**2
    sd=sin(sd+sd)
    kk=1
    i=i+1
    if(nfac(i) .ne. 2) go to 400
    !  transform for factor of 2 (including rotation factor)
    kspan=kspan/2
    k1=kspan+2
    210 k2=kk+kspan
    ak=a(k2)
    bk=b(k2)
    a(k2)=a(kk)-ak
    b(k2)=b(kk)-bk
    a(kk)=a(kk)+ak
    b(kk)=b(kk)+bk
    kk=k2+kspan
    if(kk .le. nn) go to 210
    kk=kk-nn
    if(kk .le. jc) go to 210
    if(kk .gt. kspan) go to 800
    220 c1=1.0-cd
    s1=sd
    230 k2=kk+kspan
    ak=a(kk)-a(k2)
    bk=b(kk)-b(k2)
    a(kk)=a(kk)+a(k2)
    b(kk)=b(kk)+b(k2)
    a(k2)=c1*ak-s1*bk
    b(k2)=s1*ak+c1*bk
    kk=k2+kspan
    if(kk .lt. nt) go to 230
    k2=kk-nt
    c1=-c1
    kk=k1-k2
    if(kk .gt. k2) go to 230
    ak=c1-(cd*c1+sd*s1)
    s1=(sd*c1-cd*s1)+s1
    c1=2.0-(ak**2+s1**2)
    s1=c1*s1
    c1=c1*ak
    kk=kk+jc
    if(kk .lt. k2) go to 230
    k1=k1+inc+inc
    kk=(k1-kspan)/2+jc
    if(kk .le. jc+jc) go to 220
    go to 100
    !  transform for factor of 3 (optional code)
    320 k1=kk+kspan
    k2=k1+kspan
    ak=a(kk)
    bk=b(kk)
    aj=a(k1)+a(k2)
    bj=b(k1)+b(k2)
    a(kk)=ak+aj
    b(kk)=bk+bj
    ak=-0.5*aj+ak
    bk=-0.5*bj+bk
    aj=(a(k1)-a(k2))*s120
    bj=(b(k1)-b(k2))*s120
    a(k1)=ak-bj
    b(k1)=bk+aj
    a(k2)=ak+bj
    b(k2)=bk-aj
    kk=k2+kspan
    if(kk .lt. nn) go to 320
    kk=kk-nn
    if(kk .le. kspan) go to 320
    go to 700
    !  transform for factor of 4
    400 if(nfac(i) .ne. 4) go to 600
    kspnn=kspan
    kspan=kspan/4
    410 c1=1.0
    s1=0
    420 k1=kk+kspan
    k2=k1+kspan
    k3=k2+kspan
    akp=a(kk)+a(k2)
    akm=a(kk)-a(k2)
    ajp=a(k1)+a(k3)
    ajm=a(k1)-a(k3)
    a(kk)=akp+ajp
    ajp=akp-ajp
    bkp=b(kk)+b(k2)
    bkm=b(kk)-b(k2)
    bjp=b(k1)+b(k3)
    bjm=b(k1)-b(k3)
    b(kk)=bkp+bjp
    bjp=bkp-bjp
    if(isn .lt. 0) go to 450
    akp=akm-bjm
    akm=akm+bjm
    bkp=bkm+ajm
    bkm=bkm-ajm
    if(s1 .eq. 0) go to 460
    430 a(k1)=akp*c1-bkp*s1
    b(k1)=akp*s1+bkp*c1
    a(k2)=ajp*c2-bjp*s2
    b(k2)=ajp*s2+bjp*c2
    a(k3)=akm*c3-bkm*s3
    b(k3)=akm*s3+bkm*c3
    kk=k3+kspan
    if(kk .le. nt) go to 420
    440 c2=c1-(cd*c1+sd*s1)
    s1=(sd*c1-cd*s1)+s1
    c1=2.0-(c2**2+s1**2)
    s1=c1*s1
    c1=c1*c2
    c2=c1**2-s1**2
    s2=2.0*c1*s1
    c3=c2*c1-s2*s1
    s3=c2*s1+s2*c1
    kk=kk-nt+jc
    if(kk .le. kspan) go to 420
    kk=kk-kspan+inc
    if(kk .le. jc) go to 410
    if(kspan .eq. jc) go to 800
    go to 100
    450 akp=akm+bjm
    akm=akm-bjm
    bkp=bkm-ajm
    bkm=bkm+ajm
    if(s1 .ne. 0) go to 430
    460 a(k1)=akp
    b(k1)=bkp
    a(k2)=ajp
    b(k2)=bjp
    a(k3)=akm
    b(k3)=bkm
    kk=k3+kspan
    if(kk .le. nt) go to 420
    go to 440
    !  transform for factor of 5 (optional code)
    510 c2=c72**2-s72**2
    s2=2.0*c72*s72
    520 k1=kk+kspan
    k2=k1+kspan
    k3=k2+kspan
    k4=k3+kspan
    akp=a(k1)+a(k4)
    akm=a(k1)-a(k4)
    bkp=b(k1)+b(k4)
    bkm=b(k1)-b(k4)
    ajp=a(k2)+a(k3)
    ajm=a(k2)-a(k3)
    bjp=b(k2)+b(k3)
    bjm=b(k2)-b(k3)
    aa=a(kk)
    bb=b(kk)
    a(kk)=aa+akp+ajp
    b(kk)=bb+bkp+bjp
    ak=akp*c72+ajp*c2+aa
    bk=bkp*c72+bjp*c2+bb
    aj=akm*s72+ajm*s2
    bj=bkm*s72+bjm*s2
    a(k1)=ak-bj
    a(k4)=ak+bj
    b(k1)=bk+aj
    b(k4)=bk-aj
    ak=akp*c2+ajp*c72+aa
    bk=bkp*c2+bjp*c72+bb
    aj=akm*s2-ajm*s72
    bj=bkm*s2-bjm*s72
    a(k2)=ak-bj
    a(k3)=ak+bj
    b(k2)=bk+aj
    b(k3)=bk-aj
    kk=k4+kspan
    if(kk .lt. nn) go to 520
    kk=kk-nn
    if(kk .le. kspan) go to 520
    go to 700
    !  transform for odd factors
    600 k=nfac(i)
    kspnn=kspan
    kspan=kspan/k
    if(k .eq. 3) go to 320
    if(k .eq. 5) go to 510
    if(k .eq. jf) go to 640
    jf=k
    s1=rad/float(k)
    c1=cos(s1)
    s1=sin(s1)
    if(jf .gt. maxf) go to 998
    ck(jf)=1.0
    sk(jf)=0.0
    j=1
    630 ck(j)=ck(k)*c1+sk(k)*s1
    sk(j)=ck(k)*s1-sk(k)*c1
    k=k-1
    ck(k)=ck(j)
    sk(k)=-sk(j)
    j=j+1
    if(j .lt. k) go to 630
    640 k1=kk
    k2=kk+kspnn
    aa=a(kk)
    bb=b(kk)
    ak=aa
    bk=bb
    j=1
    k1=k1+kspan
    650 k2=k2-kspan
    j=j+1
    at(j)=a(k1)+a(k2)
    ak=at(j)+ak
    bt(j)=b(k1)+b(k2)
    bk=bt(j)+bk
    j=j+1
    at(j)=a(k1)-a(k2)
    bt(j)=b(k1)-b(k2)
    k1=k1+kspan
    if(k1 .lt. k2) go to 650
    a(kk)=ak
    b(kk)=bk
    k1=kk
    k2=kk+kspnn
    j=1
    660 k1=k1+kspan
    k2=k2-kspan
    jj=j
    ak=aa
    bk=bb
    aj=0.0
    bj=0.0
    k=1
    670 k=k+1
    ak=at(k)*ck(jj)+ak
    bk=bt(k)*ck(jj)+bk
    k=k+1
    aj=at(k)*sk(jj)+aj
    bj=bt(k)*sk(jj)+bj
    jj=jj+j
    if(jj .gt. jf) jj=jj-jf
    if(k .lt. jf) go to 670
    k=jf-j
    a(k1)=ak-bj
    b(k1)=bk+aj
    a(k2)=ak+bj
    b(k2)=bk-aj
    j=j+1
    if(j .lt. k) go to 660
    kk=kk+kspnn
    if(kk .le. nn) go to 640
    kk=kk-nn
    if(kk .le. kspan) go to 640
    !  multiply by rotation factor (except for factors of 2 and 4)
    700 if(i .eq. m) go to 800
    kk=jc+1
    710 c2=1.0-cd
    s1=sd
    720 c1=c2
    s2=s1
    kk=kk+kspan
    730 ak=a(kk)
    a(kk)=c2*ak-s2*b(kk)
    b(kk)=s2*ak+c2*b(kk)
    kk=kk+kspnn
    if(kk .le. nt) go to 730
    ak=s1*s2
    s2=s1*c2+c1*s2
    c2=c1*c2-ak
    kk=kk-nt+kspan
    if(kk .le. kspnn) go to 730
    c2=c1-(cd*c1+sd*s1)
    s1=s1+(sd*c1-cd*s1)
    c1=2.0-(c2**2+s1**2)
    s1=c1*s1
    c2=c1*c2
    kk=kk-kspnn+jc
    if(kk .le. kspan) go to 720
    kk=kk-kspan+jc+inc
    if(kk .le. jc+jc) go to 710
    go to 100
    !  permute the results to normal order---done in two stages
    !  permutation for square factors of n
    800 np(1)=ks
    if(kt .eq. 0) go to 890
    k=kt+kt+1
    if(m .lt. k) k=k-1
    j=1
    np(k+1)=jc
    810 np(j+1)=np(j)/nfac(j)
    np(k)=np(k+1)*nfac(j)
    j=j+1
    k=k-1
    if(j .lt. k) go to 810
    k3=np(k+1)
    kspan=np(2)
    kk=jc+1
    k2=kspan+1
    j=1
    if(n .ne. ntot) go to 850
    !  permutation for single-variate transform (optional code)
    820 ak=a(kk)
    a(kk)=a(k2)
    a(k2)=ak
    bk=b(kk)
    b(kk)=b(k2)
    b(k2)=bk
    kk=kk+inc
    k2=kspan+k2
    if(k2 .lt. ks) go to 820
    830 k2=k2-np(j)
    j=j+1
    k2=np(j+1)+k2
    if(k2 .gt. np(j)) go to 830
    j=1
    840 if(kk .lt. k2) go to 820
    kk=kk+inc
    k2=kspan+k2
    if(k2 .lt. ks) go to 840
    if(kk .lt. ks) go to 830
    jc=k3
    go to 890
    !  permutation for multivariate transform
    850 k=kk+jc
    860 ak=a(kk)
    a(kk)=a(k2)
    a(k2)=ak
    bk=b(kk)
    b(kk)=b(k2)
    b(k2)=bk
    kk=kk+inc
    k2=k2+inc
    if(kk .lt. k) go to 860
    kk=kk+ks-jc
    k2=k2+ks-jc
    if(kk .lt. nt) go to 850
    k2=k2-nt+kspan
    kk=kk-nt+jc
    if(k2 .lt. ks) go to 850
    870 k2=k2-np(j)
    j=j+1
    k2=np(j+1)+k2
    if(k2 .gt. np(j)) go to 870
    j=1
    880 if(kk .lt. k2) go to 850
    kk=kk+jc
    k2=kspan+k2
    if(k2 .lt. ks) go to 880
    if(kk .lt. ks) go to 870
    jc=k3
    890 if(2*kt+1 .ge. m) return
    kspnn=np(kt+1)
    !  permutation for square-free factors of n
    j=m-kt
    nfac(j+1)=1
    900 nfac(j)=nfac(j)*nfac(j+1)
    j=j-1
    if(j .ne. kt) go to 900
    kt=kt+1
    nn=nfac(kt)-1
    if(nn .gt. maxp) go to 998
    jj=0
    j=0
    go to 906
    902 jj=jj-k2
    k2=kk
    k=k+1
    kk=nfac(k)
    904 jj=kk+jj
    if(jj .ge. k2) go to 902
    np(j)=jj
    906 k2=nfac(kt)
    k=kt+1
    kk=nfac(k)
    j=j+1
    if(j .le. nn) go to 904
    !  determine the permutation cycles of length greater than 1
    j=0
    go to 914
    910 k=kk
    kk=np(k)
    np(k)=-kk
    if(kk .ne. j) go to 910
    k3=kk
    914 j=j+1
    kk=np(j)
    if(kk .lt. 0) go to 914
    if(kk .ne. j) go to 910
    np(j)=-j
    if(j .ne. nn) go to 914
    maxf=inc*maxf
    !  reorder a and b, following the permutation cycles
    go to 950
    924 j=j-1
    if(np(j) .lt. 0) go to 924
    jj=jc
    926 kspan=jj
    if(jj .gt. maxf) kspan=maxf
    jj=jj-kspan
    k=np(j)
    kk=jc*k+ii+jj
    k1=kk+kspan
    k2=0
    928 k2=k2+1
    at(k2)=a(k1)
    bt(k2)=b(k1)
    k1=k1-inc
    if(k1 .ne. kk) go to 928
    932 k1=kk+kspan
    k2=k1-jc*(k+np(k))
    k=-np(k)
    936 a(k1)=a(k2)
    b(k1)=b(k2)
    k1=k1-inc
    k2=k2-inc
    if(k1 .ne. kk) go to 936
    kk=k2
    if(k .ne. j) go to 932
    k1=kk+kspan
    k2=0
    940 k2=k2+1
    a(k1)=at(k2)
    b(k1)=bt(k2)
    k1=k1-inc
    if(k1 .ne. kk) go to 940
    if(jj .ne. 0) go to 926
    if(j .ne. 1) go to 924
    950 j=k3+1
    nt=nt-kspnn
    ii=nt-inc+1
    if(nt .ge. 0) go to 924
    return
    !  error finish, insufficient array storage
    998 isn=0
    print 999
    stop
    999 format(44h0array bounds exceeded within subroutine fft)

end subroutine fft

subroutine tic(strin)

    implicit none
    character(len=*),intent(in) :: strin
    character(len=70)           :: printstr
    character(len=64)           :: str
    
    if (.not.allocated(tstart)) then
        write(*,'(A)') '----------------------------------------------------------------------------------'
        allocate(tstart(2))
        call system_clock(tstart(1),trate)
        tindex = 0
    end if
    
    if (tindex==2) then
        call toc
    end if
    
    tindex = tindex+1
    
    call system_clock(tstart(tindex))
    
    if (tindex==2) then
        if (len(strin)>64) then
            str = 'load file ...'//strin(len(strin)-50:len(strin))
        else
            str = trim(adjustl(strin))
        end if
        if (index(str,'=')==0) then
            write(printstr,'(3A)') '+ ',trim(str),' ...'
            write(*,'(A)',advance='no') adjustl(printstr)
            time_measurement = .true.
        else
            write(printstr,'(2A)') '+ ',trim(str)
            write(*,'(A)') adjustl(printstr)
            time_measurement = .false.
        end if
    end if
    
    if (tindex==1) then
        write(*,'(A)') trim(adjustl(strin))
    end if

end subroutine tic

subroutine toc

    implicit none
    integer :: tstop
    
    call system_clock(tstop)
    if (tindex==1) then
        write(*,'(A,F11.2,A)') '= total time taken                                                   ', &
        & real(tstop-tstart(tindex))/trate,'s'
        write(*,'(A)') '----------------------------------------------------------------------------------'
    else
        if (time_measurement) then
            write(*,'(F10.2,A)') real(tstop-tstart(tindex))/trate,'s'
        end if
    end if
    tindex = tindex-1

end subroutine toc

end module module_correlation_functions