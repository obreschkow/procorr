program procorr

    use module_correlation_functions
    
    implicit none
    
    character(*),parameter  :: version = '1.15'
    character(len=255)      :: filename,opt_type,opt_value
    logical                 :: opt_compute_p    ! isotropic power spectrum in CAMB normalization
    logical                 :: opt_compute_x    ! isotropic 2-point auto correlation
    logical                 :: opt_compute_y    ! isotropic 3-point auto correlation on a line
    logical                 :: opt_compute_b    ! isotropic bispectrum
    logical                 :: opt_compute_l    ! line correlation (with new definition)
    logical                 :: opt_compute_a    ! line correlation (with new definition), plus Edgeworth approximation
    logical                 :: opt_compute_d    ! density perturbation field delta(r) on regular grid
    logical                 :: opt_compute_e    ! phase field epsilon(r) on regular grid
    integer                 :: opt_input        ! input type
    real                    :: opt_sidelength   ! side length of periodic box, only needed if opt_input = 2
    integer                 :: opt_species      ! gadget particle species
    integer                 :: opt_ncells       ! number of cells aside used to discretize the cubic/square simulation box
    integer                 :: opt_method       ! [1,2] algorithm used for line correlation
    integer                 :: opt_interpolation! [1,2] interpolation method (1 = nearest neighbor, 2 = top-hat)
    real                    :: opt_accuracy     ! accuracy/completeness for bispectrum, line correlation and 3-point correlation
    logical                 :: opt_overwrite    ! whether to overwrite existing files
    logical                 :: opt_progress     ! whether to display counting progress on screen
    real                    :: opt_probability  ! subsampling probability
    integer*4               :: opt_seed         ! seed for random number generator used for subsampling if opt_probability<1
    logical                 :: ok
    integer                 :: i,j,n
    real                    :: L
    real,allocatable        :: scale(:),xi(:)
    real,allocatable        :: scalek(:),p(:)
    real,allocatable        :: value(:),error(:),value_edgeworth(:),error_edgeworth(:)
    real,allocatable        :: bispectrum(:,:),angle(:),kleg(:)
    integer                 :: iangle,ikleg
    real,allocatable        :: delta_r(:,:,:)
    real,allocatable        :: epsil_r(:,:,:)
    complex,allocatable     :: delta_k(:,:,:)
    complex,allocatable     :: epsil_k(:,:,:)
    real,parameter          :: value_ref(13) = (/1.46115053,1.36806214,1.01544416,0.725860417,0.563088596,0.462437660,0.365323305, &
                                                & 0.275995672,0.257817090,0.186858490,0.142859638,0.128057718,8.73243734E-02/)
    integer*8               :: performance(2)
    character(len=8)        :: str_acc
    character(len=100)      :: rwarning,kwarning
    real,parameter          :: pi = 3.141592653589793
    integer*4               :: ncells
    integer*8               :: nparticles
    
    character(len=255)      :: cwd, me, mydir, path
    integer                 :: scanVal 
    logical                 :: back=.true.
    logical                 :: file_exists
    
    ! Version verification
    if (version.ne.module_version) then
        write(*,'(A)') 'ERROR: procorr and module_correlation_functions do not have the same version.'
        stop
    end if
    
    ! find path of procorr
    call get_command_argument(0,me) 
    scanVal = scan (me, '/', back) 
    if (scanVal .gt. 0) then 
       mydir=me(2:scanVal) 
    else 
       call getcwd(mydir)    
    endif
    call getcwd(cwd)
    path = trim(cwd)//trim(mydir)
    
    ! default options
    opt_compute_p = .true.
    opt_compute_x = .false.
    opt_compute_y = .false.
    opt_compute_b = .false.
    opt_compute_l = .false.
    opt_compute_a = .false.
    opt_compute_d = .false.
    opt_compute_e = .false.
    opt_input = 2
    opt_sidelength = 0
    opt_species = 2
    opt_ncells = 0 ! means automatic cell size
    opt_method = 0
    opt_interpolation = 1
    opt_accuracy = 0
    opt_overwrite = .true.
    opt_progress = .false.
    opt_probability = 1.0
    opt_seed = 1
    
    ! handel input
    n = iargc() ! number of arguments
    if (n<1) then
        write(*,'(A)') 'ERROR: Filename missing.'
        write(*,'(A)') 'General use is ./procorr inputfilename [-option argument], for example:'
        write(*,'(A)') './procorr inputfilename -output pxl'
        write(*,'(A)') 'The avilable outputs are'
        write(*,'(A)') 'p: powerspectrum (in CAMB normalization)'
        write(*,'(A)') 'x: 2-point auto correlation'
        write(*,'(A)') 'y: 3-point auto correlation (still developer version)'
        write(*,'(A)') 'b: bispectrum'
        write(*,'(A)') 'l: line-correlation l(r) as in Obreschkow+2013, but with |k|,|q|,|k-q|<2*pi/r'
        write(*,'(A)') 'a: like l, but additionally calculate l(r) in the Edgeworth approximation'
        write(*,'(A)') 'd: density perturbation field delta(r) on regular grid'
        write(*,'(A)') 'e: phase field epsilon(r) on regular grid'
        write(*,'(A)') 'For all additional options, please see the README file.'
        stop
    end if
    
    call getarg(1,filename)
    filename = trim(nodot(filename))
    
    if (trim(filename)=='-test') then
    
        if ((n.ne.1).and.(n.ne.3)) then
            write(*,'(A)') 'ERROR: Wrong arguments in test run. Use either'
            write(*,'(A)') './procorr -test'
            write(*,'(A)') 'or'
            write(*,'(A)') './procorr -test -progress y'
            stop
        else
            if (n==3) then
                call getarg(2,opt_type)
                call getarg(3,opt_value)
                if (trim(opt_type)=='-progress') then
                    if ((opt_value=='n').or.(opt_value=='N')) then
                        opt_progress = .false.
                    else
                        if ((opt_value=='y').or.(opt_value=='Y')) then
                            opt_progress = .true.
                        else
                            write(*,'(A)') 'ERROR: the option progress only takes the arguments Y or N.'
                            stop
                        end if
                    end if
                else
                    write(*,'(A)') 'ERROR: Wrong arguments in test run. Use either'
                    write(*,'(A)') './procorr -test'
                    write(*,'(A)') 'or'
                    write(*,'(A)') './procorr -test -progress y'
                    stop
                end if
            end if
        end if
    
        call set_opt_progress(opt_progress)
        
        file_exists = .false.
        
        if (.not.file_exists) then
            filename = 'cdm_redshift0'
            inquire(file=trim(filename), exist=file_exists)
            write(*,*) trim(filename)
        end if
        if (.not.file_exists) then
            filename = trim(mydir)//'cdm_redshift0'
            inquire(file=trim(filename), exist=file_exists)
            write(*,*) trim(filename)
        end if
        if (.not.file_exists) then
            filename = '/'//trim(mydir)//'cdm_redshift0'
            inquire(file=trim(filename), exist=file_exists)
            write(*,*) trim(filename)
        end if
        if (.not.file_exists) then
            filename = trim(path)//'cdm_redshift0'
            inquire(file=trim(filename), exist=file_exists)
            write(*,*) trim(filename)
        end if
        if (.not.file_exists) then
            write(*,'(A)') 'ERROR: Could not find file cdm_redshift0.'
            write(*,'(A)') 'Try to run procorr using the full path:'
            write(*,'(A)') '> /.../.../procorr  -test'
            stop
        end if
        
        L = 0
        call make_density_perturbation_field_from_particle_file(trim(filename),4,2,-83.0,1,delta_r,L,1.0,1)
        call compute_lic_deterministic(delta_r,L,scale,value,opt_performance=performance)
        ! compare value to reference value ...
        ! write(*,*) value
        if (any(abs(value-value_ref)/value_ref>1e-3)) then
            write(*,'(A)') 'ERROR: TEST FAILED. PLEASE CONTACT danail.obreschkow@icrar.org FOR SUPPORT.'
        else
            write(*,'(A)') 'TEST SUCCESSFUL.'
        end if
        write(*,'(A,A)') '+ Total number of terms in sum_(k,q):   ',longint2str(performance(1))
        write(*,'(A,A)') '+ Number of terms evaluated per second: ',longint2str(performance(2))
        write(*,'(A)') '----------------------------------------------------------------------------------'
        stop
    end if
    
    if ((trim(filename)=='-version').or.(trim(filename)=='-v')) then
        write(*,'(A,A,A)') 'This is ProCorr Version ',version
        stop
    end if
    
    if (n>1) then
        if (n/2.eq.n/2.0) then
            write(*,'(A)') 'ERROR: Every option "-option" must have exactly one value.'
            stop
        end if
        do i = 2,n,2
            call getarg(i,opt_type)
            call getarg(i+1,opt_value)
            select case (trim(opt_type))
                case ('-input')
                    if (opt_value(1:1)=='m') then
                        read(opt_value(2:len(trim(opt_value))),*) opt_input
                        opt_input = -opt_input
                    else 
                        read(opt_value,*) opt_input
                    end if
                    if (.not.((opt_input==1).or.(abs(opt_input)==2).or.&
                    &(abs(opt_input)==3).or.(abs(opt_input)==4).or.(opt_input==5))) then
                        write(*,'(A)') 'ERROR: input must be 1, 2, 3, 4, 5, m2, m3, m4.'
                        stop
                    end if
                case ('-sidelength')
                    read(opt_value,*) opt_sidelength
                    if (opt_sidelength<=0) then
                        write(*,'(A)') 'ERROR: sidelength must be larger than 0.'
                        stop
                    end if
                case ('-species')
                    read(opt_value,*) opt_species
                case ('-ncells')
                    if (opt_input==1) then
                        write(*,'(A)') 'ERROR: ncells cannot be specified for input type 1.'
                        stop                    
                    end if
                    read(opt_value,*) opt_ncells
                    if ((opt_ncells<5).or.(opt_ncells>4096)) then
                        write(*,'(A)') 'ERROR: ncells must be between 5 and 4096 (on a normal desktop <= 512).'
                        stop
                    end if
                case ('-method')
                    read(opt_value,*) opt_method
                    if ((opt_method.ne.1).and.(opt_method.ne.2)) then
                        write(*,'(A)') 'ERROR: method must be 1 or 2, for stochastic or deterministic algorithm, respectively.'
                        stop
                    end if
                case ('-interpolation')
                    if (opt_value(1:1)=='n') then
                        read(opt_value(2:len(trim(opt_value))),*) opt_interpolation
                        opt_interpolation = -opt_interpolation
                    else 
                        read(opt_value,*) opt_interpolation
                    end if
                    if ((abs(opt_interpolation).ne.1).and.(abs(opt_interpolation).ne.2)) then
                        write(*,'(A)') 'ERROR: interpolation must be 1 (n1) or 2 (n2), for'
                        write(*,'(A)') 'nearest neighbor or top-hat interpolation, respectively.'
                        stop
                    end if
                case ('-accuracy')
                    read(opt_value,*) opt_accuracy
                    if (opt_accuracy<=0) then
                        write(*,'(A)') 'ERROR: accuracy must be greater than 0.'
                        stop
                    end if
                case ('-probability')
                    if (opt_input==1) then
                        write(*,'(A)') 'ERROR: probability cannot be specified for input type 1.'
                        stop                    
                    end if
                    read(opt_value,*) opt_probability
                    if (opt_probability<=0) then
                        write(*,'(A)') 'ERROR: probability must be greater than 0.'
                        stop
                    end if
                    if (opt_probability>10) then
                        write(*,'(A)') 'ERROR: probability cannot be larger than 1.'
                        stop
                    end if
                case ('-seed')
                    if (opt_input==1) then
                        write(*,'(A)') 'ERROR: seed cannot be specified for input type 1.'
                        stop                    
                    end if
                    read(opt_value,*) opt_seed
                    if (opt_seed<0) then
                        write(*,'(A)') 'ERROR: seed must be a natural integer.'
                        stop
                    end if
                case ('-output')
                    opt_compute_p = .false.
                    opt_value = adjustl(opt_value)
                    do j = 1,len(trim(opt_value))
                        select case(opt_value(j:j))
                            case('p')
                                opt_compute_p = .true.
                            case('x')
                                opt_compute_x = .true.
                            case('y')
                                opt_compute_y = .true.
                            case('b')
                                opt_compute_b = .true.
                            case('l')
                                opt_compute_l = .true.
                            case('a')
                                opt_compute_a = .true.
                            case('d')
                                opt_compute_d = .true.
                            case('e')
                                opt_compute_e = .true.
                            case default
                                write(*,'(A)') 'ERROR: '//opt_value(j:j)//' is an unknown output specifier. Use a subset of pxblde.'
                                stop
                            end select
                    end do
                case ('-overwrite')
                    if ((opt_value=='n').or.(opt_value=='N')) then
                        opt_overwrite = .false.
                    else
                        if ((opt_value=='y').or.(opt_value=='Y')) then
                            opt_overwrite = .true.
                        else
                            write(*,'(A)') 'ERROR: the option overwrite only takes the arguments Y or N.'
                            stop
                        end if
                    end if
                case ('-progress')
                    if ((opt_value=='n').or.(opt_value=='N')) then
                        opt_progress = .false.
                    else
                        if ((opt_value=='y').or.(opt_value=='Y')) then
                            opt_progress = .true.
                        else
                            write(*,'(A)') 'ERROR: the option progress only takes the arguments Y or N.'
                            stop
                        end if
                    end if
                case default
                    write(*,'(A)') 'ERROR: '//trim(opt_type)//' is an unknown option.'
                    stop
            end select
        end do
    end if
    
    call set_opt_progress(opt_progress)
    
    if ((.not.opt_compute_l).and.(opt_method.ne.0)) then
        write(*,'(A)') 'ERROR: option -method can only be used in conjunction with -output l.'
        stop
    end if
    
    if ((.not.opt_compute_l).and.(.not.opt_compute_b).and.(.not.opt_compute_y).and.&
        &(.not.opt_compute_a).and.(opt_accuracy.ne.0)) then
        write(*,'(A)') 'ERROR: option -accuracy can only be used in conjunction with -output b, y, l, or a.'
        stop
    end if
    
    ! check if files already exist
    if (.not.opt_overwrite) then
        if (opt_compute_p) then
            inquire(file=trim(filename)//'_p.txt',exist=ok)
            opt_compute_p = .not.ok
        end if
        if (opt_compute_x) then
            inquire(file=trim(filename)//'_x.txt',exist=ok)
            opt_compute_x = .not.ok
        end if
        if (opt_compute_y) then
            inquire(file=trim(filename)//'_y.txt',exist=ok)
            opt_compute_y = .not.ok
        end if
        if (opt_compute_b) then
            inquire(file=trim(filename)//'_b.txt',exist=ok)
            opt_compute_b = .not.ok
        end if
        if (opt_compute_l) then
            inquire(file=trim(filename)//'_l.txt',exist=ok)
            opt_compute_l = .not.ok
        end if
        if (opt_compute_a) then
            inquire(file=trim(filename)//'_a.txt',exist=ok)
            opt_compute_a = .not.ok
        end if
        if (opt_compute_d) then
            inquire(file=trim(filename)//'_d.txt',exist=ok)
            opt_compute_d = .not.ok
        end if
        if (opt_compute_e) then
            inquire(file=trim(filename)//'_e.txt',exist=ok)
            opt_compute_e = .not.ok
        end if
    end if

    ! convert Gadget output into a density-field on a regular grid
    if ((opt_compute_p).or.(opt_compute_x).or.(opt_compute_y).or.(opt_compute_b).or. &
        (opt_compute_l).or.(opt_compute_a).or.(opt_compute_d).or.(opt_compute_e)) then
        kwarning = 'none'
        rwarning = 'none'
        if (opt_input==1) then
            if (opt_sidelength<=0) then
                write(*,'(A)') 'ERROR: for this input type, a box side length must be specified using -sidelength ###.'
                stop
            end if
            L = opt_sidelength
            call load_density_perturbation_from_file(trim(filename),delta_r)
            ncells = size(delta_r,1)
        else
           call make_density_perturbation_field_from_particle_file(trim(filename),opt_input,opt_species,-real(opt_ncells),&
           &opt_interpolation,delta_r,opt_sidelength,opt_probability,opt_seed,nparticles=nparticles)
           L = opt_sidelength
           ncells = size(delta_r,1)
           if (nparticles<ncells) then
               write(kwarning,'(A,Es9.3)') 'oversampling, power unreliable for k>',pi*nparticles/L
               write(rwarning,'(A,Es9.3)') 'oversampling, correlations unreliable for r<',L/nparticles
           end if 
           if (opt_interpolation==2) then
               write(kwarning,'(A,Es9.3)') 'some power lost by smoothing for k>',min(ncells,int(nparticles,4))/L
           end if
        end if
    end if
    
    ! output density field
    if (opt_compute_d) then
        call tic('SAVE DENSITY FIELD')
        open(1,file=trim(filename)//'_d.bin',action='write',status='replace',form='unformatted',access='stream')
        write(1) real(delta_r,4)
        close(1)
        call toc()
    end if

    ! output phase field
    if (opt_compute_e) then
        call tic('SAVE PHASE FIELD')
        call make_or_load_DFT(delta_r,delta_k)
        call compute_phase_factors(delta_k,epsil_k)
        deallocate(delta_k)
        call IDFT(epsil_k,epsil_r)
        deallocate(epsil_k)
        open(1,file=trim(filename)//'_e.bin',action='write',status='replace',form='unformatted',access='stream')
        write(1) real(epsil_r,4)
        close(1)
        deallocate(epsil_r)
        call toc()
    end if
    
    ! compute power spectrum / 2-pt correlation
    if ((opt_compute_p).or.(opt_compute_x)) then
        call compute_xi2(delta_r,L,scale,xi,scalek,p)
        if (opt_compute_p) then
            open(1,file=trim(filename)//'_p.txt',action='write',status='replace',form='formatted')
            call write_header('power spectrum','deterministic, exact',kwarning)
            write(1,'(A)') 'Wave vector k [1/input length unit], power spectrum p(k)'
            do i = lbound(scalek,1),ubound(scalek,1)
                write(1,'(Es10.3,A,Es10.3)') scalek(i),',',p(i)
            end do
            close(1)
        end if
        if (opt_compute_x) then
            open(1,file=trim(filename)//'_x.txt',action='write',status='replace',form='formatted')
            call write_header('2-point correlation','deterministic, exact',rwarning)
            write(1,'(A)') 'Correlation scale r [input length unit], 2-point correlation xi(r)'
            do i = 1,size(scale,1)
                write(1,'(Es10.3,A,Es10.3)') scale(i),',',xi(i)
            end do
            close(1)
        end if
    end if
    
    ! 3-pt correlation
    if (opt_compute_y) then
        if (opt_accuracy==0) then
            call compute_xi3(delta_r,L,scale,xi)
            write(str_acc,'(F8.1)') 1.0
        else
            call compute_xi3(delta_r,L,scale,xi,opt_accuracy=opt_accuracy)
            write(str_acc,'(F8.1)') opt_accuracy
        end if
        open(1,file=trim(filename)//'_y.txt',action='write',status='replace',form='formatted')
        call write_header('3-point correlation','stochastic, accuracy = '//trim(adjustl(str_acc)),rwarning)
        write(1,'(A)') 'Correlation scale r [input length unit], 3-point correlation xi3(r)'
        do i = 1,size(scale,1)
            write(1,'(Es10.3,A,Es10.3)') scale(i),',',xi(i)
        end do
        close(1)
    end if
    
    ! compute bispectrum
    if (opt_compute_b) then
        if (opt_accuracy==0) then
            call compute_bispectrum(delta_r,L,bispectrum,angle,kleg)
            write(str_acc,'(F8.1)') 1.0
        else
            call compute_bispectrum(delta_r,L,bispectrum,angle,kleg,opt_accuracy=opt_accuracy)
            write(str_acc,'(F8.1)') opt_accuracy
        end if
        open(1,file=trim(filename)//'_b.txt',action='write',status='replace',form='formatted')
        call write_header('bispectrum','stochastic, accuracy = '//trim(adjustl(str_acc)),kwarning)
        write(1,'(A)') 'angle(k,q) [deg], |k|=|q| [1/input length unit], bispectrum'
        do iangle = 0,ubound(angle,1)
            do ikleg = 1,ubound(kleg,1)
                if (bispectrum(ikleg,iangle)>0) then
                    write(1,'(F5.1,A,Es10.3,A,Es10.3)') angle(iangle),',',kleg(ikleg),',',bispectrum(ikleg,iangle)
                end if
            end do
        end do
        close(1)
    end if
            
    ! compute line correlation
    if (opt_compute_l) then
        if (opt_method==0) then
            if (size(delta_r,1)<=128) then
                opt_method = 2
            else
                opt_method = 1
            end if
        end if
        if (opt_method==1) then
            if (opt_accuracy==0) then
                call compute_lic(delta_r,L,scale,value,error)
                write(str_acc,'(F8.1)') 1.0
            else
                call compute_lic(delta_r,L,scale,value,error,opt_accuracy=opt_accuracy)
                write(str_acc,'(F8.1)') opt_accuracy
            end if
            open(1,file=trim(filename)//'_l.txt',action='write',status='replace',form='formatted')
            call write_header('line correlation','stochastic, accuracy = '//trim(adjustl(str_acc)),rwarning)
            write(1,'(A)') 'Correlation length r [input length unit], line correlation l(r), '//&
            &'1-sigma numerical uncertainty of l(r)'
            do i = 1,size(scale,1)
                write(1,'(Es10.3,A,Es10.3,A,Es10.3)') scale(i),',',value(i),',',error(i)
            end do
            close(1)
        else if (opt_method==2) then
            if (opt_accuracy==0) then
                call compute_lic_deterministic(delta_r,L,scale,value)
            else
                call compute_lic_deterministic(delta_r,L,scale,value,opt_accuracy=opt_accuracy)
            end if
            open(1,file=trim(filename)//'_l.txt',action='write',status='replace',form='formatted')
            if ((opt_accuracy==0).or.(opt_accuracy>=1)) then
                call write_header('line correlation','deterministic, exact',rwarning)
            else
                write(str_acc,'(F8.6)') opt_accuracy
                call write_header('line correlation','deterministic, completeness = '//trim(adjustl(str_acc)),rwarning)
            end if
            write(1,'(A)') 'Correlation length r, line correlation l(r)'
            do i = 1,size(scale,1)
                write(1,'(Es10.3,A,Es10.3)') scale(i),',',value(i)
            end do
            close(1)
        end if
    end if
    
    ! compute special line correlation
    if (opt_compute_a) then
        if (opt_accuracy==0) then
            call compute_lic(delta_r,L,scale,value_edgeworth,error_edgeworth,opt_edgeworth=.true.)
            write(str_acc,'(F8.1)') 1.0
        else
            call compute_lic(delta_r,L,scale,value_edgeworth,error_edgeworth,opt_accuracy=opt_accuracy,&
            &opt_edgeworth=.true.)
            write(str_acc,'(F8.1)') opt_accuracy
        end if
        open(1,file=trim(filename)//'_a.txt',action='write',status='replace',form='formatted')
        call write_header('line correlation in Edgeworth expansion','stochastic, accuracy = '//trim(adjustl(str_acc)),rwarning)
        write(1,'(A)') "Correlation length r, line correlation l'(r), 1-sigma numerical uncertainty of l'(r)"
        do i = 1,size(scale,1)
            write(1,'(Es10.3,A,Es10.3,A,Es10.3)') &
            & scale(i),',',value_edgeworth(i),',',error_edgeworth(i)
        end do
        close(1)
    end if
    
contains
    
    character(len=18) function longint2str(i)
    
        implicit none
        integer*8,intent(in) :: i
        character(len=14) :: pre_str
        character(len=18) :: end_str
        integer    :: j,k
        
        end_str = ''
        write(pre_str,'(I14)') i
        k = 18
        do j = 14,1,-1
            if (pre_str(j:j)==' ') exit
            if (mod(j,3) == 0) then
                end_str(k-1:k) = ','//pre_str(j:j)
                k = k-2
            else
                end_str(k:k) = pre_str(j:j)
                k = k-1
            end if
        end do
        longint2str = end_str
        
    end function longint2str
    
    subroutine write_header(fct,method,warning)
    
        implicit none
        character(*),intent(in) :: fct,method,warning
        character(len=5) :: string
        character(len=50) :: interpolation
        
        if (opt_interpolation==1) then
            interpolation = '(nearest neighbor with shift)'
        else if (opt_interpolation==-1) then
            interpolation = '(nearest neighbor without shift)'
        else if (opt_interpolation==2) then
            interpolation = '(top-hat smoothing with shift)'
        else if (opt_interpolation==-2) then
            interpolation = '(top-hat smoothing without shift)'
        else
            write(*,*) 'Error: interpolation method'
            stop
        end if
        write(string,'(I5)') size(delta_r,1)
        write(1,'(A)') 'Function: '//fct
        write(1,'(A)') 'Method: '//method
        write(1,'(A)') 'Number of grid cells: '//trim(adjustl(string))//'^3 '//trim(interpolation)
        write(1,'(A)') 'Warning: '//trim(warning)
        write(1,'(A)') 'Procorr version '//version
        write(1,'(A)') '--------------------------------------------------------------'
        
    end subroutine write_header
            
end program procorr