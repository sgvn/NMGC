module photorates


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author
!> Sacha Gavino & Valentine Wakelam
!
!> @date 2019
!
! DESCRIPTION: all the subroutines to computes the new photorate
!              calculation.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


contains


subroutine compute_molopacity
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! description: computes molecular opacities due to main species
!              contributing to absorption. Read tables composed
!              of absorption, dissociation and ionization cross
!              sections. And uses molecular densities tables. 
!              --> EQUATION (18) in documentation. VERTICAL
!    
! parameters: -mol_opacity   is \tau_m^V(\nu,r,z)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    use utilities
    use global_variables
    implicit none
 

    mol_opacity(:) = 0.D0

    !!!!!! Computation of the opacity due to molecules H2, CO and N2
    ! NH2, NCO, NN2 are global variables and the calculation is done at each altitude so that the column density is computed
    ! from the top to the bottom (most dense part).

! write(*,*) NH2,NCO,NN2

    mol_opacity(:) = NH2*(cross_sections(:,id_photo_H2,2) + cross_sections(:,id_photo_H2,3) + cross_sections(:,id_photo_H2,4))
    mol_opacity(:) = mol_opacity(:) + NCO*(cross_sections(:,id_photo_CO,2) + cross_sections(:,id_photo_CO,3) + &
                    cross_sections(:,id_photo_CO,4))
    mol_opacity(:) = mol_opacity(:) + NN2*(cross_sections(:,id_photo_N2,2) + cross_sections(:,id_photo_N2,3) + &
                    cross_sections(:,id_photo_N2,4))

 !   write(*,*) 'cross_H2:', cross_H2(1,2)
 !   write(*,*) 'cross_CO:', cross_CO(1,2)
 !   write(*,*) 'cross_N2:', cross_N2(1,2)
  
end subroutine compute_molopacity


subroutine compute_dustopacity()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! description: computes dust opacities due to all grain species
!              contributing to absorption. Read tables composed
!              of dust densities along the altitude. VERTICAL
!              --> EQUATION (24) in documentation.
!
! parameters: -dust_opacity   is \tau_d^V(\nu,r,z)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use utilities
    use global_variables
    implicit none
    
    
    !----------------------- VARIABLES -----------------------
    integer :: ios, i, j, k, l
    integer :: nb_param, nb_sizes, size_1D, nb_wv_comments!, spatial_nb=20
    integer :: read_param=9, read_grain_sizes=10, read_static=15, read_wv=16
    
    real(8) :: autocm = 1.49597e13 ! conversion factor from AU to cm.
    real(8) :: z, n_H, T, Av, diff_coef, Td, ab_inv, Av_conv, size
    real(8) :: dust_opacity_line, dz, z_table(2)
    real(8) :: wavelength_nm, wavelength_cm, wavelength_c, Q, Qc = 8 ! this Qc must be changed accordingly to the model, so we must find a way for this routine to get this value somewhere.
    real(8), allocatable :: grain_line(:)
    
    character(len=50) :: density_file, grain_size_file, param_file, static_file, wv_file
    character(len=80) :: param_line, wv_line
    character(len=250) :: static_comments, wv_comments, grain_comments
       
    density_file='dust_density.dat'
    grain_size_file='1D_grain_sizes.in'
    static_file='1D_static.dat'
    param_file='parameters.in'
    wv_file='flux/cut_standard_DRAINE_ISRF_1978.txt'
    !---------------------------------------------------------
    
    
    !------OPEN parameters.in here------
    open(unit=read_param, file=param_file, status="old", action="read", iostat=ios)
    if (ios /= 0) stop "Error opening file parameters.in. Check if parameters.in is there."
    
    !------OPEN 1D_grain_sizes.in here------
    open(unit=read_grain_sizes, file=grain_size_file, status="old", action="read", iostat=ios)
    if (ios /= 0) stop "Error opening file 1D_grain_sizes.in. Check if 1D_grain_sizes.in is there."
    
    !------OPEN 1D_static.dat here------
    open(unit=read_static, file=static_file, status="old", action="read", iostat=ios)
    if (ios /= 0) stop "Error opening file 1D_static.dat. Check if 1D_static.dat is there."
    
    !------OPEN cut_standard_DRAIN_ISRF_1978.txt here------
    open(unit=read_wv, file=wv_file, status="old", action="read", iostat=ios)
    if (ios /= 0) stop "Error opening file cut_standard_DRAINE_ISRF_1978.txt. Check if cut_standard_DRAINE_ISRF_1978.txt is there."
    
    
    !------GET number of lines in parameters.in.------
    nb_param = 0 
    DO 
        READ (read_param,*, END=11) 
        nb_param = nb_param + 1 
    END DO 
    11 CLOSE (1)
  !  print *, "number of lines in parameters.in:", nb_param
    !-------------------------------------------------
    
    
    !------READ parameters.in in order to get the number of grain sizes.------
    REWIND(read_param)
    
    do i = 1,nb_param
        read(read_param,'(A)') param_line
        if (param_line(1:14) == "nb_grains_1D =") then
            read(param_line(15:17),*) nb_sizes
        end if
    enddo
   ! write(*,*) "number of grain sizes parameters.in:", nb_sizes
    !-------------------------------------------------------------------------
    
    
    !------READ 1D_static.in to get dz.------
    do i = 1,16 ! we can take random lines since the difference is constant along the altitude.
        read(read_static,'(a80)') static_comments
    enddo
    
    do i = 1,2
        read(read_static, *) z, n_H, T, Av, diff_coef, Td, ab_inv, Av_conv, size
        z_table(i) = z
    end do 
    
    dz = (z_table(1) - z_table(2))*autocm
    ! "difference between two adajcent altitudes dz (cm):", dz
    !----------------------------------------
    
    
    !------GET number of comments from cut_standard_DRAINE_ISRF_1978.txt.------
    nb_wv_comments = 0 
    DO 
        READ (read_wv,*, END=12) wv_line
        if (wv_line == "#") then
            nb_wv_comments = nb_wv_comments + 1
        end if 
    END DO 
    12 CLOSE (1)
   ! print *, "number of computed wavelength:", nb_line_table_flux
   ! print *, "number of comments in cut_standard_DRAINE_ISRF_1978.txt:", nb_wv_comments
    !----------------------------------------------------------------------------
    
    !------ALLOCATE the size of dust_opacity------
    !allocate(dust_opacity(nb_line_table_flux)) ! This is important. Each line will be the opacity for a wavelength.
    !---------------------------------------------
    
    !---------------------------------------
    !      COMPUTATION OF DUST OPACITY     !
    !---------------------------------------
    
    !------LOOP on wavelength from cut_standard_DRAINE_ISRF_1978.txt.------
    REWIND(read_wv)
    
    IF (nb_sizes.eq.1) then
        write(*,*) 'SINGLE GRAIN'
        do i = 1,nb_wv_comments 
            read(read_wv,'(a80)') wv_comments
        end do
    
        DO i = 1,nb_line_table_flux ! we loop on the wavelengths here. So that we get one opacity per wavelength in the dust_opacity list.
            READ (read_wv,*) wavelength_nm  ! The value here is in nanometer. We convert it in cm afterwards. 
            wavelength_cm = wavelength_nm*1e-7
            REWIND(read_static)
            dust_opacity_line = 0.D0
            do j = 1,12 
                read(read_static,'(a80)') static_comments
            enddo
    
            do j = 1,x_i
                read(read_static, *) z, n_H, T, Av, diff_coef, Td, ab_inv, Av_conv, size
    
                !-----COMPUTE Q-----
                wavelength_c = 2*PI*size
    
                if (wavelength_cm.le.PI*size) then
                    Q = 1
                else if (PI*size.lt.wavelength_cm  .and. wavelength_cm.lt.2*PI*size) then
                    Q = Qc*((wavelength_cm)/wavelength_c)**(log10(Qc)/log10(2.))
                else if (wavelength_cm.ge.2*PI*size) then
                    Q = Qc*((wavelength_cm)/wavelength_c)**(-2)
                end if
                !-----END compute Q----
            
                ! dust_opacity_line = dust_opacity_line + (n_H/ab_inv)*Q*dz*PI*size**2
                dust_opacity_line = dust_opacity_line + (n_H/ab_inv)*Q*dz*PI*size**2


            end do 
  
            dust_opacity(i) = dust_opacity_line

            !write(*,100) "wavelength(nm): ", wavelength_nm, " opacity: ", dust_opacity(i)
            !100 FORMAT(A,ES12.5, A, ES12.5)
        END DO
    
    
    !----------ELSE-----------
    ELSE
        write(*,*) 'MULTIGRAIN. --> number of sizes:', nb_sizes
        size_1D = 2*nb_sizes ! see comment right below

        !------ALLOCATE the size of grain_size_line, i.e. 2xnb_sizes------
        allocate(grain_line(size_1D)) ! This is mportant. x2 because we store the sizes and the densities in this list.
    
    
        do i = 1,nb_wv_comments ! we skip comments
            read(read_wv,'(a80)') wv_comments
        end do
    
        DO i = 1,nb_line_table_flux ! we loop on the wavelengths here. So that we get one opacity per wavelength in the dust_opacity list.
            READ (read_wv,*) wavelength_nm ! The value here is in nanometer. We convert it in cm afterwards. 
            wavelength_cm = wavelength_nm*1e-7

            REWIND(read_static)
            REWIND(read_grain_sizes)
            dust_opacity_line = 0.D0
            do j = 1,12 
                read(read_static,'(a80)') static_comments
            enddo
    
            do k = 1, 2
                read(read_grain_sizes, '(a80)') grain_comments
            end do
            
    
            do k = 1, x_i
                read(read_static, *) z, n_H, T, Av, diff_coef, Td, ab_inv, Av_conv, size
                read(read_grain_sizes,*) grain_line(:)
                do l = 1, nb_sizes
    
                    !-----COMPUTE Q-----
                    wavelength_c = 2*PI*grain_line(l)
    
                    if (wavelength_cm.le.PI*grain_line(l)) then
                        Q = 1
                    else if (PI*grain_line(l).lt.wavelength_cm  .and. wavelength_cm.lt.2*PI*grain_line(l)) then
                        Q = Qc*((wavelength_cm)/wavelength_c)**(log10(Qc)/log10(2.))
                    else if (wavelength_cm.ge.2*PI*grain_line(l)) then
                        Q = Qc*((wavelength_cm)/wavelength_c)**(-2)
                    end if
                    !-----END compute Q----
        
                    dust_opacity_line = dust_opacity_line + (n_H/grain_line(l+nb_sizes))*Q*dz*PI*grain_line(l)**2
                    
                end do ! k
            end do ! l 
    
            dust_opacity(i) =  dust_opacity_line
            !write(*,100) "wavelength(nm): ", wavelength_nm, " opacity: ", dust_opacity(i)
            !100 FORMAT(A,ES12.5, A, ES12.5)
        END DO ! loop on wavelength
    
    END IF ! end of condition on number of sizes. 
    
    !write(*,103) "dust opacity", dust_opacity(31)
    !print *, "dz: ", dz
    !102 FORMAT(ES12.5)
    !103 FORMAT(A,(ES12.5))
    
    !-----CLOSE open files-----
    close(read_param)
    close(read_grain_sizes)
    close(read_static)
    close(read_wv)
 !   write(*,*) 'x_i',x_i
end subroutine compute_dustopacity


subroutine compute_local_flux()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! description: exctract the ISRF and stellar flux spectr from
!              tables then use them to compute the local flux 
!              at each coordinates.
!              --> EQUATION (14) in documentation.
!
! parameters: -Table_local_dust   is I_d(\nu,r,z)
!             -table_ISRF         is I_{ISRF}(\nu)
!             -dust_opacity       is \tau_d^V(\nu,r,z)
!             -Table_Istar        is I_*(\nu,r,z) or I_*^0 (R_0^2/(r^2+z^2))
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use utilities
    use global_variables
    implicit none


! creer une autre subroutine pour lire les differents fichiers et flux (uy compris 
!  table_ISRF) en dehors du pas de temps OK
! table_UV_local est une table donnee par Julia et qui sera lue dans un fichier dependant 
! de la position, independante du temps, dimension(309,64)
! dust_opacity doit etre calculee OK (subroutine compute_dustopacity)
! table_Istar a lire aussi OK
! 

    real(8), dimension(nb_line_table_flux) :: RV


    RV(:) = 0.D0
    RV(:) =  -log ( ( local_flux_dust(:,x_i) - table_ISRF(:,2)*exp(-dust_opacity(:)) ) / table_Istar(:,2) ) &
    / dust_opacity(:) ! --> EQUATION (16)

    !write(*,*) RV(:)
    !print *, "\n"

    local_flux_2(:) = table_ISRF(:,2)*exp(-(dust_opacity(:)+mol_opacity(:))) + &
        &   table_Istar(:,2) * exp(-RV*(dust_opacity(:)+mol_opacity(:)))  ! --> EQUATION (14)

    local_flux(:) = local_flux_dust(:,x_i)*exp(-(mol_opacity(:)))
    !stellar_flux(:) = local_flux_dust(:,1)

    !UV_FLUX =  sum(local_flux_dust(1:wv_max_factor,1))/sum(flux_isrf_dust(1:wv_max_factor,1)) ! UV factor for the original rates calculation.
    INT_LOCAL_FLUX =  sum(local_flux(1:wv_max_factor))/sum(flux_isrf_dust(1:wv_max_factor,1)) ! UV factor for the original rates calculation.

    !write(*,*) 'INT_LOCAL_FLUX: ', INT_LOCAL_FLUX
    !print *, "flux 1"
    !write(*,34) local_flux(:)
    !print *, "\n"

    !if (x_i == 64) then   ! to print only the flux in the midplane.
        !!write(*,*) "FLUX at: ", x_i
        !!write(*,34) local_flux(:)
        !!print *, "\n"
    !end if
        !!34 FORMAT(ES12.5)

end subroutine compute_local_flux



subroutine compute_photorates()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! description: compute the rates (diss and ion) using the 
!              values of the flux and cross sections. 
!              --> EQUATION (6) in documentation.
!
! parameters: -molecule_type    is something
!             -cross-section    is sigma_diss(X,\nu)
!             -local_flux       is I_L(\nu,r,z)
!             -rate_diss        is K_d(X)   
!             -rate_ion         is K_i(X)   
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    use utilities
    use global_variables
    implicit none

    !----------------------- VARIABLES -----------------------
    integer :: i, j
    integer :: nb_cross_comments
    integer :: read_molecule=9

    real(8) :: d_lambda = 1.D0
    !---------------------------------------------------------


    !------INITIALIZE------    
    nb_cross_comments = 0
    photorates%name = 'empty'
    photorates%k_diss = 0.D0
    photorates%k_ion = 0.D0
    !----------------------

    !------LOOP on name of the molecules in cross-sections folder------
    DO i = 1, nb_line_spec_photorates ! number of molecules with cross-sections
   

        !-----COMPUTE photorates.
        photorates(i)%name = spec_photo(i)
        photorates(i)%k_diss = sum(local_flux(:)*cross_sections(:,i,3))*d_lambda
        photorates(i)%k_ion = sum(local_flux(:)*cross_sections(:,i,4))*d_lambda

        !do j = 1,nb_line_photorates ! number of photo reactions of molecules with cross-sections
           ! if ((reaction_compounds_names_photo(j,1).eq.spec_photo(i)).and.&
                    !(reaction_compounds_names_photo(j,5).eq."e-         ")) then
               ! write(*,*) 'photo_ionization ok', reaction_compounds_names_photo(j,1)
                    !table_photo_final(j) = photorates(i)%k_ion
            !ELSE if ((reaction_compounds_names_photo(j,1).eq.spec_photo(i)).and.&
                   ! (reaction_compounds_names_photo(j,5).ne."e-         "))&
            !then
                !table_photo_final(j) = br_photo(j)*photorates(i)%k_diss
                !write(*,*) 'photodissociation ok', reaction_compounds_names_photo(j,1)
           ! end if
        !end do

    END DO ! end loop on molecule names.
    !-------------------------------------------------
    
    !do i = 1, nb_line_spec_photorates
    !    write(*,100) photorates(i)%name, photorates(i)%k_diss, photorates(i)%k_ion
    !end do
    !100 FORMAT(A15, 2ES12.5)

    !-----CLOSE open files-----
    close(read_molecule)
end subroutine compute_photorates



end module photorates
