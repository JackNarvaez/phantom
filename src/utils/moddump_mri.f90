!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Adds magnetic field
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: part, setup_params, timestep
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,         only:Bxyz,mhd,rhoh,igas
 use setup_params, only:ihavesetupB
 use timestep,     only:overcleanfac
 use eos,          only:qfacdisc
  use io,          only:fatal
 integer, intent(in)    :: npartoftype(:)
 real,    intent(in)    :: massoftype(:)
 integer, intent(inout) :: npart
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: npart_start_count,npart_tot,igeom,i
 real    :: Bzero,pmassii,phi
 real    :: pindex,qindex,cs0,betaP
 real    :: r,r2,cs02,cs2,pressure
 real    :: vphiold2,vphiold,vadd,vphiadd2,corrf
 logical :: reverse_field_dir

 ! Must check that this is the same as in setup_disc
 cs0 = 0.05
 cs02 = cs0*cs0

 qindex = qfacdisc
 pindex = 1

 npart_start_count=1
 npart_tot=npart

!--add magnetic field
!--only 'toroidal' and 'vertical' geometries supported
!--field set using constant plasma beta and isothermal pressure
 if (mhd) then
    ihavesetupB = .true.
    overcleanfac  = 1.0
    igeom = 1  ! 1: toroidal; 2: vertical
    betaP = 25.
    reverse_field_dir = .false. ! if true: clockwise

    select case(igeom)
    case(1)
       corrf = -(pindex+qindex-0.5)/betaP
    case(2)
       corrf = -(1.5+pindex+qindex)/betaP
    case default
       call fatal('set_Bfield','unknown field geometry')
    end select

! Set up a magnetic field just in Bphi
    do i = npart_start_count,npart_tot
       r2 = xyzh(1,i)**2 + xyzh(2,i)**2
       r = sqrt(r2)
       phi = atan2(xyzh(2,i),xyzh(1,i))
       cs2 = cs02*r2**(-qindex)
       pmassii = massoftype(igas)
       pressure = cs2*rhoh(xyzh(4,i),pmassii)
       Bzero = sqrt(2.*pressure/betaP)

       if (reverse_field_dir) Bzero = -Bzero

       select case(igeom)
       ! toroidal magnetic field (Bphi)
       case(1)
        Bxyz(1,i) = -Bzero*sin(phi)
        Bxyz(2,i) = Bzero*cos(phi)
        Bxyz(3,i) = 0.0d0

       ! vertical magnetic field (Bz)
       case(2)
        Bxyz(1,i) = 0.0d0
        Bxyz(2,i) = 0.0d0
        Bxyz(3,i) = Bzero

       case default
        call fatal('set_Bfield','unknown field geometry')
       end select

       ! Calculate correction in v_phi due to B
        vphiold = (-xyzh(2,i)*vxyzu(1,i) + xyzh(1,i)*vxyzu(2,i))/r
        vphiold2 = vphiold*vphiold
        vphiadd2 = vphiold2 + corrf*cs2
        if (vphiadd2<0) vphiadd2 = vphiold2
        vadd = sqrt(vphiadd2)
        vxyzu(1,i) = vxyzu(1,i) + sin(phi)*(vphiold - vadd)
        vxyzu(2,i) = vxyzu(2,i) - cos(phi)*(vphiold - vadd)
    enddo

    print*,'Magnetic field added.'

    print*, ""
    print*,'|------------ PARAMETERS ------------|'
    print*, ""

    print '(A,F12.4)',' pindex     = ', pindex
    print '(A,F12.4)',' qindex     = ', qindex
    print '(A,F12.4)',' cs0        = ', cs0

    print*,' ---------------- MHD --------------- '
    print '(A,L5)'   ,' mhd        = ', mhd
    print '(A,I12)'  ,' geometry   = ', igeom
    print '(A,F12.4)',' beta_mag   = ', betaP
    print '(A,L5)'   ,' orientaton = ', reverse_field_dir

    print*, ""
    print*,'|---------- END SETUP FILE ----------|'
    print*, ""
 endif

end subroutine modify_dump

end module moddump
