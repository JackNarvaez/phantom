!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! rndisc moddump routine: adds a warp in the disc/adds magnetic field
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: part, physcon, setdisc
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use physcon, only:pi
 use part,    only:Bxyz,mhd,rhoh,igas
 use setup_params,   only:ihavesetupB
 integer, intent(in)    :: npartoftype(:)
 real,    intent(in)    :: massoftype(:)
 integer, intent(inout) :: npart
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: npart_start_count,npart_tot,ii,i
 real    :: beta,Bzero,pmassii,phi
 real    :: r2,r,omega,cs,HonR,pressure
 real    :: vphiold2,vphiold,vadd,vphicorr2

 ihavesetupB=.true.
 HonR = 0.05       ! Must check that this is the same as in setup_rndisc

! Similar to that in set_disc
 npart_start_count=1
 npart_tot=npart
!

! Add magnetic field
 if (mhd) then
    beta=25.

! Set up a magnetic field just in Bphi
    do ii = 1,npart
       r2 = xyzh(1,ii)**2 + xyzh(2,ii)**2 + xyzh(3,ii)**2
       r = sqrt(r2)
       phi = atan2(xyzh(2,ii),xyzh(1,ii))
       omega = r**(-1.5)
       cs = HonR*r*omega
       pmassii = massoftype(igas)
       pressure = cs**2*rhoh(xyzh(4,ii),pmassii)
       Bzero = sqrt(2.*pressure/beta)
       Bxyz(1,ii) = -Bzero*sin(phi)
       Bxyz(2,ii) = Bzero*cos(phi)

       ! Calculate correction in v_phi due to B
       vphiold = (-xyzh(2,ii)*vxyzu(1,ii) + xyzh(1,ii)*vxyzu(2,ii))/r
       vphiold2 = vphiold**2
       vphicorr2 = -2.*cs**2
!    if (vphicorr2 > vphi
       vadd = sqrt(vphiold2 + vphicorr2)
       vxyzu(1,ii) = vxyzu(1,ii) + sin(phi)*(vphiold - vadd)
       vxyzu(2,ii) = vxyzu(2,ii) - cos(phi)*(vphiold - vadd)
    
    enddo

    Bxyz(3,:) = 0.0

    print*,'Toroidal Magnetic field added.'
 endif

end subroutine modify_dump

end module moddump
