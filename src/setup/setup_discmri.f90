!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! This module sets up accretion discs. The goal is to simulate 
! the magnetorotational instability for the global disc. 
!
! :References: None
!
! :Owner: Jacksen Narvaez
!
! :Runtime parameters: None
!
! :Dependencies: io, part, options, physcon, setdisc, units
!                timestep
!
 implicit none
 public :: setpart

 private

contains

!----------------------------------------------------------------
!
! This subroutine is a utility for setting up discs
!
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setdisc,   only:set_disc
 use units,     only:set_units
 use part,      only:xyzmh_ptmass,vxyz_ptmass,ihacc,ihsoft,nptmass !,Bxyz,mhd
 use physcon,   only:solarm,au,pi
 use io,        only:master
 use options,   only:alpha, ieos
 use timestep,  only:tmax,dtmax

 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real :: R_in, R_out, R_ref
 real :: pindex, qindex, H_R
 real :: alphaSS
 real :: posangl, incl
 real :: Mdisc,Mstar
 real :: sig_ref, Q_min
 real :: period, deltat
 real :: accr1
 integer :: norbits
 integer :: icentral
 integer :: nsinks
 logical :: ismoothgas

! set code units
 call set_units(dist=au,mass=solarm,G=1.d0)

!--central object(s)
 icentral   = 1

!--sink particle(s)
 nsinks     = 1

!--gas disc
 R_in       = 1.
 R_out      = 150.
 R_ref      = 10.
 pindex     = 1.
 qindex     = 0.25
 alphaSS    = 0.005
 posangl    = 0.
 incl       = 0.
 H_R        = 0.05
 Mstar      = 1.0
 Mdisc      = 0.05
 sig_ref    = 1.e-02
 Q_min      = 1.0
 ismoothgas = .true.

 !--simulation time
 deltat  = 0.1
 norbits = 1

 gamma = 1.0
 npart = 1e5
 npartoftype(1) = npart
 hfact = 1.2
 time  = 0.
 accr1 = R_in

 print*,' Mstar is ', Mstar, ' in code units'
 print*,' Mdisc is ', Mdisc, ' in code units'

!--setup disc(s)
 call set_disc(id,master        = master,           &
               npart            = npartoftype(1),   &
               rref             = R_ref,            &
               rmin             = R_in,             &
               rmax             = R_out,            &
               p_index          = pindex,           &
               q_index          = qindex,           &
               HoverR           = H_R,              &
               disc_mass        = Mdisc,            &
               star_mass        = Mstar,            &
               gamma            = gamma,            &
               particle_mass    = massoftype(1),    &
               ismooth          = ismoothgas,       &
               hfact            = hfact,            &
               xyzh             = xyzh,             &
               vxyzu            = vxyzu,            &
               polyk            = polyk,            &
               alpha            = alpha,            &
               prefix           = fileprefix)

!--single star
 print "(/,a)",' Central object represented by a sink at the system origin'
 print "(a,g10.3)",'   Object mass:      ', Mstar
 print "(a,g10.3)",'   Accretion Radius: ', accr1
 nptmass                      = 1
 xyzmh_ptmass(:,:)            = 0.
 xyzmh_ptmass(1:3,nptmass)    = 0.
 xyzmh_ptmass(4,nptmass)      = Mstar
 xyzmh_ptmass(ihacc,nptmass)  = accr1
 xyzmh_ptmass(ihsoft,nptmass) = 0.
 vxyz_ptmass                  = 0.

!--setup equation of state
 ieos   = 3

!--outer disc orbital period
 period = sqrt(4.*pi**2*R_out**3/Mstar)
 dtmax  = deltat*period
 tmax   = norbits*period



end subroutine setpart

end module setup
