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
! :Dependencies: centreofmass, io, part, options, physcon, setdisc, units
!                timestep
!

 use setdisc,   only:set_disc
 use units,     only:set_units
 use part,      only:xyzmh_ptmass,vxyz_ptmass,ihacc,ihsoft,nptmass !,Bxyz,mhd
 use physcon,   only:solarm,au,pi
 use io,        only:master
 use options,   only:alpha, ieos
 use timestep,  only:tmax,dtmax
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

!--set time
 time  = 0.

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
 deltat  = 0.01
 norbits = 1

!--setup equation of state
 ieos   = 3
 gamma = 1.0
 
!--resolution
 npart = 1e5
 npartoftype(1) = npart
 hfact = 1.2

 accr1 = R_in
 alpha = alphaSS

 print*,' Mstar is ', Mstar, ' in code units'
 print*,' Mdisc is ', Mdisc, ' in code units'
 print*,' RdAcc is ', accr1, ' in code units'

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

!--reset centre of mass to the origin
 call set_centreofmass(npart,xyzh,vxyzu)

!--single star
!--Central object represented by a sink at the system origin
 nptmass                      = 1
 xyzmh_ptmass(:,:)            = 0.
 xyzmh_ptmass(1:3,nptmass)    = 0.
 xyzmh_ptmass(4,nptmass)      = Mstar
 xyzmh_ptmass(ihacc,nptmass)  = accr1
 xyzmh_ptmass(ihsoft,nptmass) = 0.
 vxyz_ptmass                  = 0.

!--outer disc orbital period
 period = sqrt(4.*pi**2*R_out**3/Mstar)
 dtmax  = deltat*period
 tmax   = norbits*period

end subroutine setpart

!--------------------------------------------------------------------------
!
!  Reset centre of mass to origin
!
!--------------------------------------------------------------------------
subroutine set_centreofmass(npart,xyzh,vxyzu)
    use centreofmass, only:reset_centreofmass
    integer, intent(in)    :: npart
    real,    intent(inout) :: xyzh(:,:)
    real,    intent(inout) :: vxyzu(:,:)
   
    call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)
   
   end subroutine set_centreofmass

end module setup