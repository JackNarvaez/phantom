!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! This module sets up accretion discs. The goal is to simulate 
! the magnetorotational instability (MRI) for the global disc. 
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
 
 use io,                only:master
 use eos,               only:qfacdisc
 use options,           only:alpha,ieos,nfulldump,overcleanfac
 use part,              only:xyzmh_ptmass,vxyz_ptmass,ihacc,ihsoft,nptmass,Bxyz,mhd,rhoh,igas
 use physcon,           only:solarm,au,pi 
 use setdisc,           only:set_disc
 use setup_params,      only:ihavesetupB
 use timestep,          only:tmax,dtmax
 use units,             only:set_units

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
 real :: R_in,R_out,R_ref
 real :: pindex,qindex,H_R
 real :: alphaSS
 real :: posangl,incl
 real :: Mdisc,Mstar
 real :: period,deltat
 real :: accr1
 real :: pmassi
 real    :: beta,Bzero,phi
 real    :: R2,R,vkep2,cs2,pressure
 real    :: vnew
 integer :: norbits
 integer :: icentral
 integer :: nsinks
 integer :: i
 logical :: ismoothgas

! set code units
 call set_units(dist=au,mass=solarm,G=1.d0)

!--set time
 time  = 0.0

!--central object(s)
 icentral   = 1
 nsinks     = 1

!--gas disc
 R_in       = 1.
 R_out      = 10.
 R_ref      = 1.
 pindex     = 1.
 qindex     = 0.25
 alphaSS    = 0.005
 posangl    = 0.
 incl       = 0.
 H_R        = 0.05
 Mstar      = 1.0
 Mdisc      = 0.05
 ismoothgas = .true.

!--simulation time
 deltat     = 0.01
 norbits    = 10
 nfulldump  = 1

!--setup equation of state
 ieos     = 3
 gamma    = 1.0
 qfacdisc = qindex

!--resolution
 npart = 1e6
 npartoftype(igas) = npart
 hfact = 1.2

 accr1 = R_in
 alpha = alphaSS

 print '(A,F12.4,A)',' Mstar is ', Mstar, ' in code units'
 print '(A,F12.4,A)',' Mdisc is ', Mdisc, ' in code units'
 print '(A,F12.4,A)',' RdAcc is ', accr1, ' in code units'

!--setup disc(s)
 call set_disc(id,master        = master,               &
               npart            = npartoftype(igas),    &
               rref             = R_ref,                &
               rmin             = R_in,                 &
               rmax             = R_out,                &
               p_index          = pindex,               &
               q_index          = qindex,               &
               HoverR           = H_R,                  &
               disc_mass        = Mdisc,                &
               star_mass        = Mstar,                &
               gamma            = gamma,                &
               particle_mass    = massoftype(igas),     &
               ismooth          = ismoothgas,           &
               hfact            = hfact,                &
               xyzh             = xyzh,                 &
               vxyzu            = vxyzu,                &
               polyk            = polyk,                &
               alpha            = alpha,                &
               prefix           = fileprefix)

!--central object represented by a sink at the system origin
 nptmass                      = nsinks
 xyzmh_ptmass(:,:)            = 0.
 xyzmh_ptmass(1:3,nptmass)    = 0.
 xyzmh_ptmass(4,nptmass)      = Mstar
 xyzmh_ptmass(ihacc,nptmass)  = accr1
 xyzmh_ptmass(ihsoft,nptmass) = 0.01*accr1
 vxyz_ptmass                  = 0.

!--outer disc orbital period
 period = sqrt(4.*pi**2*R_out**3/Mstar)
 dtmax  = deltat*period
 tmax   = norbits*period

!--add magnetic field
!--set magnetic field using plasma beta 
 if (mhd) then
  ihavesetupB   = .true.
  overcleanfac  = 1.0
  beta = 1000.

  ! toroidal field
  ! set up a magnetic field just in Bphi
  do i = 1,npart
   R2        = xyzh(1,i)**2 + xyzh(2,i)**2
   R         = sqrt(R2)
   phi       = atan2(xyzh(2,i),xyzh(1,i))
   vkep2     = Mstar/R
   cs2       = polyk*R2**(-qindex)
   pmassi    = massoftype(igas)
   pressure  = cs2*rhoh(xyzh(4,i),pmassi)
   Bzero     = sqrt(2.*pressure/beta)
   Bxyz(1,i) = -Bzero*sin(phi)
   Bxyz(2,i) = Bzero*cos(phi)

   ! calculate correction in v_phi due to B
   vnew = vkep2-cs2*(1.5+pindex+qindex)*(1.+1./beta)
   vxyzu(1,i) = -sqrt(vnew)*sin(phi)
   vxyzu(2,i) =  sqrt(vnew)*cos(phi)
   vxyzu(3,i) = 0.0d0
  enddo
  Bxyz(3,:) = 0.0d0
 endif

!--reset centre of mass to the origin
 call set_centreofmass(npart,xyzh,vxyzu)

! call set_centreofmass(npart,xyzh,vxyzu)

 print*, ""
 print*,'|----------- END SETUP FILE ---------|'
 print*, ""

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
