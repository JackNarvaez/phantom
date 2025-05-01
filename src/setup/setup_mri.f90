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
 
 use dim,               only:maxalpha,maxp,nalpha
 use io,                only:master
 use eos,               only:qfacdisc
 use options,           only:alpha,alphamax,ieos,nfulldump,overcleanfac
 use part,              only:xyzmh_ptmass,vxyz_ptmass,ihacc,ihsoft,nptmass,Bxyz,mhd,rhoh,igas
 use physcon,           only:solarm,au,pi 
 use setdisc,           only:set_disc
 use setup_params,      only:ihavesetupB
 use timestep,          only:tmax,dtmax
 use units,             only:set_units
 use viscosity,         only:irealvisc,shearparam,bulkvisc

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
 real :: pindex,qindex,q_z,H_R
 real :: alphaSS,alphaMX
 real :: posangl,incl
 real :: Mdisc,Mstar
 real :: period,deltat
 real :: accr1
 real :: pmassi
 real    :: beta,Bzero,phi
 real    :: R2,R,z2,vkep2,cs2,pressure
 real    :: vnew2
 integer :: norbits
 integer :: icentral
 integer :: nsinks
 integer :: i
 integer :: visc
 logical :: ismoothgas,shearz
 logical :: nsvisc
 character(len=16) :: geometry

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
 alphaMX    = 0.5
 posangl    = 0.
 incl       = 0.
 H_R        = 0.05
 Mstar      = 1.0
 Mdisc      = 0.05
 ismoothgas = .true.
 shearz     = .true.
 nsvisc     = .true.  ! SS viscosity

!--simulation time
 deltat     = 0.1
 norbits    = 100
 nfulldump  = 1

!--setup equation of state
 ieos     = 3
 gamma    = 1.0
 qfacdisc = qindex

!--resolution
 npart = 1e6
 npartoftype(igas) = npart
 hfact = 1.2

!--accretion radius
 accr1 = R_in

!--viscosity
! Disc viscosity LP10. Switches are turned off.
 if (maxalpha==0) then
  alpha = alphaSS
  visc  = 3
 elseif (maxalpha==maxp) then
  alphamax = alphaMX
  !alpha = 0.0
! Add Shakura/Sunyaev disc viscosity via real Navier-Stokes viscosity  
  if (nsvisc) then
   irealvisc = 2
   shearparam = alphaSS
   bulkvisc = 0
   visc = 4
! No Disc viscosity
  else
   visc = 0
   irealvisc = 0
  endif
  ! Artificial viscosity switches
  if (nalpha >= 2) then
  ! Cullen-Dehnen switch
    if (visc == 0) then
     visc = 1
    endif
  else
  ! Morris-Monaghan switch
    if (visc == 0) then
     visc = 2
    endif
  endif
 endif

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

!--add vertical shear
 if (shearz) then
  q_z = qindex
 else
  q_z = 0.0
 endif

!--add magnetic field
!--only 'toroidal' and 'vertical' geometries supported
!--field set using constant plasma beta and isothermal pressure
 if (mhd) then
  ihavesetupB   = .true.
  overcleanfac  = 2.0
  geometry = 'vertical'
  beta = 1000.

  do i = 1,npart
   R2        = xyzh(1,i)**2 + xyzh(2,i)**2
   z2        = xyzh(3,i)**2
   R         = sqrt(R2)
   phi       = atan2(xyzh(2,i),xyzh(1,i))
   vkep2     = Mstar/R
   cs2       = polyk*R2**(-qindex)
   pmassi    = massoftype(igas)
   pressure  = cs2*rhoh(xyzh(4,i),pmassi)
   Bzero     = sqrt(2.*pressure/beta)
  
  ! toroidal magnetic field (Bphi)
  if (geometry == 'toroidal') then
   Bxyz(1,i) = -Bzero*sin(phi)
   Bxyz(2,i) = Bzero*cos(phi)
   Bxyz(3,i) = 0.0d0
  ! vertical magnetic field (Bz)
  elseif (geometry == 'vertical') then
   Bxyz(1,i) = 0.0d0
   Bxyz(2,i) = 0.0d0
   Bxyz(3,i) = Bzero
  else
   print *, 'Error: Unknown magnetic field geometry: ', trim(geometry)
   stop
  endif

   ! calculate correction in v_phi due to B
   vnew2      = vkep2-(cs2*(1.5+pindex+qindex)+q_z*vkep2*z2/R2)*(1.+1./beta);
   vxyzu(1,i) = -sqrt(vnew2)*sin(phi)
   vxyzu(2,i) =  sqrt(vnew2)*cos(phi)
   vxyzu(3,i) = 0.0d0
  enddo
!--if vertical shear is on, then add dependency on z for ang.vel
 elseif (shearz) then
  do i=1,npart
   R2        = xyzh(1,i)**2 + xyzh(2,i)**2
   z2        = xyzh(3,i)**2
   R         = sqrt(R2)
   phi       = atan2(xyzh(2,i),xyzh(1,i))
   vkep2     = Mstar/R
   cs2       = polyk*R2**(-qindex)

   ! calculate correction in v_phi due to B
   vnew2      = vkep2-(cs2*(1.5+pindex+qindex)+q_z*vkep2*z2/R2);
   vxyzu(1,i) = -sqrt(vnew2)*sin(phi)
   vxyzu(2,i) =  sqrt(vnew2)*cos(phi)
   vxyzu(3,i) = 0.0d0
  enddo
 endif

!--reset centre of mass to the origin
 call set_centreofmass(npart,xyzh,vxyzu)

 print*, ""
 print*,'|------------ PARAMETERS ------------|'
 print*, ""

 print '(A,F12.4)',' Mstar      = ', Mstar
 print '(A,F12.4)',' Mdisc      = ', Mdisc
 print '(A,F12.4)',' RdAcc      = ', accr1
 print*,' ---------------- MHD --------------- '
 print '(A,L5)'   ,' mhd        = ', mhd
 print '(A,A)'    ,' geometry   = ', trim(geometry)
 print '(A,F12.4)',' beta_mag   = ', beta
 print*,' ---------------- VISC --------------- '
 print '(A,F12.4)',' alpha      = ', alpha
 print '(A,I12)'  ,' irealvisc  = ', irealvisc
 print '(A,F12.4)',' shearparam = ', shearparam
 print '(A,F12.4)',' bulkvisc   = ', bulkvisc
 print '(A,I12,A)',' visc       = ', visc        ,' (1: AVC; 2: AVM; 3: DVA; 4:DVN)'
 print '(A,I12)'  ,' maxalpha   = ', maxalpha
 print '(A,I12,A)',' nalpha     = ', nalpha      ,' (0: none; 1: Morris-Monaghan; 3: Cullen-Dehnen)'
 print '(A,L5)'   ,' zth-shear  = ', shearz

 print*, ""
 print*,'|---------- END SETUP FILE ----------|'
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
