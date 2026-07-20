!WRF:DRIVER_LAYER:CONFIGURATION
!

MODULE module_scalar_tables
  USE module_driver_constants
  USE module_state_description
  USE module_domain_type, ONLY : streamrec
#include "scalar_tables.inc"


CONTAINS


  SUBROUTINE init_module_scalar_tables
     INTEGER i , j
     DO j = 1, max_domains
#include "scalar_tables_init.inc"
     END DO
  END SUBROUTINE init_module_scalar_tables
END MODULE module_scalar_tables

#ifdef WRF_CHEM
#ifdef WRF_KPP
MODULE module_irr_diag

   INTEGER, parameter :: max_eqn = 1200
   INTEGER :: max_dom
   INTEGER, ALLOCATABLE :: irr_diag_cnt(:)
   INTEGER, ALLOCATABLE :: irr_diag_ndx(:,:)
   LOGICAL, ALLOCATABLE :: irr_option(:)
   LOGICAL, ALLOCATABLE :: irr_active(:)


CONTAINS



  SUBROUTINE init_module_irr_diag
    
    INTEGER :: astat

    IF( .not. ALLOCATED( irr_option ) ) THEN
      CALL nl_get_max_dom( 1, max_dom )
      ALLOCATE( irr_option(max_dom),irr_active(max_dom),irr_diag_ndx(max_eqn,max_dom),irr_diag_cnt(max_dom),stat=astat )
      IF( astat /= 0 ) THEN
        CALL wrf_error_fatal( "init_module_irr_diag: Failed to allocate irr_option ... irr_diag_cnt" )
      ENDIF
      irr_option(:) = .false.
      irr_active(:) = .false.
    ENDIF

  END SUBROUTINE init_module_irr_diag

END MODULE module_irr_diag
#endif
#endif

MODULE module_configure

   USE module_driver_constants
   USE module_state_description
   USE module_wrf_error

#ifdef WRF_CFBM
   ! fire_behavior
   use namelist_mod, only : namelist_t
#endif

   TYPE model_config_rec_type
      SEQUENCE
! Statements that declare namelist variables are in this file
! Note that the namelist is SEQUENCE and generated such that the first item is an
! integer, first_item_in_struct and the last is an integer last_item_in_struct
! this provides a way of converting this to a buffer for passing to and from
! the driver.
#include "namelist_defines.inc"
   END TYPE model_config_rec_type

   TYPE grid_config_rec_type
#include "namelist_defines2.inc"
   END TYPE grid_config_rec_type

   TYPE(model_config_rec_type) :: model_config_rec

!#include "scalar_tables.inc"

! special entries (put here but not enshrined in Registry for one reason or other)

!   CHARACTER (LEN=256) :: mminlu = ' '             ! character string for landuse table
#ifdef WRF_CFBM
   TYPE(namelist_t) :: cfbm_namelist
#endif

   INTEGER, PARAMETER :: MAX_SPRAY_SOURCES = 20

   INTEGER, SAVE :: spray_source_domain(MAX_SPRAY_SOURCES)

   REAL, SAVE :: spray_source_lat(MAX_SPRAY_SOURCES)
   REAL, SAVE :: spray_source_lon(MAX_SPRAY_SOURCES)
   REAL, SAVE :: spray_number_rate(MAX_SPRAY_SOURCES)
   REAL, SAVE :: spray_mode_diameter(MAX_SPRAY_SOURCES)
   REAL, SAVE :: spray_geometric_sigma(MAX_SPRAY_SOURCES)
   REAL, SAVE :: spray_particle_density(MAX_SPRAY_SOURCES)
   INTEGER, SAVE :: spray_source_i(MAX_SPRAY_SOURCES)
   INTEGER, SAVE :: spray_source_j(MAX_SPRAY_SOURCES)
   integer, save :: num_spray_sources1
   LOGICAL, SAVE :: spray_on_the_task(MAX_SPRAY_SOURCES)

CONTAINS

SUBROUTINE locate_spray_sources(                         &
     current_grid, config_flags,                         &
     num_spray_sources,                                  &
     spray_source_lat, spray_source_lon,                 &
     spray_source_domain,                                &
     spray_source_i, spray_source_j,                     &
     spray_on_the_task, ips,ipe,jps,jpe )

   USE module_domain_type, ONLY : domain

   USE module_llxy, ONLY :                               &
        proj_info, map_set, latlon_to_ij,                &
        PROJ_LATLON, PROJ_LC, PROJ_PS, PROJ_MERC

   IMPLICIT NONE

   TYPE(domain), INTENT(IN) :: current_grid
   TYPE(grid_config_rec_type), INTENT(IN) :: config_flags

   INTEGER, INTENT(IN) :: num_spray_sources
   INTEGER, INTENT(IN) :: ips, ipe, jps, jpe
   REAL, INTENT(IN) :: spray_source_lat(:)
   REAL, INTENT(IN) :: spray_source_lon(:)

   INTEGER, INTENT(IN) :: spray_source_domain(:)

   INTEGER, INTENT(INOUT) :: spray_source_i(:)
   INTEGER, INTENT(INOUT) :: spray_source_j(:)

   LOGICAL, INTENT(INOUT) :: spray_on_the_task(:)

   TYPE(proj_info) :: proj

   INTEGER :: n
   INTEGER :: isrc, jsrc
   INTEGER :: i_anchor, j_anchor
   INTEGER :: ips_src, ipe_src
   INTEGER :: jps_src, jpe_src

   REAL :: i_val, j_val
   REAL :: known_lat, known_lon
   REAL :: latinc, loninc

   LOGICAL :: on_patch

   ! Select a valid mass-grid point on this MPI task.
   !
   ! The anchor latitude/longitude and knowni/knownj must refer to
   ! exactly the same global WRF grid point.
   !---------------------------------------------------------------
   WRITE(0,*) '**SPRAY domain**', spray_source_domain
   i_anchor = MAX(current_grid%sp31, current_grid%sd31)
   j_anchor = MAX(current_grid%sp33, current_grid%sd33)

   IF (i_anchor > current_grid%ep31 .OR.                 &
       j_anchor > current_grid%ep33) THEN

      WRITE(0,*) 'SPRAY ERROR: no valid local anchor:',  &
                 ' domain=', current_grid%id,            &
                 ' anchor=', i_anchor, j_anchor,         &
                 ' patch=',                              &
                 current_grid%sp31, current_grid%ep31,   &
                 current_grid%sp33, current_grid%ep33

      CALL wrf_error_fatal(                              &
           'locate_spray_sources: no valid local anchor')

   END IF

   IF (i_anchor < LBOUND(current_grid%xlat,1) .OR.       &
       i_anchor > UBOUND(current_grid%xlat,1) .OR.       &
       j_anchor < LBOUND(current_grid%xlat,2) .OR.       &
       j_anchor > UBOUND(current_grid%xlat,2)) THEN

      WRITE(0,*) 'SPRAY ERROR: anchor outside XLAT:',    &
                 ' domain=', current_grid%id,            &
                 ' anchor=', i_anchor, j_anchor,         &
                 ' XLAT bounds=',                        &
                 LBOUND(current_grid%xlat,1),            &
                 UBOUND(current_grid%xlat,1),            &
                 LBOUND(current_grid%xlat,2),            &
                 UBOUND(current_grid%xlat,2)

      CALL wrf_error_fatal(                              &
           'locate_spray_sources: anchor outside XLAT')

   END IF

   known_lat = current_grid%xlat (i_anchor,j_anchor)
   known_lon = current_grid%xlong(i_anchor,j_anchor)

   IF (known_lat /= known_lat .OR.                       &
       known_lon /= known_lon .OR.                       &
       ABS(known_lat) > 90.0 .OR.                       &
       ABS(known_lon) > 360.0) THEN

      WRITE(0,*) 'SPRAY ERROR: invalid anchor lat/lon:', &
                 ' domain=', current_grid%id,            &
                 ' anchor=', i_anchor, j_anchor,         &
                 ' lat/lon=', known_lat, known_lon

      CALL wrf_error_fatal(                              &
           'locate_spray_sources: invalid anchor lat/lon')

   END IF

   !---------------------------------------------------------------
   ! Report the projection values actually being used.
   ! These must come from config_flags, not current_grid, because
   ! current_grid%truelat1, etc. are zero in this setup.
   !---------------------------------------------------------------

   WRITE(0,*) 'SPRAY PROJECTION SETTINGS:',              &
              ' domain=', current_grid%id,               &
              ' map_proj=', config_flags%map_proj,       &
              ' anchor i/j=', i_anchor, j_anchor,        &
              ' anchor lat/lon=', known_lat, known_lon,  &
              ' dx=', current_grid%dx,                   &
              ' truelat1/2=',                            &
              model_config_rec%truelat1(current_grid%id),                    &
              model_config_rec%truelat2(current_grid%id),                    &
              ' stand_lon=', model_config_rec%stand_lon(current_grid%id)

   !---------------------------------------------------------------
   ! Construct the projection for the current domain.
   !---------------------------------------------------------------

   SELECT CASE (config_flags%map_proj)

   CASE (PROJ_MERC)

      IF (ABS(model_config_rec%truelat1(current_grid%id)) > 90.0) THEN
         CALL wrf_error_fatal(                           &
              'Invalid truelat1 for Mercator spray mapping')
      END IF

      CALL map_set(                                      &
           PROJ_MERC, proj,                              &
           truelat1 = model_config_rec%truelat1(current_grid%id),            &
           lat1     = known_lat,                         &
           lon1     = known_lon,                         &
           knowni   = REAL(i_anchor),                    &
           knownj   = REAL(j_anchor),                    &
           dx       = current_grid%dx)

   CASE (PROJ_LC)

      IF (ABS(model_config_rec%truelat1(current_grid%id)) < 1.0E-6 .AND.      &
          ABS(model_config_rec%truelat2(current_grid%id)) < 1.0E-6) THEN

         WRITE(0,*) 'SPRAY ERROR: invalid Lambert values:', &
                    model_config_rec%truelat1(current_grid%id),               &
                    model_config_rec%truelat2(current_grid%id),               &
                    model_config_rec%stand_lon(current_grid%id)

         CALL wrf_error_fatal(                           &
              'Invalid Lambert parameters for spray mapping')

      END IF

      CALL map_set(                                      &
           PROJ_LC, proj,                                &
           truelat1 = model_config_rec%truelat1(current_grid%id),            &
           truelat2 = model_config_rec%truelat2(current_grid%id),            &
           stdlon   = model_config_rec%stand_lon(current_grid%id),           &
           lat1     = known_lat,                         &
           lon1     = known_lon,                         &
           knowni   = REAL(i_anchor),                    &
           knownj   = REAL(j_anchor),                    &
           dx       = current_grid%dx)

   CASE (PROJ_PS)

      CALL map_set(                                      &
           PROJ_PS, proj,                                &
           truelat1 = model_config_rec%truelat1(current_grid%id),            &
           stdlon   = model_config_rec%stand_lon(current_grid%id),           &
           lat1     = known_lat,                         &
           lon1     = known_lon,                         &
           knowni   = REAL(i_anchor),                    &
           knownj   = REAL(j_anchor),                    &
           dx       = current_grid%dx)

   CASE (PROJ_LATLON)

      ! module_llxy requires angular increments for PROJ_LATLON,
      ! not dx in metres. Derive them from neighbouring mass points.

      IF (i_anchor < UBOUND(current_grid%xlong,1)) THEN

         loninc = current_grid%xlong(i_anchor+1,j_anchor) &
                - current_grid%xlong(i_anchor,j_anchor)

      ELSE IF (i_anchor > LBOUND(current_grid%xlong,1)) THEN

         loninc = current_grid%xlong(i_anchor,j_anchor)   &
                - current_grid%xlong(i_anchor-1,j_anchor)

      ELSE

         CALL wrf_error_fatal(                           &
              'Cannot determine longitude increment')

      END IF

      IF (j_anchor < UBOUND(current_grid%xlat,2)) THEN

         latinc = current_grid%xlat(i_anchor,j_anchor+1)  &
                - current_grid%xlat(i_anchor,j_anchor)

      ELSE IF (j_anchor > LBOUND(current_grid%xlat,2)) THEN

         latinc = current_grid%xlat(i_anchor,j_anchor)    &
                - current_grid%xlat(i_anchor,j_anchor-1)

      ELSE

         CALL wrf_error_fatal(                           &
              'Cannot determine latitude increment')

      END IF

      IF (ABS(latinc) < 1.0E-12 .OR.                    &
          ABS(loninc) < 1.0E-12) THEN

         WRITE(0,*) 'SPRAY ERROR: invalid lat/lon increments:', &
                    latinc, loninc

         CALL wrf_error_fatal(                           &
              'Invalid increments for lat-lon spray mapping')

      END IF

      CALL map_set(                                      &
           PROJ_LATLON, proj,                            &
           lat1     = known_lat,                         &
           lon1     = known_lon,                         &
           knowni   = REAL(i_anchor),                    &
           knownj   = REAL(j_anchor),                    &
           latinc   = latinc,                            &
           loninc   = loninc)

   CASE DEFAULT

      WRITE(0,*) 'SPRAY ERROR: unsupported map projection:', &
                 config_flags%map_proj

      CALL wrf_error_fatal(                              &
           'locate_spray_sources: unsupported projection')

   END SELECT

   !---------------------------------------------------------------
   ! Processor-owned bounds for the unstaggered mass grid.
   !---------------------------------------------------------------

   ips_src = current_grid%sp31
   ipe_src = current_grid%ep31

   jps_src = current_grid%sp33
   jpe_src = current_grid%ep33

   WRITE(0,*) '============================================'
   WRITE(0,*) 'SPRAY DIAGNOSTICS: domain=', current_grid%id
   WRITE(0,*) '============================================'
   WRITE(0,*) 'Anchor point:'
   WRITE(0,*) '  i_anchor, j_anchor =', i_anchor, j_anchor
   WRITE(0,*) '  lat, lon =', known_lat, known_lon
   WRITE(0,*) ''
   WRITE(0,*) 'Domain bounds (global grid extent):'
   WRITE(0,*) '  i: ', current_grid%sd31, ' to ', current_grid%ed31-1
   WRITE(0,*) '  j: ', current_grid%sd33, ' to ', current_grid%ed33-1
   WRITE(0,*) ''
   WRITE(0,*) 'Patch bounds (this processor):'
   WRITE(0,*) '  i: ', ips_src, ' to ', ipe_src
   WRITE(0,*) '  j: ', jps_src, ' to ', jpe_src
   WRITE(0,*) ''
   WRITE(0,*) 'Spray sources to locate:'
   WRITE(0,*) '============================================'

   !---------------------------------------------------------------
   ! Locate only sources assigned to the current domain.
   !---------------------------------------------------------------

   WRITE(0,*) 'DEBUG: locate_spray_sources called for domain', current_grid%id, &
           ' on processor with patch i=', ips_src, ipe_src   !---------------------------------------------------------------

   DO n = 1, num_spray_sources

      IF (spray_source_domain(n) /= current_grid%id) CYCLE

      spray_source_i(n)    = -999
      spray_source_j(n)    = -999
      spray_on_the_task(n) = .FALSE.

      CALL latlon_to_ij(                                 &
           proj,                                         &
           spray_source_lat(n),                          &
           spray_source_lon(n),                          &
           i_val, j_val)

      ! NaN check.
      IF (i_val /= i_val .OR. j_val /= j_val) THEN

         WRITE(0,*) 'SPRAY MAPPING RETURNED NAN:',       &
                    ' source=', n,                       &
                    ' domain=', current_grid%id,         &
                    ' requested lat/lon=',               &
                    spray_source_lat(n),                 &
                    spray_source_lon(n),                 &
                    ' anchor i/j=', i_anchor, j_anchor,  &
                    ' anchor lat/lon=',                  &
                    known_lat, known_lon,                &
                    ' map_proj=', config_flags%map_proj, &
                    ' dx=', current_grid%dx,             &
                    ' truelat1/2=',                      &
                    model_config_rec%truelat1(current_grid%id),               &
                    model_config_rec%truelat2(current_grid%id),               &
                    ' stand_lon=',                       &
                    model_config_rec%stand_lon(current_grid%id)

         CALL wrf_error_fatal(                           &
              'locate_spray_sources: latlon_to_ij returned NaN')

      END IF

      isrc = NINT(i_val)
      jsrc = NINT(j_val)

      spray_source_i(n) = isrc
      spray_source_j(n) = jsrc

      ! Check against the physical unstaggered domain.
      IF (isrc < current_grid%sd31 .OR.                  &
          isrc > current_grid%ed31-1 .OR.                &
          jsrc < current_grid%sd33 .OR.                  &
          jsrc > current_grid%ed33-1) THEN

         WRITE(0,*) 'SPRAY SOURCE OUTSIDE DOMAIN:',      &
                    ' source=', n,                       &
                    ' domain=', current_grid%id,         &
                    ' requested lat/lon=',               &
                    spray_source_lat(n),                 &
                    spray_source_lon(n),                 &
                    ' real i/j=', i_val, j_val,          &
                    ' integer i/j=', isrc, jsrc,         &
                    ' domain bounds=',                   &
                    current_grid%sd31,                   &
                    current_grid%ed31-1,                 &
                    current_grid%sd33,                   &
                    current_grid%ed33-1

         CALL wrf_error_fatal(                           &
              'Spray source lies outside assigned domain')

      END IF

      ! Only the task owning the rounded mass-grid cell emits it.
      on_patch =                                         &
           isrc >= ips_src .AND. isrc <= ipe_src .AND.   &
           jsrc >= jps_src .AND. jsrc <= jpe_src

      spray_on_the_task(n) = on_patch

      WRITE(0,'(A,I3,A,I2)')                            &
           'Source ', n, ' domain ', current_grid%id
      WRITE(0,'(A,2F11.5)')                             &
           '  Requested lat/lon: ', spray_source_lat(n), spray_source_lon(n)
      WRITE(0,'(A,2F10.3)')                             &
           '  Real i/j: ', i_val, j_val
      WRITE(0,'(A,2I6)')                                &
           '  Integer i/j: ', isrc, jsrc
      WRITE(0,'(A,4I6)')                                &
           '  Patch bounds (i,j): ', ips_src, ipe_src, jps_src, jpe_src
      WRITE(0,'(A,L1)')                                 &
           '  On this task? ', on_patch
      WRITE(0,*) ''


   END DO

END SUBROUTINE locate_spray_sources

SUBROUTINE add_spray_increment(                          &
     qscalar, is, grid_id, num_spray_sources, dt_model,  &
     alt, ph, phb, msftx, msfty, dx, dy,                 &
     ids, ide, jds, jde, kds, kde,                       &
     ims, ime, jms, jme, kms, kme,                       &
     its, ite, jts, jte, kts, kte )

   USE module_state_description, ONLY :                  &
        P_QNSPRAER, P_QASPRAER, P_QMSPRAER

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: is
   INTEGER, INTENT(IN) :: grid_id
   INTEGER, INTENT(IN) :: num_spray_sources

   INTEGER, INTENT(IN) :: ids, ide, jds, jde, kds, kde
   INTEGER, INTENT(IN) :: ims, ime, jms, jme, kms, kme
   INTEGER, INTENT(IN) :: its, ite, jts, jte, kts, kte

   REAL, INTENT(IN) :: dt_model
   REAL, INTENT(IN) :: dx, dy

   REAL, DIMENSION(ims:ime,kms:kme,jms:jme),             &
         INTENT(INOUT) :: qscalar

   REAL, DIMENSION(ims:ime,kms:kme,jms:jme),             &
         INTENT(IN) :: alt

   REAL, DIMENSION(ims:ime,kms:kme,jms:jme),             &
         INTENT(IN) :: ph, phb

   REAL, DIMENSION(ims:ime,jms:jme), INTENT(IN) ::        &
         msftx, msfty

   INTEGER :: n
   INTEGER :: i, j, k

   REAL :: ndot
   REAL :: adot
   REAL :: mdot
   REAL :: source_rate

   REAL :: dg
   REAL :: sigma_g
   REAL :: ln_sigma_g
   REAL :: rho_particle

   REAL :: cell_area
   REAL :: dz
   REAL :: rho_dry
   REAL :: dry_air_mass
   REAL :: delta_q

   REAL, PARAMETER :: pi_spray = 3.14159265358979323846
   REAL, PARAMETER :: gravity  = 9.80665

   IF (num_spray_sources <= 0) RETURN

   DO n = 1, num_spray_sources

      ! Source is emitted only on its assigned WRF domain.
      IF (spray_source_domain(n) /= grid_id) CYCLE

      ! Only one MPI task should own this source.
      IF (.NOT. spray_on_the_task(n)) CYCLE

      IF (spray_number_rate(n) <= 0.0) CYCLE

      i = spray_source_i(n)
      j = spray_source_j(n)

      ! One MPI patch may contain multiple OpenMP tiles.
      ! Only the tile containing this point should update it.
      IF (i < its .OR. i > ite) CYCLE
      IF (j < jts .OR. j > jte) CYCLE

      ! Check the physical unstaggered mass-grid domain.
      IF (i < ids .OR. i > ide-1) CYCLE
      IF (j < jds .OR. j > jde-1) CYCLE

      ! Check that the point is available in this task's memory.
      IF (i < ims .OR. i > ime) THEN
         CALL wrf_error_fatal(                           &
              'Spray source i index outside local memory')
      END IF

      IF (j < jms .OR. j > jme) THEN
         CALL wrf_error_fatal(                           &
              'Spray source j index outside local memory')
      END IF

      ! Emit into the lowest unstaggered model level.
      k = MAX(kts,kds)

      IF (k > MIN(kte,kde-1)) CYCLE

      !------------------------------------------------------------
      ! Calculate total number, surface-area and dry-mass rates.
      !------------------------------------------------------------

      ndot         = spray_number_rate(n)
      dg           = spray_mode_diameter(n)
      sigma_g      = spray_geometric_sigma(n)
      rho_particle = spray_particle_density(n)

      IF (dg <= 0.0) THEN
         CALL wrf_error_fatal(                           &
              'Active spray source has non-positive diameter')
      END IF

      IF (sigma_g < 1.0) THEN
         CALL wrf_error_fatal(                           &
              'Active spray source has sigma_g below one')
      END IF

      IF (rho_particle <= 0.0) THEN
         CALL wrf_error_fatal(                           &
              'Active spray source has non-positive density')
      END IF

      ln_sigma_g = LOG(sigma_g)

      ! Total geometric surface area emitted per second:
      ! units m2 s-1.
      adot = ndot * pi_spray * dg**2 *                   &
             EXP(2.0 * ln_sigma_g**2)

      ! Total dry aerosol mass emitted per second:
      ! units kg s-1.
      mdot = ndot * (pi_spray / 6.0) * rho_particle *    &
             dg**3 * EXP(4.5 * ln_sigma_g**2)

      IF (is == P_QNSPRAER) THEN

         ! particles s-1
         source_rate = ndot

      ELSE IF (is == P_QASPRAER) THEN

         ! m2 s-1
         source_rate = adot

      ELSE IF (is == P_QMSPRAER) THEN

         ! kg s-1
         source_rate = mdot

      ELSE

         CYCLE

      END IF

      !------------------------------------------------------------
      ! Dry-air mass in the receiving WRF grid cell.
      !------------------------------------------------------------

      ! Physical horizontal area, accounting for map factors.
      cell_area = dx * dy /                              &
                  (msftx(i,j) * msfty(i,j))

      IF (cell_area <= 0.0) THEN
         CALL wrf_error_fatal(                           &
              'Non-positive spray-source cell area')
      END IF

      ! PH and PHB are geopotential on vertically staggered levels.
      dz = ( ph(i,k+1,j) + phb(i,k+1,j)                  &
           - ph(i,k  ,j) - phb(i,k  ,j) ) / gravity

      IF (dz <= 0.0) THEN
         CALL wrf_error_fatal(                           &
              'Non-positive spray-source layer depth')
      END IF

      ! ALT is inverse dry-air density, m3 kg-1.
      rho_dry = 1.0 / alt(i,k,j)

      IF (rho_dry <= 0.0) THEN
         CALL wrf_error_fatal(                           &
              'Non-positive dry-air density in spray cell')
      END IF

      dry_air_mass = rho_dry * cell_area * dz

      IF (dry_air_mass <= 0.0) THEN
         CALL wrf_error_fatal(                           &
              'Non-positive dry-air mass in spray cell')
      END IF

      !------------------------------------------------------------
      ! Convert the absolute rate into a mixing-ratio increment.
      !
      ! QNSPRAER: particles kg-1
      ! QASPRAER: m2 kg-1
      ! QMSPRAER: kg kg-1
      !------------------------------------------------------------

      delta_q = source_rate * dt_model / dry_air_mass

      qscalar(i,k,j) = MAX(0.0,qscalar(i,k,j)) + delta_q

      ! Temporary conservation diagnostic for the number field.
      IF (is == P_QNSPRAER) THEN

         WRITE(0,*) 'SPRAY APPLY:',                      &
                    ' source=', n,                       &
                    ' domain=', grid_id,                 &
                    ' i/j/k=', i, j, k,                  &
                    ' dt=', dt_model,                    &
                    ' ndot=', ndot,                      &
                    ' air mass=', dry_air_mass,          &
                    ' expected particles=', ndot*dt_model, &
                    ' recovered particles=',             &
                    delta_q*dry_air_mass

      END IF

   END DO

END SUBROUTINE add_spray_increment


   SUBROUTINE initialise_spray_defaults()

      spray_source_domain(:)      = 0
      spray_source_lat(:)         = 0.0
      spray_source_lon(:)         = 0.0
      spray_number_rate(:)        = 0.0
      spray_mode_diameter(:)      = 0.0
      spray_geometric_sigma(:)    = 1.0
      spray_particle_density(:)   = 2160.0
      spray_source_i(:)    = -999
      spray_source_j(:)    = -999
      spray_on_the_task(:)    = .false.

   END SUBROUTINE initialise_spray_defaults

! Model layer, even though it does I/O -- special case of namelist I/O.

   SUBROUTINE initial_config
!<DESCRIPTION>
! This routine reads in the namelist.input file and sets
! module_config_rec, a structure of TYPE(model_config_rec_type), which is is seen via USE association by any
! subprogram that uses module_configure.  The module_config_rec structure
! contains all namelist settings for all domains.  Variables that apply
! to the entire run and have only one value regardless of domain are
! scalars.  Variables that allow different settings for each domain are
! defined as arrays of dimension max_domains (defined in
! frame/module_driver_constants.F, from a setting passed in from
! configure.wrf). There is another type in WRF, TYPE(grid_config_rec_type), in which
! all fields pertain only to a single domain (and are all scalars). The subroutine
! model_to_grid_config_rec(), also in frame/module_configure.F, is used to retrieve
! the settings for a given domain from a TYPE(module_config_rec_type) and put them into
! a TYPE(grid_config_rec_type), variables of which type are often called <em>config_flags</em>
! in the WRF code.
!
! Most of the code in this routine is generated from the Registry file
! rconfig entries and included from the following files (found in the inc directory):
!
! <pre>
! namelist_defines.inc	declarations of namelist variables (local to this routine)
! namelist_statements.inc	NAMELIST statements for each variable
! namelist_defaults.inc	assignment to default values if specified in Registry
! config_reads.inc		read statements for each namelist record
! config_assigns.inc	assign each variable to field in module_config_rec
! </pre>
!
!NOTE: generated subroutines from Registry file rconfig entries are renamed nl_
! instead of rconfig_ due to length limits for subroutine names.
!
! Note for version WRF 2.0: there is code here to force all domains to
! have the same mp_physics setting. This is because different mp_physics
! packages have different numbers of tracers but the nest forcing and
! feedback code relies on the parent and nest having the same number and
! kind of tracers. This means that the microphysics option
! specified on the highest numbered domain is the microphysics
! option for <em>all</em> domains in the run. This will be revisited.
!
!</DESCRIPTION>
      IMPLICIT NONE

      INTEGER              :: io_status
      INTEGER              :: i

      LOGICAL              :: nml_read_error

      CHARACTER (LEN=1024) :: nml_name

      INTEGER, PARAMETER :: nml_write_unit= 9
      INTEGER, PARAMETER :: nml_read_unit = 10

      CHARACTER (LEN=1024) :: entire_line


! define as temporaries
#include "namelist_defines.inc"


! Statements that specify the namelists
#include "namelist_statements.inc"

     NAMELIST /physics/                 spray_source_domain                           
     NAMELIST /physics/                 spray_source_lat                             
     NAMELIST /physics/                 spray_source_lon                           
     NAMELIST /physics/                 spray_number_rate                           
     NAMELIST /physics/                 spray_mode_diameter                        
     NAMELIST /physics/                 spray_geometric_sigma                       
     NAMELIST /physics/                 spray_particle_density

     CALL initialise_spray_defaults()

      OPEN ( UNIT   = nml_read_unit    ,      &
             FILE   = "namelist.input" ,      &
             FORM   = "FORMATTED"      ,      &
             STATUS = "OLD"            ,      &
             IOSTAT = io_status         )

      IF ( io_status .NE. 0 ) THEN
        CALL WRF_ERROR_FATAL ( 'ERROR OPENING namelist.input' )
      ENDIF

#ifndef NO_NAMELIST_PRINT
      OPEN ( UNIT   = nml_write_unit    ,      &
#if (DA_CORE == 1)
             FILE   = "namelist.output.da" ,      &
#else
             FILE   = "namelist.output" ,      &
#endif
             FORM   = "FORMATTED"      ,      &
             STATUS = "REPLACE"        ,      &
             IOSTAT = io_status         )

      IF ( io_status .NE. 0 ) THEN
#if (DA_CORE == 1)
        CALL WRF_ERROR_FATAL ( 'ERROR OPENING namelist.output.da' )
#else
        CALL WRF_ERROR_FATAL ( 'ERROR OPENING namelist.output' )
#endif
      ENDIF
#endif

! Statements that set the namelist vars to default vals
#  include "namelist_defaults.inc"

#if (DA_CORE == 1)
! Override the default values, because we can not assigned a arrary with different values in registry.

      as1(1:3) = (/ 0.25, 1.0, 1.5 /)
      as2(1:3) = (/ 0.25, 1.0, 1.5 /)
      as3(1:3) = (/ 0.25, 1.0, 1.5 /)
      as4(1:3) = (/ 0.25, 1.0, 1.5 /)
      as5(1:3) = (/ 0.25, 1.0, 1.5 /)
#endif

! Statements that read the namelist are in this file
#  include "config_reads.inc"

! 2004/04/28  JM (with consensus by the group of developers)
! This is needed to ensure that nesting will work, since
! different mp_physics packages have different numbers of
! tracers. Basically, this says that the microphysics option
! specified on the highest numbered domain *is* the microphysics
! option for the run. Not the best solution but okay for 2.0.
!

      DO i = 1, max_dom
         mp_physics(i) = mp_physics(max_dom)
      ENDDO

! Statements that assign the variables to the cfg record are in this file
! except the namelist_derived variables where are assigned below
#undef SOURCE_RECORD
#undef DEST_RECORD
#undef SOURCE_REC_DEX
#define SOURCE_RECORD
#define DEST_RECORD model_config_rec %
#define SOURCE_REC_DEX
#include "config_assigns.inc"


      CLOSE ( UNIT = nml_read_unit , IOSTAT = io_status )

      IF ( io_status .NE. 0 ) THEN
        CALL WRF_ERROR_FATAL ( 'ERROR CLOSING namelist.input' )
      ENDIF
        num_spray_sources1 = model_config_rec%num_spray_sources

#ifndef NO_NAMELIST_PRINT
      CLOSE ( UNIT = nml_write_unit , IOSTAT = io_status )

      IF ( io_status .NE. 0 ) THEN
        CALL WRF_ERROR_FATAL ( 'ERROR CLOSING namelist.output' )
      ENDIF
#endif

#ifdef _WIN32
      model_config_rec%nocolons = .TRUE.   ! always no colons for Windows
#endif

#ifdef WRF_CFBM
      do i = 1, max_dom
        if (model_config_rec%ifire(i) == 1) then
          call cfbm_namelist % Initialization( "namelist.cfbm" )
          exit
        end if
      end do
#endif

      RETURN

   END SUBROUTINE initial_config

#if 1
   SUBROUTINE get_config_as_buffer( buffer, buflen, ncopied )
! note that model_config_rec_type must be defined as a sequence derived type
      INTEGER,   INTENT(INOUT) ::  buffer(*)
      INTEGER,   INTENT(IN)    ::  buflen
      INTEGER,   INTENT(OUT)   ::  ncopied
!      TYPE(model_config_rec_type) :: model_config_rec
      INTEGER :: nbytes
      CALL wrf_num_bytes_between ( model_config_rec%last_item_in_struct ,   &
                                   model_config_rec%first_item_in_struct ,  &
                                   nbytes )
!      nbytes = loc(model_config_rec%last_item_in_struct) - &
!               loc(model_config_rec%first_item_in_struct)
      IF ( nbytes .gt. buflen ) THEN
        CALL wrf_error_fatal( &
        "get_config_rec_as_buffer: buffer size too small for config_rec" )
      ENDIF
      CALL wrf_mem_copy( model_config_rec, buffer, nbytes )
      ncopied = nbytes
      RETURN
   END SUBROUTINE get_config_as_buffer

   SUBROUTINE set_config_as_buffer( buffer, buflen )
! note that model_config_rec_type must be defined as a sequence derived type
      INTEGER,   INTENT(INOUT) ::  buffer(*)
      INTEGER,   INTENT(IN)    ::  buflen
!      TYPE(model_config_rec_type) :: model_config_rec
      INTEGER :: nbytes
      CALL wrf_num_bytes_between ( model_config_rec%last_item_in_struct ,  &
                                   model_config_rec%first_item_in_struct , &
                                   nbytes )
!      nbytes = loc(model_config_rec%last_item_in_struct) - &
!               loc(model_config_rec%first_item_in_struct)
      IF ( nbytes .gt. buflen ) THEN
        CALL wrf_error_fatal( &
        "set_config_rec_as_buffer: buffer length too small to fill model config record" )
      ENDIF
      CALL wrf_mem_copy( buffer, model_config_rec, nbytes )
      RETURN
   END SUBROUTINE set_config_as_buffer
#else
   SUBROUTINE get_config_as_buffer( buffer, buflen, ncopied )
! note that model_config_rec_type must be defined as a sequence derived type
      INTEGER*1, INTENT(INOUT) ::  buffer(*)
      INTEGER,   INTENT(IN)    ::  buflen
      INTEGER,   INTENT(OUT)   ::  ncopied
!      TYPE(model_config_rec_type) :: model_config_rec
      INTEGER :: nbytes
      nbytes = loc(model_config_rec%last_item_in_struct) - &
               loc(model_config_rec%first_item_in_struct)
      IF ( nbytes .gt. buflen ) THEN
        CALL wrf_error_fatal( &
        "get_config_rec_as_buffer: buffer size too small for config_rec" )
      ENDIF
      CALL wrf_mem_copy( model_config_rec, buffer, nbytes )
      ncopied = nbytes
      RETURN
   END SUBROUTINE get_config_as_buffer

   SUBROUTINE set_config_as_buffer( buffer, buflen )
! note that model_config_rec_type must be defined as a sequence derived type
      INTEGER*1, INTENT(INOUT) ::  buffer(*)
      INTEGER,   INTENT(IN)    ::  buflen
!      TYPE(model_config_rec_type) :: model_config_rec
      INTEGER :: nbytes
      nbytes = loc(model_config_rec%last_item_in_struct) - &
               loc(model_config_rec%first_item_in_struct)
      IF ( nbytes .gt. buflen ) THEN
        CALL wrf_error_fatal( &
        "set_config_rec_as_buffer: buffer length too small to fill model config record" )
      ENDIF
      CALL wrf_mem_copy( buffer, model_config_rec, nbytes )
      RETURN
   END SUBROUTINE set_config_as_buffer
#endif

   SUBROUTINE model_to_grid_config_rec ( id_id , model_config_rec , grid_config_rec )
      INTEGER , INTENT(IN)                         ::  id_id
      TYPE ( model_config_rec_type ) , INTENT(IN)  ::  model_config_rec
      TYPE ( grid_config_rec_type  ) , INTENT(OUT) ::  grid_config_rec
! <DESCRIPTION>
! This routine is called to populate a domain specific configuration
! record of TYPE(grid_config_rec_type) with the configuration information
! for that domain that is stored in TYPE(model_config_rec). Both types
! are defined in frame/module_configure.F.  The input argument is the
! record of type model_config_rec_type contains the model-wide
! configuration information (that is, settings that apply to the model in
! general) and configuration information for each individual domain.  The
! output argument is the record of type grid_config_rec_type which
! contains the model-wide configuration information and the
! domain-specific information for this domain only.  In the
! model_config_rec, the domain specific information is arrays, indexed by
! the grid id's.  In the grid_config_rec the domain-specific information
! is scalar and for the specific domain.  The first argument to this
! routine is the grid id (top-most domain is always 1) as specified in
! the domain-specific namelist variable grid_id.
!
! The actual assignments form the model_config_rec_type to the
! grid_config_rec_type are generate from the rconfig entries in the
! Registry file and included by this routine from the file
! inc/config_assigns.inc.
!
!NOTE: generated subroutines from Registry file rconfig entries are renamed nl_
! instead of rconfig_ due to length limits for subroutine names.
!
!
! </DESCRIPTION>
#undef SOURCE_RECORD
#undef SOURCE_REC_DEX
#undef DEST_RECORD
#define SOURCE_RECORD model_config_rec %
#define SOURCE_REC_DEX (id_id)
#define DEST_RECORD   grid_config_rec %
#include "config_assigns.inc"
   END SUBROUTINE model_to_grid_config_rec


   FUNCTION in_use_for_config ( id, vname ) RESULT ( in_use )
     INTEGER, INTENT(IN) :: id
     CHARACTER*(*), INTENT(IN) :: vname
     LOGICAL in_use
     INTEGER uses

     uses = 0
     in_use = .TRUE.

     IF      ( vname(1:1) .EQ. 'a' ) THEN
#    include "in_use_for_config_a.inc"
     ELSE IF ( vname(1:1) .EQ. 'b' ) THEN
#    include "in_use_for_config_b.inc"
     ELSE IF ( vname(1:1) .EQ. 'c' ) THEN
#    include "in_use_for_config_c.inc"
     ELSE IF ( vname(1:1) .EQ. 'd' ) THEN
#    include "in_use_for_config_d.inc"
     ELSE IF ( vname(1:1) .EQ. 'e' ) THEN
#    include "in_use_for_config_e.inc"
     ELSE IF ( vname(1:1) .EQ. 'f' ) THEN
#    include "in_use_for_config_f.inc"
     ELSE IF ( vname(1:1) .EQ. 'g' ) THEN
#    include "in_use_for_config_g.inc"
     ELSE IF ( vname(1:1) .EQ. 'h' ) THEN
#    include "in_use_for_config_h.inc"
     ELSE IF ( vname(1:1) .EQ. 'i' ) THEN
#    include "in_use_for_config_i.inc"
     ELSE IF ( vname(1:1) .EQ. 'j' ) THEN
#    include "in_use_for_config_j.inc"
     ELSE IF ( vname(1:1) .EQ. 'k' ) THEN
#    include "in_use_for_config_k.inc"
     ELSE IF ( vname(1:1) .EQ. 'l' ) THEN
#    include "in_use_for_config_l.inc"
     ELSE IF ( vname(1:1) .EQ. 'm' ) THEN
#    include "in_use_for_config_m.inc"
     ELSE IF ( vname(1:1) .EQ. 'n' ) THEN
#    include "in_use_for_config_n.inc"
     ELSE IF ( vname(1:1) .EQ. 'o' ) THEN
#    include "in_use_for_config_o.inc"
     ELSE IF ( vname(1:1) .EQ. 'p' ) THEN
#    include "in_use_for_config_p.inc"
     ELSE IF ( vname(1:1) .EQ. 'q' ) THEN
#    include "in_use_for_config_q.inc"
     ELSE IF ( vname(1:1) .EQ. 'r' ) THEN
#    include "in_use_for_config_r.inc"
     ELSE IF ( vname(1:1) .EQ. 's' ) THEN
#    include "in_use_for_config_s.inc"
     ELSE IF ( vname(1:1) .EQ. 't' ) THEN
#    include "in_use_for_config_t.inc"
     ELSE IF ( vname(1:1) .EQ. 'u' ) THEN
#    include "in_use_for_config_u.inc"
     ELSE IF ( vname(1:1) .EQ. 'v' ) THEN
#    include "in_use_for_config_v.inc"
     ELSE IF ( vname(1:1) .EQ. 'w' ) THEN
#    include "in_use_for_config_w.inc"
     ELSE IF ( vname(1:1) .EQ. 'x' ) THEN
#    include "in_use_for_config_x.inc"
     ELSE IF ( vname(1:1) .EQ. 'y' ) THEN
#    include "in_use_for_config_y.inc"
     ELSE IF ( vname(1:1) .EQ. 'z' ) THEN
#    include "in_use_for_config_z.inc"
     ENDIF

     RETURN
   END FUNCTION


! Include the definitions of all the routines that return a namelist values
! back to the driver. These are generated by the registry

   SUBROUTINE init_module_configure
     USE module_scalar_tables
     IMPLICIT NONE
     CALL init_module_scalar_tables
   END SUBROUTINE init_module_configure

   SUBROUTINE wrf_alt_nml_obsolete (nml_read_unit, nml_name)
!
!<DESCRIPTION>
! If there is an error reading the "nml_name" namelist, this routine is
! called to check for namelist variables that have been removed by the 
! developers and are still in user's namelists.
! 
! The calls to this routine are in registry-generated code: inc/config_reads.inc
!</DESCRIPTION>
!
     IMPLICIT NONE
     INTEGER, INTENT(IN)       :: nml_read_unit
     CHARACTER*(*), INTENT(IN) :: nml_name
     INTEGER                   :: nml_error

#include "namelist_defines.inc"
#include "namelist_statements.inc"

! These are the variables that have been removed
     logical , DIMENSION(max_domains) :: pd_moist, pd_chem, pd_tke, pd_scalar
     NAMELIST /dynamics/                 pd_moist, pd_chem, pd_tke, pd_scalar

     integer , DIMENSION(max_domains) :: ucmcall
     NAMELIST /physics/                  ucmcall

     integer , DIMENSION(max_domains) :: obs_nobs_prt
     NAMELIST /fdda/                     obs_nobs_prt

     LOGICAL ::         global, print_detail_airep, print_detail_timing
     NAMELIST /wrfvar1/ global, print_detail_airep, print_detail_timing

     LOGICAL ::         write_qcw, write_qrn, write_qci, write_qsn
     NAMELIST /wrfvar2/ write_qcw, write_qrn, write_qci, write_qsn
     LOGICAL ::          write_qgr, write_filtered_obs
     NAMELIST /wrfvar2/  write_qgr, write_filtered_obs

     LOGICAL ::         use_eos_radobs, use_3dvar_phy
     NAMELIST /wrfvar4/ use_eos_radobs, use_3dvar_phy

     LOGICAL             :: use_crtm_kmatrix_fast
     NAMELIST /wrfvar14/    use_crtm_kmatrix_fast
     CHARACTER (LEN=256) :: spccoeff_file, taucoeff_file, aerosolcoeff_file
     NAMELIST /wrfvar14/    spccoeff_file, taucoeff_file, aerosolcoeff_file
     CHARACTER (LEN=256) :: cloudcoeff_file, emiscoeff_file
     NAMELIST /wrfvar14/    cloudcoeff_file, emiscoeff_file

     LOGICAL             :: alpha_vertloc
     NAMELIST /wrfvar16/    alpha_vertloc


! Read the namelist again, if it succeeds after adding the above variables,
! it probably failed because these are still in the namelist.  If it fails
! again, we will return.

     REWIND ( UNIT = nml_read_unit )

!----------------------------- dynamics ---------------------------------
     if ( TRIM(nml_name) .eq. "dynamics" ) then

        READ   ( UNIT = nml_read_unit , NML = dynamics , iostat=nml_error )

        IF ( nml_error .EQ. 0 ) then    ! Successul, rm variables must be problem
           CALL wrf_debug(0, "-- Are pd_moist, pd_chem, pd_tke, or pd_scalar still in your "// &
                               TRIM(nml_name)//" namelist?")
           CALL wrf_debug(0, "-- Replace them with moist_adv_opt, chem_adv_opt, tke_adv_opt "// &
                             " and scalar_adv_opt, respectively.")
        ENDIF

!---------------------------------- physics -----------------------------
     else if ( TRIM(nml_name) .eq. "physics" ) then

        READ   ( UNIT = nml_read_unit , NML = physics , iostat=nml_error )

        IF ( nml_error .EQ. 0 ) then    ! Successul, rm variables must be problem
           CALL wrf_debug(0, "-- Is ucmcall still in your "// TRIM(nml_name)//" namelist?")
           CALL wrf_debug(0, "-- Replace it with sf_urban_physics")
        ENDIF

!---------------------------------- fdda --------------------------------
     else if ( TRIM(nml_name) .eq. "fdda" ) then

        READ   ( UNIT = nml_read_unit , NML = fdda , iostat=nml_error )

        IF ( nml_error .EQ. 0 ) then    ! Successul, rm variables must be problem
           CALL wrf_debug(0, "-- Is obs_nobs_prt still in your "// TRIM(nml_name)//" namelist?")
           CALL wrf_debug(0, "-- Replace it with obs_prt_max")
        ENDIF

!---------------------------------- wrfvar1 -----------------------------
     else if ( TRIM(nml_name) .eq. "wrfvar1" ) then

        READ   ( UNIT = nml_read_unit , NML = wrfvar1 , iostat=nml_error )

        IF ( nml_error .EQ. 0 ) then    ! Successul, rm variables must be problem
           CALL wrf_debug(0, "-- Are global, print_detail_airep, print_detail_timing still in your "// &
                              TRIM(nml_name)//" namelist?")
           CALL wrf_debug(0, "-- Remove global, print_detail_airep, print_detail_timing "// &
                             "from wrfvar1 namelist as they are obsolete.")
        ENDIF

!---------------------------------- wrfvar2 -----------------------------
     else if ( TRIM(nml_name) .eq. "wrfvar2" ) then

        READ   ( UNIT = nml_read_unit , NML = wrfvar2 , iostat=nml_error )

        IF ( nml_error .EQ. 0 ) then    ! Successul, rm variables must be problem
           CALL wrf_debug(0, "-- Are write_qcw, write_qrn, write_qci, write_qsn, write_qgr, "// &
                             "write_filtered_obs still in your "// &
                              TRIM(nml_name)//" namelist?")
           CALL wrf_debug(0, "-- Remove write_qcw, write_qrn, write_qci, write_qsn, write_qgr, "// &
                             "write_filtered_obs as they are obsolete.")
        ENDIF

!---------------------------------- wrfvar4 -----------------------------
     else if ( TRIM(nml_name) .eq. "wrfvar4" ) then

        READ   ( UNIT = nml_read_unit , NML = wrfvar4 , iostat=nml_error )

        IF ( nml_error .EQ. 0 ) then    ! Successul, rm variables must be problem
           CALL wrf_debug(0, "-- Is use_3dvar_phy or use_eos_radobs still in your "// &
                              TRIM(nml_name)//" namelist?")
           CALL wrf_debug(0, "-- Remove use_3dvar_phy, use_eos_radobs as they are obsolete.")
        ENDIF

!---------------------------------- wrfvar14 -----------------------------
     else if ( TRIM(nml_name) .eq. "wrfvar14" ) then

     READ   ( UNIT = nml_read_unit , NML = wrfvar14 , iostat=nml_error )

        IF ( nml_error .EQ. 0 ) then    ! Successul, rm variables must be problem
           CALL wrf_debug(0, "-- Are use_crtm_kmatrix_fast, spccoeff_file, taucoeff_file, "// &
                             "aerosolcoeff_file, cloudcoeff_file, emiscoeff_file still in your "// &
                              TRIM(nml_name)//" namelist?")
           CALL wrf_debug(0, "-- Remove them as they are obsolete.")
        ENDIF

!---------------------------------- wrfvar16 -----------------------------
     else if ( TRIM(nml_name) .eq. "wrfvar16" ) then

     READ   ( UNIT = nml_read_unit , NML = wrfvar16 , iostat=nml_error )

        IF ( nml_error .EQ. 0 ) then    ! Successul, rm variables must be problem
           CALL wrf_debug(0, "-- Is alpha_vertloc still in your "// &
                              TRIM(nml_name)//" namelist?")
           CALL wrf_debug(0, "-- Remove it as it is obsolete.")
           CALL wrf_debug(0, "-- For EnVar DA (ensdim_alpha>0), set proper alpha_vertloc_opt instead.")
        ENDIF

!---------------------------------- error -------------------------------
     else
         IF ( &
#include "namelist_nametest.inc"
              ) THEN
            nml_error = 0
         ELSE
            CALL wrf_debug(0, TRIM(nml_name)//" is not a valid namelist name")
         ENDIF
     end if

     IF ( nml_error .NE. 0 ) then    ! Still failed
        return
     ENDIF

   END SUBROUTINE wrf_alt_nml_obsolete

END MODULE module_configure

#if( WRF_CHEM == 1 && WRF_KPP == 1 )
INTEGER FUNCTION get_funit()

    IMPLICIT NONE

    INTEGER, PARAMETER :: min_unit_number = 30
    INTEGER, PARAMETER :: max_unit_number = 99 

    LOGICAL :: opened

    DO get_funit = min_unit_number, max_unit_number
       INQUIRE(UNIT=get_funit, OPENED=opened)
       IF ( .NOT. opened ) RETURN
    END DO

    get_funit = -1

END FUNCTION get_funit
#endif

SUBROUTINE set_scalar_indices_from_config ( idomain , dummy2, dummy1 )
  USE module_driver_constants
  USE module_state_description
  USE module_wrf_error
  USE module_configure, ONLY : model_config_rec
  USE module_scalar_tables
#if( WRF_CHEM == 1 && WRF_KPP == 1 )
  USE module_irr_diag
#endif
  IMPLICIT NONE
  INTEGER , INTENT(IN)  :: idomain
  INTEGER               :: dummy1
  INTEGER               :: dummy2

#if( WRF_CHEM == 1 && WRF_KPP == 1 )
  INTEGER :: astat, istat
  INTEGER :: l, m, n
  INTEGER :: ndx
  INTEGER :: irr_cnt
  INTEGER :: nFound
  INTEGER :: funit
  INTEGER :: chem_opt
  INTEGER, ALLOCATABLE :: irr_ndx(:)
  CHARACTER(LEN=128) :: err_mes
  CHARACTER(LEN=128) :: fname
  CHARACTER(LEN=64)  :: irr_rxt, lc_tmp
  CHARACTER(LEN=64), POINTER     :: rxt_lst(:)
  CHARACTER(LEN=32), ALLOCATABLE :: irr_lst(:)
  LOGICAL, ALLOCATABLE :: irr_msk(:)
  LOGICAL :: found, opened
  LOGICAL :: do_all

  INTEGER :: get_funit

#include "scalar_indices_irr_diag_decls.inc"
#endif

!<DESCRIPTION>
!This routine is called to adjust the integer variables that are defined
!in frame/module_state_description.F (Registry-generated) and that serve
!as indices into 4D tracer arrays for moisture, chemistry, etc.
!Different domains (different grid data structures) are allowed to have
!different sets of tracers so these indices can vary from domain to
!domain. However, since the indices are defined globally in
!module_state_description (a shortcoming in the current software), it is
!necessary that these indices be reset each time a different grid is to
!be computed on.
!
!The scalar idices are set according to the particular physics
!packages -- more specifically in the case of the moisture tracers, microphysics
!packages -- that are stored for each domain in model_config_rec and
!indexed by the grid id, passed in as an argument to this routine.  (The
!initial_config() routine in module_configure is what reads the
!namelist.input file and sets model_config_rec.)
!
!The actual code for calculating the scalar indices on a particular
!domain is generated from the Registry state array definitions for the
!4d tracers and from the package definitions that indicate which physics
!packages use which tracers.
!
!</DESCRIPTION>

#include "scalar_indices.inc"
#if( WRF_CHEM == 1 && WRF_KPP == 1 )
#include "scalar_indices_irr_diag.inc"

   CALL init_module_irr_diag

has_irr_diag: &
   IF( model_config_rec%irr_opt(idomain) == 1 .and. .not. irr_active(idomain) ) THEN
     do_all = .false.
     irr_active(idomain) = .true.
     CALL nl_get_chem_opt( idomain, chem_opt )
     DO chm_opt_ndx = 1,nchm_opts
       IF( chem_opt == chm_opts_ndx(chm_opt_ndx) ) THEN
         EXIT
       ENDIF
     ENDDO
     IF( chm_opt_ndx > nchm_opts ) THEN
       write(err_mes,'(''IRR not supported for chem option '',i3.3)') chem_opt
       CALL wrf_error_fatal( trim(err_mes) )
     ENDIF
     funit = get_funit()
     write(fname,'(''wrfchem_irr_diag_d'',i2.2)') idomain
     OPEN( unit = funit, file=trim(fname), status='OLD', iostat=istat )
     IF( istat /= 0 ) THEN
       write(err_mes,'(''IRR failed to open '',a)') trim(fname)
       CALL wrf_error_fatal( trim(err_mes) )
     ENDIF
     irr_cnt = 0
     DO
       read(funit,'(a)',iostat=istat) irr_rxt
       IF( istat /= 0 ) THEN
         EXIT
       ENDIF
       irr_cnt = irr_cnt + 1
       if( irr_cnt == 1 ) then
         CALL change_to_lower_case( irr_rxt, lc_tmp )
         if( trim(lc_tmp) == 'all' ) then
           do_all  = .true.
           irr_cnt = chm_opts_cnt(chm_opt_ndx) 
           write(err_mes,'(''IRR do_all, irr_cnt = '',l,i4)') do_all, irr_cnt
           CALL wrf_debug( 0,trim(err_mes) )
           EXIT
         endif
       endif
     ENDDO

has_irr_input: &
     IF( irr_cnt == 0 ) THEN
       CALL nl_set_irr_opt( idomain, 0 )
     ELSE has_irr_input
       ALLOCATE( irr_lst(irr_cnt),irr_msk(irr_cnt),irr_ndx(irr_cnt),stat=astat )
       IF( astat /= 0 ) THEN
         write(err_mes,'(''IRR failed to allocate irr_lst .. irr_ndx; error = '',i8)') astat
         CALL wrf_error_fatal( trim(err_mes) )
       ENDIF
       if( .not. do_all ) THEN
         REWIND( funit )
       ENDIF
       nFound = 0
       irr_msk(:) = .false.
       irr_ndx(:) = 0
       rxt_lst => rxtsym(:,chm_opt_ndx)
input_irr_loop: &
       DO n = 1,irr_cnt
         IF( .not. do_all ) THEN
           read(funit,'(a)',iostat=istat) irr_lst(n)
           IF( istat /= 0 ) THEN
             write(err_mes,'(''IRR failed to read '',a)') trim(fname)
             CALL wrf_error_fatal( trim(err_mes) )
           ENDIF
           found = .false.
         ELSE
           irr_lst(n) = rxt_lst(n)
           found = .false.
         ENDIF
rxt_loop: &
         DO l = 1,chm_opts_cnt(chm_opt_ndx)
have_match: &
           IF( trim(rxt_lst(l)) == trim(irr_lst(n)) ) THEN
             found = .true.
             nFound = nFound + 1
             irr_ndx(nFound) = l
             irr_msk(nFound) = .true.
             irr_rxt = rxt_lst(l)
             do m = 1,len_trim( irr_rxt )
               if( irr_rxt(m:m) == '+' ) then
                 irr_rxt(m:m) = '_'
               endif
             enddo
             SELECT CASE( chm_opts_ndx(chm_opt_ndx) )
               CASE( MOZCART_KPP )
                 irr_diag_mozcart_num_table(idomain) = irr_diag_mozcart_num_table(idomain) + 1
                 ndx = irr_diag_mozcart_num_table(idomain)
                 irr_diag_mozcart_boundary_table(idomain,ndx) = .FALSE.
                 irr_diag_mozcart_dname_table(idomain,ndx) = trim( irr_rxt ) // '_IRR'
                 irr_diag_mozcart_desc_table(idomain,ndx)  = trim(rxt_lst(l) ) // ' Integrated Reaction Rate'
                 irr_diag_mozcart_units_table(idomain,ndx) = 'ppmv'
                 irr_diag_mozcart_streams_table(idomain,ndx)%stream(1) = 512
                 irr_diag_mozcart_streams_table(idomain,ndx)%stream(2) = 2097152
               CASE( T1_MOZCART_KPP )
                 irr_diag_t1_mozcart_num_table(idomain) = irr_diag_t1_mozcart_num_table(idomain) + 1
                 ndx = irr_diag_t1_mozcart_num_table(idomain)
                 irr_diag_t1_mozcart_boundary_table(idomain,ndx) = .FALSE.
                 irr_diag_t1_mozcart_dname_table(idomain,ndx) = trim( irr_rxt ) // '_IRR'
                 irr_diag_t1_mozcart_desc_table(idomain,ndx)  = trim(rxt_lst(l) ) // ' Integrated Reaction Rate'
                 irr_diag_t1_mozcart_units_table(idomain,ndx) = 'ppmv'
                 irr_diag_t1_mozcart_streams_table(idomain,ndx)%stream(1) = 512
                 irr_diag_t1_mozcart_streams_table(idomain,ndx)%stream(2) = 2097152
               CASE( MOZART_MOSAIC_4BIN_KPP )
                 irr_diag_mozart_mosaic_4bin_num_table(idomain) = irr_diag_mozart_mosaic_4bin_num_table(idomain) + 1
                 ndx = irr_diag_mozart_mosaic_4bin_num_table(idomain)
                 irr_diag_mozart_mosaic_4bin_boundary_table(idomain,ndx) = .FALSE.
                 irr_diag_mozart_mosaic_4bin_dname_table(idomain,ndx) = trim( irr_rxt ) // '_IRR'
                 irr_diag_mozart_mosaic_4bin_desc_table(idomain,ndx)  = trim(rxt_lst(l) ) // ' Integrated Reaction Rate'
                 irr_diag_mozart_mosaic_4bin_units_table(idomain,ndx) = 'ppmv'
                 irr_diag_mozart_mosaic_4bin_streams_table(idomain,ndx)%stream(1) = 512
                 irr_diag_mozart_mosaic_4bin_streams_table(idomain,ndx)%stream(2) = 2097152
               CASE( MOZART_MOSAIC_4BIN_AQ_KPP )
                 irr_diag_mozart_mosaic_4bin_aq_num_table(idomain) = irr_diag_mozart_mosaic_4bin_aq_num_table(idomain) + 1
                 ndx = irr_diag_mozart_mosaic_4bin_aq_num_table(idomain)
                 irr_diag_mozart_mosaic_4bin_aq_boundary_table(idomain,ndx) = .FALSE.
                 irr_diag_mozart_mosaic_4bin_aq_dname_table(idomain,ndx) = trim( irr_rxt ) // '_IRR'
                 irr_diag_mozart_mosaic_4bin_aq_desc_table(idomain,ndx)  = trim( rxt_lst(l) ) // ' Integrated Reaction Rate'
                 irr_diag_mozart_mosaic_4bin_aq_units_table(idomain,ndx) = 'ppmv'
                 irr_diag_mozart_mosaic_4bin_aq_streams_table(idomain,ndx)%stream(1) = 512
                 irr_diag_mozart_mosaic_4bin_aq_streams_table(idomain,ndx)%stream(2) = 2097152
             END SELECT
             EXIT rxt_loop
           ENDIF have_match
         ENDDO rxt_loop
       ENDDO input_irr_loop
       IF( nFound > 0 ) THEN
         irr_option(idomain) = .true.
         irr_diag_ndx(:nFound,idomain) = irr_ndx(:nFound)
         irr_diag_cnt(idomain) = nFound
       ENDIF
     ENDIF has_irr_input

     INQUIRE( unit=funit,opened=opened )
     IF( opened ) THEN
       CLOSE( funit )
     ENDIF

     IF( ALLOCATED( irr_lst ) ) THEN
       DEALLOCATE( irr_lst )
     ENDIF
     IF( ALLOCATED( irr_msk ) ) THEN
       DEALLOCATE( irr_msk )
     ENDIF
     IF( ALLOCATED( irr_ndx ) ) THEN
       DEALLOCATE( irr_ndx )
     ENDIF
   ENDIF has_irr_diag


#endif
#include "scalar_indices_init.inc"
  RETURN
END SUBROUTINE set_scalar_indices_from_config

#if ( DA_CORE != 1 )
   SUBROUTINE change_to_lower_case(instr,outstr)
     CHARACTER*(*) ,INTENT(IN)  :: instr
     CHARACTER*(*) ,INTENT(OUT) :: outstr
!Local
     CHARACTER*1                :: c
     INTEGER       ,PARAMETER   :: upper_to_lower =IACHAR('a')-IACHAR('A')
     INTEGER                    :: i,n,n1
!
     outstr = ' '
     N  = len(instr)
     N1 = len(outstr)
     N  = MIN(N,N1)
     outstr(1:N) = instr(1:N)
     DO i=1,N
       c = instr(i:i)
       if('A'<=c .and. c <='Z') outstr(i:i)=achar(iachar(c)+upper_to_lower)
     ENDDO

   END SUBROUTINE change_to_lower_case
#endif
