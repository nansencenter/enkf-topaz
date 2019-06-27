










module mod_measurement

  integer, parameter, public :: OBSTYPESTRLEN = 5

  type measurement
     real d                       ! Measurement value
     real var                     ! Error variance of measurement
     character(len=OBSTYPESTRLEN) id ! Type, can be one of those:
                                  ! 'SST' 'SLA' 'ICEC' 'SAL' 'TEM'
                                  ! 'GSAL' 'GTEM' 'TSLA'
     real lon                     ! Longitude position
     real lat                     ! Latitude position
     real depth                   ! depths of position 
     integer ipiv                 ! i-pivot point in grid
     integer jpiv                 ! j-pivot point in grid
     integer ns                   ! representativity in mod cells (meas. support)
                                  ! ns=0 means: point measurements
                                  ! used in m_Generate_element_Sij.F90
     real a1                      ! bilinear coefficient (for ni=0)
     real a2                      ! bilinear coefficient
     real a3                      ! bilinear coefficient
     real a4                      ! bilinear coefficient
     logical status               ! active or not
     integer i_orig_grid          ! KAL - orig grid index for ice drift
                                  ! processing
     integer j_orig_grid          ! orig grid index
     real h                       ! PS - layer thickness, sorry for that
     integer date                 ! FanF - age of the data 
     integer orig_id              ! PS - used in superobing
  end type measurement

  type measurement_info
     character(len=OBSTYPESTRLEN) id ! Type, can be one of those:
                                  ! 'SST' 'SLA' 'ICEC' 'SAL' 'TEM'
     real jday,jtime              ! julian date and time in the day
     real lon                     ! Longitude position
     real lat                     ! Latitude position
     real mindep                  ! shallowest depth, -999 is default 
     real maxdep                  ! deepest depth, -999 is default 

     integer ipiv                 ! i-pivot point in grid
     integer jpiv                 ! j-pivot point in grid
     
     character*16  platnum        ! platform number (like profile ID) 
     character*100 inifile        ! file source

     integer signal               ! 0 out of domain,or other reasons:  
                                  !-1: used for location/time QC
                                  !-2: pre/sal/tem QC 
                                  !-3: Sal and density check 
                                  ! 1 ~ 3 quality order;
  end type measurement_info


end module mod_measurement
