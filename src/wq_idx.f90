module wq_idx
    !! Named indices for water-quality state variables stored in `c(i, :, j)`.
    use m_constants, only: nv
    implicit none

    integer, parameter :: IDX_CONDUCTIVITY          = 1
    integer, parameter :: IDX_INORG_SUSP_SOLIDS     = 2
    integer, parameter :: IDX_DISSOLVED_OXYGEN      = 3
    integer, parameter :: IDX_CBOD_SLOW             = 4
    integer, parameter :: IDX_CBOD_FAST             = 5
    integer, parameter :: IDX_ORGANIC_N             = 6
    integer, parameter :: IDX_AMMONIA_N             = 7
    integer, parameter :: IDX_NOX_N                 = 8
    integer, parameter :: IDX_ORGANIC_P             = 9
    integer, parameter :: IDX_SOLUBLE_REACTIVE_P    = 10
    integer, parameter :: IDX_PHYTOPLANKTON         = 11
    integer, parameter :: IDX_PARTICULATE_ORG_MAT   = 12
    integer, parameter :: IDX_PATHOGEN              = 13
    integer, parameter :: IDX_GENERIC_CONST         = 14
    integer, parameter :: IDX_ALKALINITY            = nv - 2
    integer, parameter :: IDX_TOTAL_INORGANIC_C     = nv - 1
    integer, parameter :: IDX_BENTHIC_ALGAE         = nv
    integer, parameter :: IDX_HYPORHEIC_BIOFILM     = IDX_BENTHIC_ALGAE

end module wq_idx
