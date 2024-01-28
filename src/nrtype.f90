
!Data type definition
!public constants

MODULE nrtype
    !Symbolic names for kind types of 4-, 2-, and 1-byte integers:
    INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
    INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
    INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
    !Symbolic names for kind types of single- and double-precision reals:
    INTEGER, PARAMETER :: SP = KIND(1.0)
    INTEGER, PARAMETER :: DP = KIND(1.0D0)
! nrtype.f90
! definition of all data types and constants

    !Symbolic names for kind types of single- and double-precision complex:
    INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
    INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))
    !Symbolic name for kind type of default logical:
    INTEGER, PARAMETER :: LGT = KIND(.true.)
    !Frequently used mathematical constants (with precision to spare):
    !REAL(SP), PARAMETER :: PI=3.141592653589793238462643383279502884197_sp
    !REAL(SP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_sp
    !REAL(SP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_sp
    !REAL(SP), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_sp
    !REAL(SP), PARAMETER :: EULER=0.5772156649015328606065120900824024310422_sp

    ! constants used in Q2K model
    !REAL(DP), PARAMETER :: PII=3.141592653589793238462643383279502884197_dp
    REAL(DP), PARAMETER :: PII=3.14159265358979_DP
    REAL(DP), PARAMETER :: PIO2_D=1.57079632679489661923132169163975144209858_dp
    REAL(DP), PARAMETER :: TWOPI_D=6.283185307179586476925286766559005768394_dp

    REAL(DP), PARAMETER :: CON = 1.74532925199433E-02_dp
    REAL(DP), PARAMETER :: RHOW = 1.0_DP, CPW = 1.0_DP
    REAL(DP), PARAMETER :: ACOEFF = 0.6_DP, RL = 0.03_DP, BOWEN = 0.47_DP
    REAL(DP), PARAMETER :: EPS = 0.97_DP, SIGMA = 0.000000117_DP
    REAL(DP), PARAMETER :: ALPHAS = 0.0035_DP * 86400.0_DP / 10000.0_DP, HSED = 0.1_DP
    REAL(DP), PARAMETER :: ADAM = 1.25_DP, BDAM= 0.90_DP
    REAL(DP), PARAMETER :: GRAV = 9.81_DP
    REAL(DP), PARAMETER :: es = 0.001_DP, e=2.302585093_DP
    REAL(DP), PARAMETER :: W0 = 1367.0_DP

    INTEGER(I4B), PARAMETER :: IMAX = 13
    INTEGER(I4B), PARAMETER ::HRSDAY = 24					!total hours per day
    !gp 30-Nov-04 INTEGER(I4B), PARAMETER :: NV = 16, NL = 2
    INTEGER(I4B), PARAMETER :: NV = 17, NL = 2		!gp 30-Nov-04 add generic constituent

    !gp 30-Jan-06 define null value for identifying missing reach-specific rates
    REAL(DP), PARAMETER :: NULL_VAL = -999_dp

    TYPE outputSummary_type
        REAL(DP) min, max, avg
    END TYPE outputSummary_type
END MODULE nrtype
