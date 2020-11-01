! -----------------------------------------------------------
!     Module for using the spatial flux variation method
! -----------------------------------------------------------
!
!                    Public subroutines
!
! predict_spatial_flux -> main subroutine that accepts data at
!   known and unknown points, applies the perturbation theory
!   approach and predicts the flux distribution at the unknown
!   point.
!
! -----------------------------------------------------------
module sfvmod
   implicit none
   integer, parameter :: dp = kind(1.d0)
   real(dp), parameter :: unity = 1.0_dp
   integer, parameter :: sfv_v_major=0
   integer, parameter :: sfv_v_minor=3
   integer, parameter,dimension(2) :: sfv_v = (/ sfv_v_major, sfv_v_minor /)
!f2py integer,parameter,dimension(2) :: sfv_v
!f2py integer,parameter :: sfv_v_major
!f2py integer,parameter :: sfv_v_minor

   public :: sfv_v_major, sfv_v_minor, sfv_v
   public :: predict_spatial_flux
   private :: unity
   private :: get_expansion_coeffs
   private :: compute_delta_lambda_fopt

contains

   ! ------------------------------------------------------------
   ! subroutine predict_spatial_flux
   ! ------------------------------------------------------------
   ! Computes the change in spatial flux distribution
   ! and applies it onto phi0
   ! Values are taken between two states: a known and predicted
   ! state
   ! ------------------------------------------------------------
   ! nMat [integer, in] - number of materials or nodes to be
   !   considered.
   ! nMode [integer, in] - number of modes of the flux in the
   !   forward and adjoint vectors
   ! abs0 [double (nMat), in] - macroscopic absorption cross
   !   section in each node for the known state
   ! abs1 [double (nMat), in] - macroscopic absorption cross
   !   section in each node for the predicted / unknown state
   ! nsf0 [double (nMat), in] - macroscopic nu-sigma fission
   !   cross section in each node for the known state
   ! nsf1 [double (nMat), in] - macroscopic nu-sigma fission
   !   cross section in each node for the unknown state
   ! keff0 [double, in] - multiplication factor for the known state
   ! keff1 [double, in] - multiplication factor for the unknown
   !   state. If less than zero, then use first order perturbation
   !   theory and macroscopic cross sections to predict the
   !   change in reactivity and compute keff1
   ! adjFxMom [double (nMat, nMode), in] - adjoint flux modes in
   !   each node. adjFxMom(:, i) is the i-th left eigenvector
   !   that corresponds to eigW(i)
   ! fwdFxMom [double (nMat, nMode), in] - forward flux modes in
   !   each node. fwdFxMom(:, i) is the i-th right eigenvector
   !   that corresponds to eigW(i)
   ! eigW [double (nMode), in] - k-eigenvalues corresponding to
   !   each eigenvector eigW(i) in fwdFxMom(:, i) and adjFxMom(:, i)
   ! phi0 [double (nMode), inout] - Flux at the known point in each
   !   node. Will be overwritten with the prediction of the flux
   !   at the unknown point
   ! ------------------------------------------------------------
   subroutine predict_spatial_flux(nMat, nMode, abs0, abs1, &
                                   nsf0, nsf1, keff0, keff1, &
                                   adjFxMom, fwdFxMom, &
                                   eigW, phi0)
      implicit none
      integer, intent(in) :: nMat, nMode
      real(dp), intent(in) :: abs0(nMat), abs1(nMat)
      real(dp), intent(in) :: nsf0(nMat), nsf1(nMat)
      real(dp), intent(in) :: keff0, keff1
      real(dp), intent(in) :: adjFxMom(nMat, nMode)
      real(dp), intent(in) :: fwdFxMom(nMat, nMode)
      real(dp), intent(in) :: eigW(nMode)
      real(dp), intent(inout) :: phi0(nMat)
!f2py integer,intent(hide),depend(abs0) :: nMat=len(abs0)
!f2py integer,intent(hide),depend(eigW) :: nMode=len(eigW)
!f2py real(dp),optional :: keff1

      integer :: j
      real(dp) :: keff_inv
      real(dp) :: dsiga(nMat), dnsf(nMat)
      real(dp) :: deltaLambda, expCoeffs(nMode - 1)

      keff_inv = unity/keff0

      ! compute differences in cross sections
      dsiga = abs1 - abs0
      dnsf = nsf1 - nsf0

      if (keff1 .le. 0.0_dp) then

         ! compute delta lambda using first order perturbation theory
         call compute_delta_lambda_fopt(nMat, keff_inv, fwdFxMom(:, 1), &
                                        phi0, abs1, dSiga, &
                                        dnsf, deltaLambda)

      else

         deltaLambda = (keff0 - keff1)/(keff0*keff1)

      endif

      call get_expansion_coeffs(nMat, nMode, fwdFxMom, &
                                adjFxMom, keff_inv, dSiga, &
                                nsf0, dnsf, deltaLambda, &
                                eigW, phi0, expCoeffs)
      do j = 2, nMode
         call daxpy(nMat, expCoeffs(j - 1), fwdFxMom(:, j), 1, &
                    phi0, 1)
      end do
      call dscal(nMat, unity/sum(phi0), phi0, 1)

   end subroutine

   ! Compute the change in lambda eigenvalue using first
   ! order perturbation theory
   subroutine compute_delta_lambda_fopt(nMat, lambda0, funAdjMom, &
                                        flux0, abs1, dsiga, &
                                        dnsf, deltaLambda)
      integer, intent(in) :: nMat
      real(dp), intent(in) :: lambda0
      real(dp), intent(in) :: funAdjMom(nMat), flux0(nMat)
      real(dp), intent(in) :: abs1(nMat)
      real(dp), intent(in) :: dsiga(nMat), dnsf(nMat)
      real(dp), intent(out) :: deltaLambda

      real(dp) :: work(nMat)
      real(dp) :: ddot

      call dcopy(nMat, dnsf, 1, work, 1)
      call daxpy(nMat, -lambda0, dsiga, 1, work, 1)
      deltaLambda = ddot(nMat, funAdjMom, 1, work * flux0, 1) / &
         ddot(nMat, funAdjMom, 1, abs1 * flux0, 1)
   end subroutine

   ! ------------------------------------------------------------
   ! Routine for computing the expansion coefficents needed to
   ! compute the change in flux
   ! ------------------------------------------------------------
   ! nMats [int, in] - Number of materials / nodes to be considered
   ! nModes [int, in] - Number of modes of the forward and adjoint
   !   flux that are coming in
   ! fwdFlx0 [double (nMats, nModes), in] - Moments of the forward
   !   flux at the known
   ! adjFlx0 [double (nMats, nModes), in] - Moments of the adjoint
   !   flux at the known
   ! inv_keff0 [double, in] - Lambda {1 / k_eff} eigenvalue of the
   !   system
   ! dsiga [double (nMats), in] - Change in macroscopic absorption
   !   cross section in each node from the unknown point to the
   !   known point (siga_1 - siga_1)
   ! nsf0 [double (nMats), in] - nu-sigma fission cross section
   !   in each node at the known point
   ! dnsf [double (nMats), in] - Change in macroscopic nu-sigma
   !   fission cross section in each node from the unknown point
   !    to the known point (siga_1 - siga_1)
   ! deltaLambda [double, in] - Change in lambda eigenvalues
   !   between the unknown and known point
   ! eig0 [double (nModes), in] - Eigenvalues from the same
   !   eigensystem as fwdFlx0 and adjFlx0
   ! phi0 [double (nMats), in] - Forward flux at the known point
   ! expCoeffs [double (nMats - 1), in] - Expansion coefficients
   !   a_m to compute the change in flux distribution
   ! ------------------------------------------------------------
   subroutine get_expansion_coeffs(nMats, nModes, fwdFlx0, &
                                   adjFlx0, inv_keff0, dsiga, &
                                   nsf0, dnsf, &
                                   deltaLambda, &
                                   eig0, phi0, expCoeffs)
      implicit none
      integer, intent(in) :: nMats, nModes
      real(dp), intent(in) :: inv_keff0, deltaLambda
      real(dp), intent(in), dimension(nMats) :: dsiga, &
                                                nsf0, dnsf
      real(dp), intent(in) :: eig0(nModes)
      real(dp), intent(in), dimension(nMats, nModes) :: &
         fwdFlx0, adjFlx0
      real(dp), intent(in), dimension(nMats) :: phi0
      real(dp), intent(out) :: expCoeffs(nModes - 1)

      real(dp) :: fundEig, deltaEig, t0, t1, denom
      integer :: modeIndex, matIndex
      ! temp vectors
      real(dp) :: dnsf_o_k(nMats)

      fundEig = eig0(1)

      call dcopy(nMats, dnsf, 1, dnsf_o_k, 1)
      call dscal(nMats, inv_keff0, dnsf_o_k, 1)

      do modeIndex = 2, nModes
         deltaEig = (eig0(modeIndex) - fundEig)/ &
                    (eig0(modeIndex)*fundEig)
         t0 = 0.0_dp
         t1 = 0.0_dp
         denom = 0.0_dp
         do matIndex = 1, nMats
            t0 = t0 + (adjFlx0(matIndex, modeIndex)* &
                       (dsiga(matIndex) - dnsf_o_k(matIndex))* &
                       phi0(matIndex))
            t1 = t1 + (adjFlx0(matIndex, modeIndex)* &
                       nsf0(matIndex)*phi0(matIndex))
            denom = denom + (nsf0(matIndex)* &
                             adjFlx0(matIndex, modeIndex)* &
                             fwdFlx0(matIndex, modeIndex))
         end do
         expCoeffs(modeIndex - 1) = (t0 - deltaLambda*t1) &
                                    /(denom*deltaEig)
      end do
   end subroutine
end module
