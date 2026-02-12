module module_mesh
  !! Mesh generation utilities (cell-centered).
  !!
  !! Progressive geometric mesh with refinement around BASE (s = Lp):
  !!  - Pipe   : [0, Lp] with sizes decreasing towards Lp
  !!  - Riser  : [Lp, Lp+Lr] with sizes increasing away from Lp
  !!
  !! Output coordinates are non-dimensionalized by Lr (same convention as main).

  use mod_kinds, only: wp, i4
  implicit none
  private
  public :: build_progressive_mesh

contains

  pure real(wp) function geom_sum(r, n) result(s)
    real(wp), intent(in) :: r
    integer(i4), intent(in) :: n
    if (n <= 0) then
      s = 0.0_wp
    else if (abs(r - 1.0_wp) < 1.0e-14_wp) then
      s = real(n, wp)
    else
      s = (r**real(n, wp) - 1.0_wp) / (r - 1.0_wp)
    end if
  end function geom_sum

  subroutine split_cells_sym_base(nCells, Lp, Lr, r, Np, Nr)
    !! Choose (Np,Nr) so that the **minimum** cell size adjacent to the BASE
    !! is as symmetric as possible on both sides.
    !!
    !! Why this exists:
    !! If we split cells only by length (Lp/(Lp+Lr)), then with a geometric
    !! ratio r>1 the pipe (many cells) becomes extremely refined near BASE,
    !! while the riser (few cells) stays coarse. That produces a large jump
    !! in theta(s) right after BASE.
    !!
    !! Here we pick Np,Nr such that:
    !!   dx0_pipe  = Lp / geom_sum(r, Np)
    !!   dx0_riser = Lr / geom_sum(r, Nr)
    !! are as close as possible (dx0_pipe ~= dx0_riser).
    integer(i4), intent(in)  :: nCells
    real(wp),    intent(in)  :: Lp, Lr, r
    integer(i4), intent(out) :: Np, Nr

    integer(i4) :: Np_try, Nr_try
    real(wp)    :: dxp, dxr, obj, obj_best

    if (nCells < 4) then
      Np = max(1, nCells-1)
      Nr = nCells - Np
      return
    end if

    obj_best = huge(1.0_wp)
    Np = max(2, min(nCells-2, nint(real(nCells, wp) * (Lp / max(Lp+Lr, 1.0e-30_wp)) )))
    Nr = nCells - Np

    do Np_try = 2, nCells-2
      Nr_try = nCells - Np_try
      if (Nr_try < 2) cycle

      dxp = Lp / max(geom_sum(r, Np_try), 1.0e-30_wp)
      dxr = Lr / max(geom_sum(r, Nr_try), 1.0e-30_wp)

      ! symmetric refinement around BASE -> dxp ~ dxr
      obj = abs(log(max(dxp,1.0e-30_wp)/max(dxr,1.0e-30_wp)))

      if (obj < obj_best) then
        obj_best = obj
        Np = Np_try
        Nr = Nr_try
      end if
    end do
  end subroutine split_cells_sym_base

  subroutine build_progressive_mesh(nCells, Lp, Lr, mesh_ratio, Np, Nr, s_faces_nd, s_centers_nd)
    integer(i4), intent(in)  :: nCells
    real(wp),    intent(in)  :: Lp, Lr, mesh_ratio
    integer(i4), intent(out) :: Np, Nr
    real(wp), allocatable, intent(out) :: s_faces_nd(:), s_centers_nd(:)

    integer(i4) :: i
    real(wp) :: Ltot, r
    real(wp) :: dx0_pipe, dx0_riser
    real(wp), allocatable :: faces_dim(:)

    r = max(mesh_ratio, 1.0_wp)   ! ensure >=1
    Ltot = Lp + Lr

    ! --- Split cells to make refinement around BASE symmetric (pipe vs riser) ---
    call split_cells_sym_base(nCells, Lp, Lr, r, Np, Nr)

    allocate(faces_dim(nCells+1))
    faces_dim = 0.0_wp

    ! --- PIPE faces: smallest cell next to BASE ---
    dx0_pipe  = Lp / max( geom_sum(r, Np), 1.0e-30_wp )

    faces_dim(1) = 0.0_wp
    do i = 1, Np
      faces_dim(i+1) = faces_dim(i) + dx0_pipe * r**real(Np-i, wp)
    end do
    faces_dim(Np+1) = Lp

    ! --- RISER faces: smallest cell next to BASE ---
    dx0_riser = Lr / max( geom_sum(r, Nr), 1.0e-30_wp )
    do i = 1, Nr
      faces_dim(Np+1+i) = faces_dim(Np+i) + dx0_riser * r**real(i-1, wp)
    end do
    faces_dim(nCells+1) = Ltot

    ! --- Non-dimensionalize by Lr ---
    allocate(s_faces_nd(nCells+1), s_centers_nd(nCells))
    do i = 1, nCells+1
      s_faces_nd(i) = faces_dim(i) / max(Lr, 1.0e-30_wp)
    end do

    ! centers (cell-centered)
    do i = 1, nCells
      s_centers_nd(i) = 0.5_wp * (faces_dim(i) + faces_dim(i+1)) / max(Lr, 1.0e-30_wp)
    end do

    deallocate(faces_dim)
  end subroutine build_progressive_mesh

end module module_mesh
