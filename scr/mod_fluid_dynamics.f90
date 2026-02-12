module mod_fluid_dynamics
  !! Thin compatibility wrapper.
  !!
  !! The core implementation is in `module_find` (faithful translation of find.py).
  !! This module re-exports the main entry point with a name that matches the
  !! newer code organization.
  use module_find, only: compute_primitive_variables
  implicit none
  public :: compute_primitive_variables
end module mod_fluid_dynamics
