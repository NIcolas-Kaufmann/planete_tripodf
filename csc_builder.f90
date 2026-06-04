
! =============================================================================
! csc_from_coo.f90
!
! Module  : csc_builder
! Routine : coo_to_csc
!
! PURPOSE
!   Convert a sparse matrix from COO (Coordinate / triplet) format
!
!       (coo_val, coo_row, coo_col)   nnz_in entries, 1-based indices
!
!   to CSC (Compressed Sparse Column) format
!
!       col_ptr(1:ncol+1)   column pointer array (1-based)
!       row_ind(1:nnz_csc)  row index of each stored value
!       csc_val(1:nnz_csc)  stored values
!
!   Every pair of entries sharing the SAME (row, col) position are SUMMED
!   into a single stored value (standard FEM / matrix-assembly convention).
!
! ALGORITHM
!   1. Build a sort permutation ordering COO entries by (col ASC, row ASC)
!      via two-pass LSD radix / counting sort  ->  O(nnz + ncol + nrow).
!   2. Walk the sorted sequence once, accumulating duplicates in-place.
!   3. Build col_ptr from per-column unique-entry counts (prefix sum).
!
! INTERFACE
!   subroutine coo_to_csc(nrow, ncol, nnz_in,
!                         coo_val, coo_row, coo_col,
!                         col_ptr, row_ind, csc_val, nnz_csc)
!
!   IN
!     integer(ip) nrow, ncol        matrix dimensions
!     integer(ip) nnz_in            number of COO triplets (duplicates OK)
!     real(dp)    coo_val(nnz_in)   numerical values
!     integer(ip) coo_row(nnz_in)   1-based row indices    (1..nrow)
!     integer(ip) coo_col(nnz_in)   1-based column indices (1..ncol)
!
!   OUT
!     integer(ip) col_ptr(ncol+1)         CSC column pointers (1-based)
!     integer(ip) row_ind(:) allocatable  CSC row indices
!     real(dp)    csc_val(:) allocatable  CSC values (duplicates summed)
!     integer(ip) nnz_csc                 number of unique stored entries
!
! COMPILATION
!   gfortran -O2 -std=f2008 csc_from_coo.f90 -o test_csc && ./test_csc
! =============================================================================

module csc_builder
  use iso_fortran_env, only: dp => real64, ip => int32
  implicit none
  private
  public :: coo_to_csc

contains

subroutine coo_to_csc(nrow, ncol, nnz_in,       &
                        coo_val, coo_row, coo_col,  &
                        col_ptr, row_ind, csc_val, nnz_csc)
  ! ===========================================================================
    integer(ip), intent(in)  :: nrow, ncol, nnz_in
    real(dp),    intent(in)  :: coo_val(nnz_in)
    integer(ip), intent(in)  :: coo_row(nnz_in), coo_col(nnz_in)
 
    integer(ip),              intent(out) :: col_ptr(ncol + 1)
    integer(ip), allocatable, intent(out) :: row_ind(:)
    real(dp),    allocatable, intent(out) :: csc_val(:)
    integer(ip),              intent(out) :: nnz_csc
 
    integer(ip) :: perm(nnz_in)       ! sort permutation
    integer(ip) :: cnt(ncol)          ! unique entries per column
    integer(ip), allocatable :: tmp_row(:)
    real(dp),    allocatable :: tmp_val(:)
    integer(ip) :: i, j, k, c, r
 
    ! ------------------------------------------------------------------
    ! Step 1 – Build sort permutation: order COO entries (col ASC, row ASC).
    ! ------------------------------------------------------------------
    do i = 1, nnz_in
      perm(i) = i
    end do
    call sort_by_col_row(nnz_in, ncol, nrow, coo_col, coo_row, perm)
 
    ! ------------------------------------------------------------------
    ! Step 2 – Walk sorted sequence; accumulate duplicate (row,col) entries.
    !
    !   Two consecutive sorted entries at positions i-1 and i are a DUPLICATE
    !   when they share both column (the sort guarantees adjacency) AND row.
    !   Detection uses:
    !     coo_col(perm(i)) == coo_col(perm(i-1))   <- always true for dups
    !     coo_row(perm(i)) == tmp_row(k)            <- row of last unique entry
    ! ------------------------------------------------------------------
    allocate(tmp_row(nnz_in), tmp_val(nnz_in))
    cnt = 0
    k   = 0
 
    do i = 1, nnz_in
      c = coo_col(perm(i))
      r = coo_row(perm(i))
 
      ! NOTE: Fortran does NOT guarantee short-circuit evaluation of .and.
      ! operands, so we must use nested if-blocks to guard array accesses
      ! that would be out-of-bounds when k=0 (tmp_row) or i=1 (perm(i-1)).
      if (k > 0) then
        if (r == tmp_row(k) .and. c == coo_col(perm(i-1))) then
          ! Duplicate: add value to the running sum of the current unique entry.
          tmp_val(k) = tmp_val(k) + coo_val(perm(i))
        else
          ! New unique (row, col) pair.
          k          = k + 1
          tmp_row(k) = r
          tmp_val(k) = coo_val(perm(i))
          cnt(c)     = cnt(c) + 1
        end if
      else
        ! First entry ever — unconditionally open the first unique slot.
        k          = k + 1
        tmp_row(k) = r
        tmp_val(k) = coo_val(perm(i))
        cnt(c)     = cnt(c) + 1
      end if
    end do
    nnz_csc = k
 
    ! ------------------------------------------------------------------
    ! Step 3 – Build col_ptr: 1-based inclusive prefix sum of cnt.
    ! ------------------------------------------------------------------
    col_ptr(1) = 1
    do j = 1, ncol
      col_ptr(j + 1) = col_ptr(j) + cnt(j)
    end do
 
    ! ------------------------------------------------------------------
    ! Step 4 – Copy compacted workspace into allocatable output arrays.
    ! ------------------------------------------------------------------
    allocate(row_ind(nnz_csc), csc_val(nnz_csc))
    row_ind(1:nnz_csc) = tmp_row(1:nnz_csc)
    csc_val(1:nnz_csc) = tmp_val(1:nnz_csc)
 
    deallocate(tmp_row, tmp_val)
 
  end subroutine coo_to_csc


  ! ===========================================================================
  ! sort_by_col_row
  !
  ! Stably sort perm(1:nnz) so that
  !   ( col(perm(i)), row(perm(i)) )  is non-decreasing (col first, row second).
  !
  ! Method: LSD radix sort – two passes of stable counting sort.
  !   Pass 1: sort by row  (secondary key)
  !   Pass 2: sort by col  (primary key; stability preserves row order within
  !           each column bucket)
  !
  ! Complexity: O(nnz + ncol + nrow) time, O(nnz + max(ncol,nrow)) extra space.
  ! ===========================================================================
  subroutine sort_by_col_row(nnz, ncol, nrow, col, row, perm)
    integer(ip), intent(in)    :: nnz, ncol, nrow
    integer(ip), intent(in)    :: col(nnz), row(nnz)
    integer(ip), intent(inout) :: perm(nnz)   ! initialised 1..nnz on entry

    integer(ip) :: buf(nnz)
    integer(ip) :: i, key

    ! ------------------------------------------------------------------
    ! Pass 1: stable counting sort by row (secondary key)
    ! ------------------------------------------------------------------
    block
      integer(ip) :: off(nrow + 1), nxt(nrow)
      ! Count occurrences
      off = 0
      do i = 1, nnz
        off(row(i) + 1) = off(row(i) + 1) + 1
      end do
      ! Exclusive prefix sum -> start offset of each bucket (0-based position,
      ! i.e. off(1) = 0 means the first element of bucket 1 goes to index 1)
      off(1) = 1                          ! convert to 1-based write position
      do key = 2, nrow + 1
        off(key) = off(key) + off(key - 1)
      end do
      ! nxt(key) = next write position for bucket key
      nxt(1:nrow) = off(1:nrow)
      ! Scatter
      do i = 1, nnz
        key = row(perm(i))
        buf(nxt(key)) = perm(i)
        nxt(key) = nxt(key) + 1
      end do
      perm = buf
    end block

    ! ------------------------------------------------------------------
    ! Pass 2: stable counting sort by col (primary key)
    ! ------------------------------------------------------------------
    block
      integer(ip) :: off(ncol + 1), nxt(ncol)
      off = 0
      do i = 1, nnz
        off(col(i) + 1) = off(col(i) + 1) + 1
      end do
      off(1) = 1
      do key = 2, ncol + 1
        off(key) = off(key) + off(key - 1)
      end do
      nxt(1:ncol) = off(1:ncol)
      do i = 1, nnz
        key = col(perm(i))
        buf(nxt(key)) = perm(i)
        nxt(key) = nxt(key) + 1
      end do
      perm = buf
    end block

  end subroutine sort_by_col_row

end module csc_builder