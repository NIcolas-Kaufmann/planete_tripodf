subroutine OPEN_OUTPUT_FILE(desc, lines, keep, truncate, base, suf_len, suf_val)

! I/O variables
integer, intent(in) :: desc
integer, intent(inout) :: lines
logical, intent(in) :: keep
logical, intent(in) :: truncate
character(len=*), intent(in) :: base
integer, intent(in) :: suf_len
integer, intent(in) :: suf_val

! Local variables
logical :: opened
integer :: iline
character(len=1) :: num
character(len=256) :: name


inquire(unit = desc, opened = opened)

if (opened) then
   ! File already open
   if (.not. keep) then
      ! Reset to start
      rewind(unit = desc)
      lines = 0
   end if

   return
end if

! Get the file name
if (suf_len == 0 .or. suf_val == 0) then
   name = base // ".outputdat"
else
   write(num, "(i1.1)") suf_len
   write(name,"(a,a,i" // num // "." // num // ",a)") base, "_", suf_val, ".outputdat"
end if

if (lines == 0 .or. .not.keep) then
   ! Restart from the beginning
   open(unit = desc, file = name, position = "rewind", action = "write")
   lines = 0
else if (.not. truncate) then
   ! Append to the file without checking its length
   open(unit = desc, file = name, position = "append", action = "write")
else
   ! Append at the correct location
   open(unit = desc, file = name, position="rewind", action="readwrite")

   ! Position after the lines that should be kept
   do iline = 1, lines
      read(desc, *)
   end do
end if

end


subroutine read_csv(filename, nrows, ncols, data_out, ierr)
  ! Reads a CSV file and fills a 2D array with the entries.
  !
  ! Parameters
  ! ----------
  ! filename : path to the CSV file
  ! nrows    : number of data rows    (excluding header)
  ! ncols    : number of columns
  !
  ! Output
  ! ------
  ! data_out : array of values read   (nrows, ncols)
  ! ierr     : 0 on success, non-zero on failure

  implicit none

  character(len=*), intent(in)  :: filename
  integer,          intent(in)  :: nrows, ncols
  real(8),          intent(out) :: data_out(nrows, ncols)
  integer,          intent(out) :: ierr

  ! Local variables
  integer            :: unit, i, j, ios
  character(len=512) :: line
  character(len=64)  :: token
  integer            :: pos, next, col

  ierr = 0
  unit = 20

  open(unit=unit, file=trim(filename), status='old', action='read', iostat=ios)
  if (ios /= 0) then
    write(*,*) 'ERROR: could not open file: ', trim(filename)
    ierr = 1
    return
  end if

  ! ----------------------------------------------------------------
  ! Skip header line
  ! ----------------------------------------------------------------
  !read(unit, '(A)', iostat=ios) line
  !if (ios /= 0) then
  !  write(*,*) 'ERROR: could not read header line'
  !  ierr = 2
  !  close(unit)
  !  return
  !end if

  ! ----------------------------------------------------------------
  ! Read data rows
  ! ----------------------------------------------------------------
  do i = 1, nrows
    read(unit, '(A)', iostat=ios) line
    if (ios /= 0) then
      write(*,'(A,I0)') 'ERROR: failed reading row ', i
      ierr = 3
      close(unit)
      return
    end if

    ! Parse comma-separated tokens from the line
    pos = 1
    col = 1
    do while (col <= ncols)
      ! Find next comma or end of string
      next = index(line(pos:), ',')

      if (next == 0) then
        ! No more commas: take the rest of the line as last token
        token = adjustl(trim(line(pos:)))
      else
        token = adjustl(trim(line(pos : pos + next - 2)))
        pos   = pos + next
      end if

      read(token, *, iostat=ios) data_out(i, col)
      if (ios /= 0) then
        write(*,'(A,I0,A,I0)') 'ERROR: could not parse value at row ', i, ', col ', col
        ierr = 4
        close(unit)
        return
      end if

      col = col + 1
      if (next == 0) exit   ! end of line reached
    end do

  end do

  close(unit)

end subroutine read_csv