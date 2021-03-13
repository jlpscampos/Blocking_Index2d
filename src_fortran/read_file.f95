PROGRAM MAIN
    IMPLICIT NONE

    CHARACTER(100) :: input_dir, output_dir, conf_file
    INTEGER        :: begin, final

    conf_file =  'config.fortran.txt'
    CALL read_file( conf_file, input_dir, output_dir, begin, final)
    
    WRITE(*,*) input_dir
    WRITE(*,*) output_dir
    WRITE(*,*) begin
    WRITE(*,*) final
    
    CONTAINS

       SUBROUTINE read_file( conf_file, input_dir, output_dir, begin, final)

          IMPLICIT NONE
          
          CHARACTER(100), INTENT(IN)  :: conf_file
          CHARACTER(100), INTENT(INOUT) :: input_dir, output_dir
          INTEGER       , INTENT(INOUT) :: begin, final

          OPEN(10, file=conf_file )
          
          READ(10,'(a)') input_dir
          READ(10,'(a)') output_dir
          READ(10,'(i4)') begin
          READ(10,'(i4)') final
    
          CLOSE(10)

       END SUBROUTINE

END PROGRAM
