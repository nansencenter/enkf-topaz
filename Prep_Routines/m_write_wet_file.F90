module m_write_wet_file

contains 

  subroutine write_wet_file(obs, nrobs)
    use mod_measurement
    implicit none

    type (measurement), intent(in) :: obs(:)

    integer, intent(in):: nrobs
    integer j, i, nrshow
    logical ex 
    character*200 tmp_str

    nrshow = max(nrobs / 10, 1)
    print *, '10 observations:'
    tmp_str='    #    obs       var    id      lon   lat  depth   ipiv  jpiv   nsup 4-bilin-coeffs    active  orig (i,j)    N    age orig_id'
    print '(a)',trim(tmp_str)

    inquire(iolength = i) obs(1)
    open (11, file = 'observations.uf', status = 'replace',&
         form = 'unformatted', access = 'direct', recl = i)

    do j = 1, nrobs
       write(11, rec = j) obs(j) 
       if (obs(j) % d > 1.01 .and. trim(obs(j) % id) == 'ICEC') then
          print *, obs(j) % lon, obs(j) % lat, obs(j) % d, obs(j) % var
       end if
       if (mod(j, nrshow) == 0) then
          print '(i6,2g10.2,a6,3f6.1,3i6,4f5.1,l5,2i7,f7.1,i5,i8)', j, obs(j)
       end if
    enddo
    close(11)
    print *, 'Observations printed to file observation.uf'
  end subroutine write_wet_file

end module m_write_wet_file

