      program main
      implicit none
      integer:: i, count
      character(len=32) :: arg, fn
      real*8 ac(8), m, k, l(2,6,6)
      count = command_argument_count()

      if (count.ne.10) then
         write(*,*) 'The number of arguments should be 10'
         stop -1
      endif

      do i=1, 8
         call get_command_argument(i,arg)
         read(arg,*) ac(i)
      enddo
      call get_command_argument(9,arg)
      read(arg,*) m
      k = 2
      call calcyld(ac,m,k,l)

      !! save l to elsewhere?
      call get_command_argument(10,fn)
      open(1,file=fn,status='unknown',form='unformatted',
     $     access='sequential')
      write(1) l
      close(1)
c$$$      !! test reading
c$$$      open(1,file=fn,status='unknown',form='unformatted',
c$$$     $     access='sequential')
c$$$      read(1) l
c$$$      write(*,*) l
c$$$
      !! test tension tests.
      call tension2d(m,k,l,100)
      return
      end program main
