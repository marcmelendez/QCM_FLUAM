      integer npmx
      parameter (npmx=5000)
      real x(npmx),y(npmx),z(npmx),sx
      real xd(npmx),yd(npmx),zd(npmx)
      character*10 TEXT
      real lbox(3)
      real dum
      integer ii,ic(npmx),istp,iout,iz
      integer diameter !particle diameter

      integer np
      real lx,ly,lz
      common/geom/lx,ly,lz

      open(47,file='data.punto')
      read(47,*) np
      read(47,*) lbox
      read(47,*) diameter
      close(47)
      if(np.gt.npmx) stop
      lx=lbox(1)
      ly=lbox(2)
      lz=lbox(3)

      read(*,*,end=100) TEXT, TEXT, NP
      write(*,*) NP

      nbitac=1
      iout=33

      call radialdistribution(0,NP,x,y,z,nbitac,iout)
      ii=0
      istp=0

 10   continue
      ii=ii+1
      read(*,*,end=100) time
      do i=1,np
         read(*,*,end=100) x(i),y(i),z(i)
         x(i)=x(i)-anint(x(i)/lbox(1))*lbox(1)
         y(i)=y(i)-anint(y(i)/lbox(2))*lbox(2)
         z(i)=z(i)-anint(z(i)/lbox(3))*lbox(3)
c     THIS IS TO DRAW PARTICLES WITH BOX IN PUNTo,  size and color
c     USING punto -r -c ..
         write(21,*) x(i),y(i),z(i),diameter
      end do
      call radialdistribution(1,NP,x,y,z,nbitac,iout)

      write(21,*) '#',time
      goto 10


 100  continue

      call radialdistribution(2,NP,x,y,z,nbitac,iout)
      end
