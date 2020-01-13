
      subroutine radialdistribution(control,NP,rx,ry,rz,nbitac,iout)
      implicit none
      include 'params_rdf.inc' !pi 
      integer control,i,j,bin,maxbin,ISTEP,nbitac,istp,iout
      parameter (maxbin=250)
      INTEGER HIST(MAXbin)
      real GR(maxBIN),DELT_R,CONST,RUP,RLOW,NIDEAL      
      real rxij,ryij,rzij,rij,DENS
      SAVE ISTEP,ISTP,DENS,DELT_R,HIST,GR
      real pi
      save pi
      integer np
      real rx(npmx),ry(npmx),rz(npmx)



      
c     INISIALIZATION ---->CONTROL =0
      IF(CONTROL.EQ.0) THEN
         pi=4.0*atan(1.0)

         open(iout,file='rdf.res')
         DELT_R=LX/REAL(MAXBIN)/2.0
         ISTEP=0
         ISTP=0
         DO BIN=1,MAXBIN
            HIST(BIN)=0
            GR(BIN)=0.0
         END DO
         DENS=REAL(NP)/(LX*LY*LZ)

c     STORE DATA----CONTROL=1

         ELSE IF(CONTROL.EQ.1) THEN

         ISTEP=ISTEP+1


         IF(MOD(ISTEP,NBITAC).eq.0) THEN
            ISTP=ISTP+1

            WRITE(*,*) 'ISTEP---RADIAL',ISTP
            do i=1,nP-1
               do j=I+1,NP
                  
                  rxij=rx(i)-rx(j)
                  ryij=ry(i)-ry(j)
                  rzij=rz(i)-rz(j) 

               RXIJ=RXIJ-ANINT(rXIJ/LX)*LX
               RyIJ=RyIJ-ANINT(ryIJ/Ly)*LY
               RzIJ=RzIJ-ANINT(rzIJ/Lz)*LZ
            
            rij=rxij*rxij+ryij*ryij+rzij*rzij
            Rij=sqrt(RIJ)

            BIN=INT(RIJ/DELT_R)+1

            IF(BIN.LE.MAXBIN) THEN
               HIST(BIN)=HIST(BIN)+2
            END IF

            end do
         end do
      END IF

c  FOR CONTROL=2 ----> NORMALIZATION 

         ELSE IF(CONTROL.EQ.2) THEN
         CONST=4.0*PI*DENS/3.0

         DO BIN=1,MAXBIN

            RLOW=REAL(BIN-1)*DELT_R
            RUP=RLOW+DELT_R
            NIDEAL=CONST*(RUP**3-RLOW**3)
            GR(BIN)=REAL(HIST(BIN))/REAL(ISTP)/REAL(NP)/NIDEAL

            WRITE(IOUT,*) (RUP+RLOW)/2.,GR(BIN)

         end do

         close(iout)

      END IF

      RETURN
      END
