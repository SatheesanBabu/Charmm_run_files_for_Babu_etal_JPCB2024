c     POT_FIELD calculates the net electrostatic potential and field
c     at the solute sotes in TIP3P water molecules
c  
c     This program is applicable if the simulation is done
c     in a CUBIC periodic boundary of chosen radius
c
c************************************************************
c     Programmed by Dr. C.S. Babu
c     Structural biology program
c     Institute of Biomedical Sciences
c     Academia Sinica, Taipei 11529, Taiwan
c*************************************************************
c
c     The program reads a cor file with # of frames given in the
c     variable 'nfrmax'. The variable 'natmax' specifies the
c     total # of atoms in the PSF file.
c*************************************************************
c     VERSION 1.0
c*************************************************************
      parameter (nfrmax=3000,ishmax=100,natmax=81393,ncsmax=100,
     * nsite=1, maxbin=500)
c      parameter (indx=13944)
      implicit real*8 (a-h,o-z)
      character*4 ftype
      character*35 title
      character*80 dcdfile
      character*16 outfile
      character*16 crdfile
      character*16 potfile1
      character*16 potfile2
      character*16 potfile3
      character*16 potfile4
      character*16 potfile5
      character*16 potfile6
      character*16 potfile7
      character*16 psffile
c     in dcd files all the x,y,z are saved as real not real*8
      real*4 x(natmax),y(natmax),z(natmax)
      real*4 sumx(natmax),sumy(natmax),sumz(natmax)
      real*4 xavg(natmax),yavg(natmax),zavg(natmax)
      real corr,r75
        dimension xn(nsite),yn(nsite),zn(nsite)
        dimension atnum(natmax),resnum(natmax),sresnum(natmax)
        dimension atnum_fix(nsite),resnum_fix(nsite),cp(natmax)
        dimension hist_o(maxbin,nsite), hist_h(maxbin,nsite)
     * ,gro(maxbin,nsite),grh(maxbin,nsite),nro(maxbin,nsite)
     * ,nrh(maxbin,nsite),charge(nsite)
        dimension enro(maxbin,nsite),enrh(maxbin,nsite)
     * , enr_tot(maxbin,nsite)
        dimension costh(maxbin,nsite),ac(maxbin,nsite)
        dimension pot_o(nsite,maxbin),pot_h(nsite,maxbin),
     *    pot_tot(nsite,maxbin),pot_p(nsite,maxbin),
     *   pot_sum(nsite,maxbin),pot_cl(nsite,maxbin),
     *   pot_cla(nsite,maxbin), pot_cal(nsite,maxbin),
     *   pot_ion(nsite,maxbin)  
        dimension fx_o(nsite),fy_o(nsite),fz_o(nsite),fx_tot(nsite)
        dimension fx_h(nsite),fy_h(nsite),fz_h(nsite)
     * ,fy_tot(nsite),fz_tot(nsite)
        character*8  resid(natmax), atid(natmax), segid(natmax)
        character*8  resid_fix(nsite), atid_fix(nsite)
        character*8 solvid1,solvid2
        real *4 nideal_o, nideal_h,nro,nrh
        integer*4 atnum,resnum,sresnum,atnum_fix,resnum_fix
      integer icntrl(20),natom,numpnt,numframe,bin
      pi=acos(-1.d0)
C
      read (5,*) nstart,numframe,natomt,nfix
      natom=natomt-nfix
      if (natomt.gt.natmax) then
        write(6,*) 'ERROR: too many atoms'
        write(6,*) 'natomt must be less than ', natmax
        stop
      endif
      if (numframe.gt.nfrmax) then
        write(6,*) 'ERROR: too many time frames'
        write(6,*) 'numframe must be less than ', nfrmax
        stop
      endif
c
      print *, 'Total number of atoms = ',natomT,'
     >         , fixed atoms = ',nfix,', free atoms =',natom
      print *, 'first time frame = ',nstart,
     >         ', number of frames = ',numframe
c
c
c Read the coordinate file for residue records
        read(5,*)boxl
        write(6,*) 'boxl=', boxl
        read(5,*) crdfile
        read(5,*) potfile1
        open(1,file=crdfile, form='formatted')
        do i=1,5
        read(1,*)
        enddo
        do 120 i=1,natomt
        read(1,300)atnum(i),resnum(i),resid(i),atid(i)
     1 , x(i),y(i),z(i),segid(i),sresnum(i)
 120    continue
        close (1)
300    format(2I10,2(2x,A8),3(f20.10),2x,A8,2x,I8)
      read(5,*) psffile
      write(6,*) psffile
c
c..Define the number of sites you need to calculate the PAIRCF
      read(5,*)ns,delgr
      write(6,*)nsite,delgr

c   Now enter the segid of solvent and the # of solvent molecules
      print *, ' Read complete' 
 301  format(2I5,2(1x,A4))

c.. Now read the charges from PSF file
      open (22,file=psffile, form='formatted')
      do i=1,natomt
      read(22,401)cp(i)
      write(*,402)cp(i)
 401  format(51x,E18.10)
 402  format(F11.6)
      enddo

C
c
c.. Initialize the bins for pair correlation functions
       do i=1,maxbin
       do j=1,nsite
       hist_o(i,j)=0.0
       hist_h(i,j)=0.0
       enddo
	enddo
c.. Initialize the bins for potentials and fields
       do i = 1,nsite
       do j = 1,maxbin
       pot_o(i,j)=0.0
       pot_h(i,j)=0.0
       pot_p(i,j)=0.0
       pot_cla(i,j)=0.0
       pot_cal(i,j)=0.0
       enddo
         enddo
c.. tiph and tipo are charges on tip3p water
       tipo= -0.834
       tiph= 0.417
c.. Initialize the bins for pcos(theta)


       do i = 1,maxbin
       do j = 1, nsite
       costh(i,j)=0.0
         enddo
          enddo

c..Compute the constants for g(r)
c.. rho is the bulk solvent density = 0.033particles/A**3
       rho= 0.0334
       const_o= 4.0 * pi * rho /3.0
       const_h= 4.0 * pi * 2.0 * rho /3.0
       r75= 4.0/3.0

c     iterating on time frames
c.. open  the trajectory file

      open(20,file= "out_file")
       nline=0
       nf=0
 4444  nline=nline+1
c       write(*,*) nline
       read(20,98) title
c       write(*,98) title
  98   format (A35)
        if (title .eq. ' CHARMM>    PRINT COOR' ) then
        nf= nf +1
c        write(*,98) title
        iline=6
        do j=1,iline
        read(20,*)
        enddo

       do i=1,natomT
        read(20,300)atnum(i),resnum(i),resid(i),atid(i)
     1 , x(i),y(i),z(i),segid(i),sresnum(i)
        enddo


c.. Now calculate the distances of water molecules from the
c.. sites
c   Locate the M20 residues M20, N18, E17 and G 16
       do i=1,natomT
       if (segid(i) .eq. 'DHFR' ) then 
        if (resnum(i) .eq. 20 ) then
       if ( atid(i) .eq. 'CA' ) then
        x19= x(i)
        y19= y(i)
        z19= z(i)
           endif
            endif
        if (resnum(i) .eq. 18 ) then
        if ( atid(i) .eq. 'CA' ) then
        x29= x(i)
        y29= y(i)
        z29= z(i)
           endif
            endif
        if (resnum(i) .eq. 17 ) then
        if ( atid(i) .eq. 'CA' ) then
        x39= x(i)
        y39= y(i)
        z39= z(i)
           endif
            endif 
        if (resnum(i) .eq. 16 ) then
        if ( atid(i) .eq. 'CA' ) then
        x49= x(i)
        y49= y(i)
        z49= z(i)
           endif
              endif
                endif

                enddo
c...    Average coordivates
         xav =  (x19 + x29 + x39 + x49) / 4.0
         yav =  (y19 + y29 + y39 + y49) / 4.0
         zav =  (z19 + z29 + z39 + z49) / 4.0

        do 190 j=1,maxbin
        dr= float(j-1)*delgr
        do 181 i= 1,natomT
 
c.. Compute the center of S-S briges in cluster



c        if (segid(i).eq. solvid ) then
        do 182 k=1,nsite
c
c... Compute the solvent contributions
         if ( atid(i) .eq. 'OH2' ) then
c  For site-O distribution
        rx = xav - x(i)
        ry = yav - y(i)
        rz = zav - z(i)
c..divide by blen
        rx=rx/boxl
        ry=ry/boxl
        rz=rz/boxl
        rx = rx -  anint (rx)
        ry = ry -  anint (ry)
        rz = rz -  anint (rz)
c..multiply back bu blen
        rx= rx*boxl
        ry= ry*boxl
        rz= rz*boxl
        rss= sqrt(rx**2 + ry**2 + rz**2 )
        if (rss. le. dr ) then

c.. NOW COMPUTE THE POTENTIALS AND FIELDS
        pot_o(k,j)= pot_o(k,j) + tipo*332.0/rss
         endif
          endif
c           endif           
       if ( atid(i) .eq. 'H1' .or. atid(i) .eq. 'H2' ) then
c
c  For site-H distribution
        rx = xav - x(i)
        ry = yav - y(i)
        rz = zav - z(i)
c..divide by blen
        rx=rx/boxl
        ry=ry/boxl
        rz=rz/boxl
        rx = rx - anint (rx)
        ry = ry - anint (ry)
        rz = rz - anint (rz)
c.. correct for the toc
c..multiply back by blen
        rx= rx*boxl
        ry= ry*boxl
        rz= rz*boxl
        rss= sqrt(rx**2 + ry**2 + rz**2 )
        if (rss. le. dr ) then

c.. NOW COMPUTE THE POTENTIALS AND FIELDS
        pot_h(k,j)= pot_h(k,j) + tiph*332.0/rss
         endif
          endif

c... Now for IONS CLA and CAL
       if ( atid(i) .eq. 'CLA' ) then
c
c  For site-H distribution
        rx = xav - x(i)
        ry = yav - y(i)
        rz = zav - z(i)
c..divide by blen
        rx=rx/boxl
        ry=ry/boxl
        rz=rz/boxl
        rx = rx - anint (rx)
        ry = ry - anint (ry)
        rz = rz - anint (rz)
c.. correct for the toc
c..multiply back by blen
        rx= rx*boxl
        ry= ry*boxl
        rz= rz*boxl
        rss= sqrt(rx**2 + ry**2 + rz**2 )
        if (rss. le. dr ) then

c.. NOW COMPUTE THE POTENTIALS AND FIELDS
        pot_cla(k,j)= pot_cla(k,j) + (-1.0*332.0)/rss
         endif
          endif

       if ( atid(i) .eq. 'CAL' ) then
c
c  For site-H distribution
        rx = xav - x(i)
        ry = yav - y(i)
        rz = zav - z(i)
c..divide by blen
        rx=rx/boxl
        ry=ry/boxl
        rz=rz/boxl
        rx = rx - anint (rx)
        ry = ry - anint (ry)
        rz = rz - anint (rz)
c.. correct for the toc
c..multiply back by blen
        rx= rx*boxl
        ry= ry*boxl
        rz= rz*boxl
        rss= sqrt(rx**2 + ry**2 + rz**2 )
        if (rss. le. dr ) then

c.. NOW COMPUTE THE POTENTIALS AND FIELDS
        pot_cal(k,j)= pot_cal(k,j) + ( 2.0*332.0 )/rss
         endif
          endif


c..Now ion_protein interaction
        if (segid(i) .eq. 'DHFR'  )  then
        rx = xav - x(i)
        ry = yav - y(i)
        rz = zav - z(i)
c..divide by blen
        rx=rx/boxl
        ry=ry/boxl
        rz=rz/boxl
        rx = rx - anint (rx)
        ry = ry - anint (ry)
        rz = rz - anint (rz)
c.. correct for the toc
c..multiply back by blen
        rx= rx*boxl
        ry= ry*boxl
        rz= rz*boxl
        rss= sqrt(rx**2 + ry**2 + rz**2 )

c.. to escape CA s of 20,18,17, 16
        if (rss .eq. 0.0 ) then
         go to 1182
          endif

        if (rss .le. dr ) then
c.. NOW COMPUTE THE POTENTIALS AND FIELDS
        pot_p(k,j)= pot_p(k,j) + cp(i)*332.0/rss
         endif
1182          endif
c               endif
c                endif


  182   continue
c           endif
  181    continue
  190     continue
         write(6,*) nf
           endif
        if (nf .eq. numframe ) then
          go to 180
            endif
        if (title .ne. ' CHARMM>    STOP') then
          go to 4444
             else
            endif
  180     continue
         write(6,*) nf
         if (nf .ne. numframe) then
         write (*,*) 'error in numframe' 
         endif
c.. Add the O and H potential and field contributions at sites
        open(15, file=potfile1 , form='formatted')
 700    format(6f20.5)
c        write(15,*) ' pot_o   pot_h   pot_tot'
        do j= 1,maxbin
        do i = 1, nsite
        pot_o(i,j)=pot_o(i,j)/real(numframe) 
        pot_h(i,j)=pot_h(i,j)/real(numframe) 
        pot_cla(i,j) = pot_cla(i,j)/real(numframe)
        pot_cal(i,j) = pot_cal(i,j)/real(numframe)
        pot_ion(i,j) = pot_cla(i,j) + pot_cal(i,j)
        pot_tot(i,j) = pot_o(i,j) + pot_h(i,j)
        pot_p(i,j)= pot_p(i,j)/real(numframe)
        pot_sum(i,j)=pot_tot(i,j)+pot_p(i,j)+ pot_ion(i,j)
        enddo
        enddo
        do j = 1,maxbin
        dr = float(j-1)*delgr
       write(15,700)dr,pot_tot(1,j),pot_p(1,j),pot_ion(1,j),pot_sum(1,j)
         enddo
c.. Now write the distribution functions

      stop
      END
