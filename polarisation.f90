      program polarisation_relaxation
      !include "mpif.h"  
      implicit none 
      include "mpif.h"
      real(8)::f,l,f0,e,r0,rt,r,summ,lsum
      integer::nn,np,nsi,i,j,k,fn,cnp
      integer,allocatable:: st(:,:,:),stold(:,:,:)
      real(8):: eb,nbe,dp 
      real(8),parameter::kt=0.026,v0=3,d=0.6 ! v0=3E+12=3Thz=3000/ns
      integer::ns 
      integer :: myid,ierr,npcs,status(MPI_STATUS_SIZE) 

      call random_seed ()
      call MPI_INIT( ierr )     
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )     
      call MPI_COMM_SIZE( MPI_COMM_WORLD, npcs, ierr )
      call init_random_seed(myid)


       open(101,file="in",status="old") 
      read(101,*) eb, nbe !  barrier, neighbouring 0.2/4
      read(101,*) dp 
      read(101,*) f0,cnp 
      read(101,*) ns 
      close(101) 
!    write(*,*) "begin allocate"
    allocate(st(cnp+2,cnp+2,ns),stold(cnp+2,cnp+2,ns))
!     write(*,*) eb,nbe
!     write(*,*) dp
!     write(*,*) f0 , cnp 
!     write(*,*) "simulation times",ns 
 
      f0=f0/300 ! (np*0.6)
      r0=exp(2*dp*f0/kt) ! np/nn
!      write(*,"(a,xf10.5)") "The initial ratio is", r0-1 
      
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Set initial states
!!!
do nsi=myid+1,ns,npcs
  do i=2,cnp+1
      do j=2,cnp+1
       CALL RANDOM_NUMBER(r)
       if (r <= r0/(1+r0)) stold(i,j,nsi)=1
       if (r > r0/(1+r0)) stold(i,j,nsi)=0 
       end do 
   end do 
   stold(1,:,nsi)   =stold(cnp+1,:,nsi)
   stold(:,1,nsi)   =stold(:,cnp+1,nsi)
   stold(cnp+2,:,nsi)=stold(2,:,nsi)
   stold(:,cnp+2,nsi)=stold(:,2,nsi)
end do 
   
!      write(*,"(a)") "Initialized"
! begin time evolvement 
! time step is 1 ns 
st=0
   do k=1,100
       lsum=0 
       do nsi=myid+1,ns,npcs 
       nn=0; np=0
       do i=2,cnp+1
           do j=2,cnp+1 
             e=nbe*(abs(stold(i,j,nsi)-stold(i-1,j,nsi))+abs(stold(i,j,nsi)-stold(i+1,j,nsi)) & 
             & +abs(stold(i,j,nsi)-stold(i,j-1,nsi))+abs(stold(i,j,nsi)-stold(i,j+1,nsi))) 
             
             rt=exp((-eb+e)/kt)  ! hopping rate 
!             rt=rhop*rt   !possibility 
             CALL RANDOM_NUMBER(r)
 !            if (r <= rt .and. 1/(rt*v0) <= 1) then 
          if (r <= rt ) then 
                st(i,j,nsi)=abs(stold(i,j,nsi)-1) 
                fn=fn+1
 !              write(*,"(i0,xi0)") stold(i,j) , st(i,j)
             else
               st(i,j,nsi)=stold(i,j,nsi)
!                write(*,*) "Not flipped"
             end if
             if (st(i,j,nsi)==1) np=np+1
              if (st(i,j,nsi)==0) nn=nn+1
            end do
       end do 
     st(1,:,nsi)   =st(cnp+1,:,nsi)
     st(:,1,nsi)   =st(:,cnp+1,nsi)
     st(cnp+2,:,nsi)=st(2,:,nsi)
     st(:,cnp+2,nsi)=st(:,2,nsi)
!       if (np/nn <= exp(f*d/kt)) exit 
       stold(:,:,nsi)=st(:,:,nsi)
!       if (mod(k,1000)==0)       
!       write(*,"(a,xi0,xi0,xi0,xf10.5)") "flipped",fn, np,nn,1.0*np/nn 
      lsum=lsum+(1.0*np/nn-1)/(r0-1)
     end do ! nsi 
!     write(*,*) "before reduce",myid,lsum
     call  MPI_Reduce(lsum,summ,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr);
!     write(*,*) "end reduce"
!call MPI_ALLREDUCE(tstd,ttstd,1,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,ierr)  
     
    
    if (myid==0) write(*,*) k, summ/ns 
   end do 
!      write(*,"(a,xi0,xi0,xf10.5)") "End loop", np,nn,1.0*np/nn 
      call MPI_Finalize(ierr)
    
   end program 
   
   
   

      SUBROUTINE init_random_seed(myid)
            INTEGER :: i, n, clock
            INTEGER, DIMENSION(:), ALLOCATABLE :: seed
            CALL RANDOM_SEED(size = n)
            ALLOCATE(seed(n))
            CALL SYSTEM_CLOCK(COUNT=clock)
            seed =clock+37*(/(i-1,i=1,n)/)+1891202*myid+1
            CALL RANDOM_SEED(PUT = seed)
            DEALLOCATE(seed)
       END SUBROUTINE
