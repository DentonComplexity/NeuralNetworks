program langevin
implicit none
real (kind=8) lam, m,th,t1,t2,x,mx,r,st
real (kind=8),dimension(:), allocatable :: y,te
real (kind=8),dimension(:,:),allocatable :: z
integer (kind=8) n,i,f,l,p,q,j,s
integer (kind=1), dimension(:), allocatable :: k
lam=4
th=.9999
st=0
!s=0
n=40
f=400000
!r=0
!write(*,*) r
allocate(y(n))
allocate(te(n))
allocate(k(n))
allocate(z(n,n))
call cpu_time(t1)
do s=1,100
!k=0
!call cpu_time(t1)
!z=0
do i=1,n
call random_number(x)
!s=s+1
y(i)=dlog(1+dmod(dlog(1-x)/lam,th))
    do j=1,n
!call random_number(x)
!write(*,*) n,i,j,mod(j+1,n),mod(j-1,n)
if (0==mod((j-1-i)*(j+1-i)*(j+2-i)*(j-2-i),n)) then
            z(i,j)=x*0*(1-th)+(s*1.1e-4)+2.25e-2
        else
            z(i,j)=0
        end if
!if (i==j) then
!z(i,j)=0.d0
!else
!z(i,j)=x*0*(1-th)+(s*2e-4)
!end if
    end do
end do
y=y-dlog(1-th)
!write(*,*) z
!call cpu_time(t2)
!write(*,*) t2-t1
i=0
q=0
!call cpu_time(t1)
do while (i<f)
    j=1
    k=0
    te=0
    mx=minval(y)
    !s=minval(minloc(y))
    p=0
    do while (j>0)
        j=0
        !s=1
        do l=1,n
            st=dexp(y(l)-mx)
!write(*,*) i,j,y
            !do s=1,n
            !    st=st+z(l,s)*k(s)
            !end do
            if ((((st-te(l)).le.1).or.(y(l)==mx)).and.(k(l).ne.1)) then
                call random_number(x)
                !s=s+1
                !if (j==0) then
                !    y(l)=y(l)+dlog(1+dmod(dlog(1-x)/lam,th))-dlog(1-th)
                !else
                    y(l)=mx+dlog(1+dmod(dlog(1-x)/lam,th))-dlog(1-th)
                !end if
                j=j+1
                i=i+1
                k(l)=1
                te=te+z(l,:)
            else
                if ((p>0).and.(k(l).ne.1)) then
                    y(l)=mx+dlog(st-te(l))
                end if
            end if
                !if ((s>2837987).and.(s<2837991)) then
                !   write(*,*) j,y,k,l,s,dlog(1+dmod(dlog(1-x)/lam,th))-dlog(1-th),te
                !end if
        end do
        !if (j>n) then
        !    write(*,*) "Too many firings"
        !    exit
        !end if
        if (p>n) then
            write(*,*) "Too many iterations"
            exit
        end if
        p=p+1
    end do
    q=q+1
end do
!call cpu_time(t2)
!write(*,*) t2-t1
write(*,*) z(2,1)*(1-th),minval(y)/(f+1),minval(y)
!write(*,*) i,q
end do
call cpu_time(t2)
write(*,*) t2-t1
end
