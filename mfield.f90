program langevin
implicit none
real (kind=8) th,c,mu,lam,dxi,mu1,t1,t2,xi
integer (kind=8) i,j,k,n,num,l
th=.999
lam=.3
dxi=1.0/10000
n=ceiling(th/dxi)
dxi=th/n
num=1000
do l=1,320
c=(l)*1e-7
k=0
mu1=log(1/(1-th))
mu=0
call cpu_time(t1)
do while ((dabs(mu1-mu)>5e-7).and.(k<10000))
mu=mu1
mu1=0
do i=1,n
    xi=(dexp(-lam*(i-1)*dxi)-dexp(-lam*i*dxi))/(1-exp(-lam*th))
    j=max(min(floor(num*dlog(((1-(i-.5)*dxi)*(1-dexp(-mu/num))+c)/((1-th)*(1-dexp(-mu/num))+c))/mu),num-1),0)
!    do while ((dexp(j*mu/num)>((1-(i-.5)*dxi)/(1-th)-c*(dexp(mu*j/num)-1)/((1-th)*(1-dexp(-mu/num))))).and.(j.ge.0))
!        j=j-1
!    end do
!    *((1-xi)**(num-1-j))
    mu1=mu1+xi*((1-xi)**(num-1-j))*(dlog((1-(i-.5)*dxi)-c*(dexp(mu*j/num)-1)/(1-dexp(-mu/num)))-log(1-th))
end do
mu1=mu1
!write(*,*) mu1,dabs(mu1-mu)
k=k+1
end do
call cpu_time(t2)
write(*,*) t2-t1,k,c,num,th,mu1,dabs(mu1-mu)
end do
end program
