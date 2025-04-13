!   时间积分
    subroutine ti(ua,va,za,ub,vb,zb,uc,vc,zc,rm,f,d,dt,zo,m,n)
    dimension ua(m,n),va(m,n),za(m,n),ub(m,n),vb(m,n),zb(m,n),    &
        uc(m,n),vc(m,n),zc(m,n),rm(m,n),f(m,n)
!   中央差格式
        c=0.25/d
        m1=m-1
        n1=n-1
        do 10 i=2,m1
        do 10 j=2,n1
        e=-c*rm(i,j)*((ub(i+1,j)+ub(i,j))*(ub(i+1,j)-ub(i,j))&
                      +(ub(i,j)+ub(i-1,j))*(ub(i,j)-ub(i-1,j))&
               +(vb(i,j-1)+vb(i,j))*(ub(i,j)-ub(i,j-1))&
               +(vb(i,j)+vb(i,j+1))*(ub(i,j+1)-ub(i,j))&
               +19.6*(zb(i+1,j)-zb(i-1,j)))+f(i,j)*vb(i,j)
        uc(i,j)=ua(i,j)+e*dt
        g=-c*rm(i,j)*((ub(i+1,j)+ub(i,j))*(vb(i+1,j)-vb(i,j))&
                     +(ub(I,j)+ub(i-1,j))*(vb(i,j)-vb(i-1,j))&
               +(vb(I,j-1)+vb(i,j))*(vb(i,j)-vb(i,j-1))&
               +(vb(I,j)+vb(i,j+1))*(vb(i,j+1)-vb(i,j))&
               +19.6*(zb(i,j+1)-zb(i,j-1)))-f(i,j)*ub(i,j)
    10  vc(i,j)=va(i,j)+g*dt
        do 20 i=2,m1
        do 20 j=2,n1
        h=-c*rm(i,j)*((ub(I+1,j)+ub(i,j))*(zb(i+1,j)-zb(i,j))&
               +(ub(I,j)+ub(i-1,j))*(zb(i,j)-zb(i-1,j))&
               +(vb(I,j-1)+vb(i,j))*(zb(i,j)-zb(i,j-1))&
               +(vb(I,j)+vb(i,j+1))*(zb(i,j+1)-zb(i,j))&
               +2.0*(zb(i,j)-zo)*(ub(i+1,j)-ub(i-1,j)+vb(i,j+1)-vb(i,j-1)))
    20  zc(i,j)=za(i,j)+h*dt
        return
        end
    