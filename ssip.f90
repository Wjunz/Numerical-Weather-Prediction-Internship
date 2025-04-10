!   五点平滑子程序
    subroutine ssip(a,w,s,m,n,c)
    dimension a(m,n),w(m,n)
    if(c==0) then !执行正平滑
        !正平滑
        do i=2,m-1
            do j=2,n-1
                w(i,j)=a(i,j)+0.25*s*(a(i+1,j)+a(i-1,j)+&
                                    a(i,j+1)+a(i,j-1)-4.0*a(i,j))
            end do
        end do
        do i=2,m-1
            do j=2,n-1
                a(i,j)=w(i,j)
            end do
        end do
    else if(c==1) then !执行正逆平滑
        !正平滑
        do i=2,m-1
            do j=2,n-1
                w(i,j)=a(i,j)+0.25*s*(a(i+1,j)+a(i-1,j)+&
                                    a(i,j+1)+a(i,j-1)-4.0*a(i,j))
            end do
        end do
        do i=2,m-1
            do j=2,n-1
                a(i,j)=w(i,j)
            end do
        end do
        !逆平滑
        do i=2,m-1
            do j=2,n-1
                w(i,j)=a(i,j)-0.25*s*(a(i+1,j)+a(i-1,j)+&
                                    a(i,j+1)+a(i,j-1)-4.0*a(i,j))
            end do
        end do
    end if
    do i=2,m-1
        do j=2,n-1
            a(i,j)=w(i,j)
        end do
    end do
    return  
    end
