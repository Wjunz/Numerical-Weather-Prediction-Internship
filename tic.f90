subroutine tic(ua,va,za,ub,vb,zb,uc,vc,zc,rm,f,d,dt,zo,m,n)
    dimension ua(m,n),va(m,n),za(m,n),ub(m,n),vb(m,n),zb(m,n),    &
              uc(m,n),vc(m,n),zc(m,n),rm(m,n),f(m,n)
    c=0.25/d
    m1=m-1
    n1=n-1
    
    ! 保存初始状态到ua,va,za
    do i=1,m
        do j=1,n
            ua(i,j) = ub(i,j)
            va(i,j) = vb(i,j)
            za(i,j) = zb(i,j)
        end do
    end do

    !第一步：Δt/2向前差
    do i=2,m1
    do j=2,n1
        ! 计算u分量
        e = -c*rm(i,j)*((ua(i+1,j)+ua(i,j))*(ua(i+1,j)-ua(i,j)) &
                    + (ua(i,j)+ua(i-1,j))*(ua(i,j)-ua(i-1,j))    &
                    + (va(i,j-1)+va(i,j))*(ua(i,j)-ua(i,j-1))    &
                    + (va(i,j)+va(i,j+1))*(ua(i,j+1)-ua(i,j))    &
                    + 19.6*(za(i+1,j)-za(i-1,j))) + f(i,j)*va(i,j)
        ! 计算v分量
        g = -c*rm(i,j)*((ua(i+1,j)+ua(i,j))*(va(i+1,j)-va(i,j)) &
                    + (ua(i,j)+ua(i-1,j))*(va(i,j)-va(i-1,j))    &
                    + (va(i,j-1)+va(i,j))*(va(i,j)-va(i,j-1))    &
                    + (va(i,j)+va(i,j+1))*(va(i,j+1)-va(i,j))    &
                    + 19.6*(za(i,j+1)-za(i,j-1))) - f(i,j)*ua(i,j)
        ! 计算z分量
        h = -c*rm(i,j)*((ua(i+1,j)+ua(i,j))*(za(i+1,j)-za(i,j)) &
                    + (ua(i,j)+ua(i-1,j))*(za(i,j)-za(i-1,j))    &
                    + (va(i,j-1)+va(i,j))*(za(i,j)-za(i,j-1))    &
                    + (va(i,j)+va(i,j+1))*(za(i,j+1)-za(i,j))    &
                    + 2.0*(za(i,j)-zo)*(ua(i+1,j)-ua(i-1,j)+va(i,j+1)-va(i,j-1)))
        
        ! 更新中间量到ub数组
        ub(i,j) = ua(i,j) + e*(dt/2)
        vb(i,j) = va(i,j) + g*(dt/2)
        zb(i,j) = za(i,j) + h*(dt/2)
    end do
    end do

    !第二步：Δt/2中央差
    do i=2,m1
    do j=2,n1
        ! 计算u分量
        e = -c*rm(i,j)*((ub(i+1,j)+ub(i,j))*(ub(i+1,j)-ub(i,j)) &
                    + (ub(i,j)+ub(i-1,j))*(ub(i,j)-ub(i-1,j))     &
                    + (vb(i,j-1)+vb(i,j))*(ub(i,j)-ub(i,j-1))     &
                    + (vb(i,j)+vb(i,j+1))*(ub(i,j+1)-ub(i,j))     &
                    + 19.6*(zb(i+1,j)-zb(i-1,j))) + f(i,j)*vb(i,j)
        ! 计算v分量
        g = -c*rm(i,j)*((ub(i+1,j)+ub(i,j))*(vb(i+1,j)-vb(i,j)) &
                    + (ub(i,j)+ub(i-1,j))*(vb(i,j)-vb(i-1,j))     &
                    + (vb(i,j-1)+vb(i,j))*(vb(i,j)-vb(i,j-1))     &
                    + (vb(i,j)+vb(i,j+1))*(vb(i,j+1)-vb(i,j))     &
                    + 19.6*(zb(i,j+1)-zb(i,j-1))) - f(i,j)*ub(i,j)
        ! 计算z分量
        h = -c*rm(i,j)*((ub(i+1,j)+ub(i,j))*(zb(i+1,j)-zb(i,j)) &
                    + (ub(i,j)+ub(i-1,j))*(zb(i,j)-zb(i-1,j))     &
                    + (vb(i,j-1)+vb(i,j))*(zb(i,j)-zb(i,j-1))     &
                    + (vb(i,j)+vb(i,j+1))*(zb(i,j+1)-zb(i,j))     &
                    + 2.0*(zb(i,j)-zo)*(ub(i+1,j)-ub(i-1,j)+vb(i,j+1)-vb(i,j-1)))
        
        ! 更新中间量到uc数组
        uc(i,j) = ua(i,j) + e*dt/2  ! 使用初始状态作为基点
        vc(i,j) = va(i,j) + g*dt/2
        zc(i,j) = za(i,j) + h*dt/2
    end do
    end do

    !第三步：Δt中央差
    do i=2,m1
    do j=2,n1
        ! 计算u分量
        e = -c*rm(i,j)*((uc(i+1,j)+uc(i,j))*(uc(i+1,j)-uc(i,j)) &
                    + (uc(i,j)+uc(i-1,j))*(uc(i,j)-uc(i-1,j))     &
                    + (vc(i,j-1)+vc(i,j))*(uc(i,j)-uc(i,j-1))     &
                    + (vc(i,j)+vc(i,j+1))*(uc(i,j+1)-uc(i,j))     &
                    + 19.6*(zc(i+1,j)-zc(i-1,j))) + f(i,j)*vc(i,j)
        ! 计算v分量
        g = -c*rm(i,j)*((uc(i+1,j)+uc(i,j))*(vc(i+1,j)-vc(i,j)) &
                    + (uc(i,j)+uc(i-1,j))*(vc(i,j)-vc(i-1,j))     &
                    + (vc(i,j-1)+vc(i,j))*(vc(i,j)-vc(i,j-1))     &
                    + (vc(i,j)+vc(i,j+1))*(vc(i,j+1)-vc(i,j))     &
                    + 19.6*(zc(i,j+1)-zc(i,j-1))) - f(i,j)*uc(i,j)
        ! 计算z分量
        h = -c*rm(i,j)*((uc(i+1,j)+uc(i,j))*(zc(i+1,j)-zc(i,j)) &
                    + (uc(i,j)+uc(i-1,j))*(zc(i,j)-zc(i-1,j))     &
                    + (vc(i,j-1)+vc(i,j))*(zc(i,j)-zc(i,j-1))     &
                    + (vc(i,j)+vc(i,j+1))*(zc(i,j+1)-zc(i,j))     &
                    + 2.0*(zc(i,j)-zo)*(uc(i+1,j)-uc(i-1,j)+vc(i,j+1)-vc(i,j-1)))
        
        ! 最终结果
        uc(i,j) = ua(i,j) + e*dt  ! 基于初始状态的全步长更新
        vc(i,j) = va(i,j) + g*dt
        zc(i,j) = za(i,j) + h*dt
    end do
    end do

    return
end