subroutine mysub(Y, Yout, alpha_1, alpha_2, alpha_3, alpha_4, h, h_angle, eps)
    implicit none

    real, parameter :: PI = 3.1415926
    double precision, intent(in) :: alpha_1, alpha_2, alpha_3, alpha_4, h, h_angle, eps
    double precision, dimension(8, 2, 2), intent(in) :: Y
    double precision, dimension(8, 2, 2), intent(out) :: Yout

    double precision, dimension(2, 2) :: Reg, Img, Regt, Imgt, Redg, Imdg, Redgt, Imdgt, s0, sx, sz, A_1, A_3
    double complex, dimension(2, 2) :: g, gt, dg, dgt, n, nt, rhs, rhst    

    Reg   = Y(1,:,:)
    Img   = Y(2,:,:)
    Regt  = Y(3,:,:)
    Imgt  = Y(4,:,:)
    Redg  = Y(5,:,:)
    Imdg  = Y(6,:,:)
    Redgt = Y(7,:,:)
    Imdgt = Y(8,:,:)
    s0    = reshape((/ 1, 0, 0, 1 /), shape(s0))
    sx    = reshape((/ 0, 1, 1, 0 /), shape(sx))
    sz    = reshape((/ 1, 0, 0,-1 /), shape(sz))

    A_1   = alpha_1*sx + alpha_2*sz ! 1st component of the spin-orbit vector potential
    A_3   = alpha_3*sx + alpha_4*sz ! 3rd component of the spin-orbit vector potential
    ! 2nd component is zero.

    g   = Reg   + (0d0,1d0)*Img
    gt  = Regt  + (0d0,1d0)*Imgt
    dg  = Redg  + (0d0,1d0)*Imdg
    dgt = Redgt + (0d0,1d0)*Imdgt

    ! These are the inverse matrices defined in the equation 3.26 in my Master's thesis

    n(1,1)  = ( 1 + g(2,1)*gt(1,2) + g(2,2)*gt(2,2) ) / ( 1 + g(1,1)*gt(1,1) + g(2,1)*gt(1,2) + g(1,2)*gt(2,1)&
    + g(1,2)*g(2,1)*gt(1,2)*gt(2,1) - g(1,1)*g(2,2)*gt(1,2)*gt(2,1) + g(2,2)*gt(2,2) - g(1,2)*g(2,1)*gt(1,1)*gt(2,2)&
    + g(1,1)*g(2,2)*gt(1,1)*gt(2,2) )

    n(1,2)  = -( g(1,1)*gt(1,2) + g(1,2)*gt(2,2) ) / ( 1 + g(1,1)*gt(1,1) + g(2,1)*gt(1,2) + g(1,2)*gt(2,1)&
    + g(1,2)*g(2,1)*gt(1,2)*gt(2,1) - g(1,1)*g(2,2)*gt(1,2)*gt(2,1) + g(2,2)*gt(2,2) - g(1,2)*g(2,1)*gt(1,1)*gt(2,2)&
    + g(1,1)*g(2,2)*gt(1,1)*gt(2,2) )

    n(2,1)  = -( g(2,1)*gt(1,1) + g(2,2)*gt(2,1) ) / ( 1 + g(1,1)*gt(1,1) + g(2,1)*gt(1,2) + g(1,2)*gt(2,1)&
    + g(1,2)*g(2,1)*gt(1,2)*gt(2,1) - g(1,1)*g(2,2)*gt(1,2)*gt(2,1) + g(2,2)*gt(2,2) - g(1,2)*g(2,1)*gt(1,1)*gt(2,2)&
    + g(1,1)*g(2,2)*gt(1,1)*gt(2,2) )

    n(2,2)  = ( 1 + g(1,1)*gt(1,1) + g(1,2)*gt(2,1) ) / ( 1 + g(1,1)*gt(1,1) + g(2,1)*gt(1,2) + g(1,2)*gt(2,1)&
    + g(1,2)*g(2,1)*gt(1,2)*gt(2,1) - g(1,1)*g(2,2)*gt(1,2)*gt(2,1) + g(2,2)*gt(2,2) - g(1,2)*g(2,1)*gt(1,1)*gt(2,2)&
    + g(1,1)*g(2,2)*gt(1,1)*gt(2,2) )
   
    nt(1,1) = ( 1 + g(1,2)*gt(2,1) + g(2,2)*gt(2,2) ) / ( 1 + g(1,1)*gt(1,1) + g(2,1)*gt(1,2) + g(1,2)*gt(2,1)&
    + g(1,2)*g(2,1)*gt(1,2)*gt(2,1) - g(1,1)*g(2,2)*gt(1,2)*gt(2,1) + g(2,2)*gt(2,2) - g(1,2)*g(2,1)*gt(1,1)*gt(2,2)&
    + g(1,1)*g(2,2)*gt(1,1)*gt(2,2) )

    nt(1,2) = -( g(1,2)*gt(1,1) + g(2,2)*gt(1,2) ) / ( 1 + g(1,1)*gt(1,1) + g(2,1)*gt(1,2) + g(1,2)*gt(2,1)&
    + g(1,2)*g(2,1)*gt(1,2)*gt(2,1) - g(1,1)*g(2,2)*gt(1,2)*gt(2,1) + g(2,2)*gt(2,2) - g(1,2)*g(2,1)*gt(1,1)*gt(2,2)&
    + g(1,1)*g(2,2)*gt(1,1)*gt(2,2) )

    nt(2,1) = -( g(1,1)*gt(2,1) + g(2,1)*gt(2,2) ) / ( 1 + g(1,1)*gt(1,1) + g(2,1)*gt(1,2) + g(1,2)*gt(2,1)&
    + g(1,2)*g(2,1)*gt(1,2)*gt(2,1) - g(1,1)*g(2,2)*gt(1,2)*gt(2,1) + g(2,2)*gt(2,2) - g(1,2)*g(2,1)*gt(1,1)*gt(2,2)&
    + g(1,1)*g(2,2)*gt(1,1)*gt(2,2) )

    nt(2,2) = ( 1 + g(1,1)*gt(1,1) + g(2,1)*gt(1,2) ) / ( 1 + g(1,1)*gt(1,1) + g(2,1)*gt(1,2) + g(1,2)*gt(2,1)&
    + g(1,2)*g(2,1)*gt(1,2)*gt(2,1) - g(1,1)*g(2,2)*gt(1,2)*gt(2,1) + g(2,2)*gt(2,2) - g(1,2)*g(2,1)*gt(1,1)*gt(2,2)&
    + g(1,1)*g(2,2)*gt(1,1)*gt(2,2) )


    ! This is the equation 4.2 in my Master's thesis. The 1st term on the LHS is not here as we want to give the solver
    ! the derivative terms.

    rhs = - (0d0,1d0)*2*eps*MATMUL(s0,g)& ! This is the energy term. 1st on the RHS.

    ! The 2nd term on the LHS coming from the normal gradient part of the covariant derivative.
    + 2*MATMUL(dg, MATMUL(gt, MATMUL(n, dg)))& 

    ! These two terms are the exchange field. With angle 0 you only get sigma_3 and with angle 90 you only get sigma_1.
    + (0d0,1d0)*h*COS(h_angle*PI/180)*(MATMUL(g, sz) - MATMUL(sz, g))&
    + (0d0,1d0)*h*SIN(h_angle*PI/180)*(MATMUL(g, sx) - MATMUL(sx, g))&

    ! These are the second order spin-orbit terms. Both A1 and A3 components from the vector potential contribute.
    + MATMUL(A_3, MATMUL(A_3, g)) - MATMUL(g, MATMUL(A_3, A_3)) + 2*MATMUL(g, MATMUL(A_3, MATMUL(nt, A_3)))&
    - 2*MATMUL(g, MATMUL(A_3, MATMUL(nt, MATMUL(gt, MATMUL(A_3, g))))) + 2*MATMUL(A_3, MATMUL(g, MATMUL(nt, A_3)))&
    - 2*MATMUL(A_3, MATMUL(g, MATMUL(nt, MATMUL(gt, MATMUL(A_3, g)))))&

    + MATMUL(A_1, MATMUL(A_1, g)) - MATMUL(g, MATMUL(A_1, A_1)) + 2*MATMUL(g, MATMUL(A_1, MATMUL(nt, A_1)))&
    - 2*MATMUL(g, MATMUL(A_1, MATMUL(nt, MATMUL(gt, MATMUL(A_1, g))))) + 2*MATMUL(A_1, MATMUL(g, MATMUL(nt, A_1)))&
    - 2*MATMUL(A_1, MATMUL(g, MATMUL(nt, MATMUL(gt, MATMUL(A_1, g)))))&

    ! This is the term where the A3 component couples to the gradient terms. Only A3 couples to the gradients
    ! since it is in the direction of the nanowire and gradients to other directions are assumed small.
    + (0d0,1d0)*2*(MATMUL(dg, MATMUL(nt, A_3)) - MATMUL(dg, MATMUL(nt, MATMUL(gt, MATMUL(A_3, g))))&
    + MATMUL(A_3, MATMUL(n, dg)) - MATMUL(g, MATMUL(A_3, MATMUL(gt, MATMUL(n, dg)))))


    ! This is the tilde counterpart to the equation 4.2. Change g <-> gt, n <-> nt and take complex conjugate
    ! of the scalars. The explanations of the terms are as above.

    rhst = -(0d0,1d0)*2*eps*MATMUL(s0,gt)&

    + 2*MATMUL(dgt, MATMUL(g, MATMUL(nt, dgt)))&

    - (0d0,1d0)*h*COS(h_angle*PI/180)*(MATMUL(gt, sz) - MATMUL(sz, gt))&
    - (0d0,1d0)*h*SIN(h_angle*PI/180)*(MATMUL(gt, sx) - MATMUL(sx, gt))&

    + MATMUL(A_3, MATMUL(A_3, gt)) - MATMUL(gt, MATMUL(A_3, A_3)) + 2*MATMUL(gt, MATMUL(A_3, MATMUL(n, A_3)))&
    - 2*MATMUL(gt, MATMUL(A_3, MATMUL(n, MATMUL(g, MATMUL(A_3, gt))))) + 2*MATMUL(A_3, MATMUL(gt, MATMUL(n, A_3)))&
    - 2*MATMUL(A_3, MATMUL(gt, MATMUL(n, MATMUL(g, MATMUL(A_3, gt)))))&

    + MATMUL(A_1, MATMUL(A_1, gt)) - MATMUL(gt, MATMUL(A_1, A_1)) + 2*MATMUL(gt, MATMUL(A_1, MATMUL(n, A_1)))&
    - 2*MATMUL(gt, MATMUL(A_1, MATMUL(n, MATMUL(g, MATMUL(A_1, gt))))) + 2*MATMUL(A_1, MATMUL(gt, MATMUL(n, A_1)))&
    - 2*MATMUL(A_1, MATMUL(gt, MATMUL(n, MATMUL(g, MATMUL(A_1, gt)))))&

    - (0d0,1d0)*2*(MATMUL(dgt, MATMUL(n, A_3)) - MATMUL(dgt, MATMUL(n, MATMUL(g, MATMUL(A_3, gt))))&
    + MATMUL(A_3, MATMUL(nt, dgt)) - MATMUL(gt, MATMUL(A_3, MATMUL(g, MATMUL(nt, dgt)))))

    Yout(1,:,:) = Redg
    Yout(2,:,:) = Imdg
    Yout(3,:,:) = Redgt
    Yout(4,:,:) = Imdgt
    Yout(5,:,:) = REALPART(rhs)
    Yout(6,:,:) = IMAGPART(rhs)
    Yout(7,:,:) = REALPART(rhst)
    Yout(8,:,:) = IMAGPART(rhst)

end subroutine mysub
