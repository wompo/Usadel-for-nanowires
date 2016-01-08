subroutine gsub(Yout, delta, eps, offset, phase)

    implicit none
    double precision, intent(in) :: delta, eps, offset, phase
    double precision, dimension(4, 2, 2), intent(out) :: Yout
    double precision, dimension(2, 2) :: Reg, Img, Regt, Imgt, sy

    sy  = reshape((/ 0, 1, -1, 0 /), shape(sy))

    Reg     = REALPART(delta * EXP((0d0,1d0) * phase) / (eps + (0d0,1d0) * SQRT(delta**2 - (eps + (0d0,1d0) * offset)**2)) * sy)
    Img     = IMAGPART(delta * EXP((0d0,1d0) * phase) / (eps + (0d0,1d0) * SQRT(delta**2 - (eps + (0d0,1d0) * offset)**2)) * sy)
    Regt    = REALPART(delta * EXP(-(0d0,1d0) * phase) / (eps + (0d0,1d0) * SQRT(delta**2 - (eps + (0d0,1d0) * offset)**2)) * sy)
    Imgt    = IMAGPART(delta * EXP(-(0d0,1d0) * phase) / (eps + (0d0,1d0) * SQRT(delta**2 - (eps + (0d0,1d0) * offset)**2)) * sy)

    Yout(1,:,:) = Reg
    Yout(2,:,:) = Img
    Yout(3,:,:) = Regt
    Yout(4,:,:) = Imgt

end subroutine gsub
