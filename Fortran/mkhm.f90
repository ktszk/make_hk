      !
      !> calculate hamiltonian_k and eigenvalue
      !
      subroutine mkhm(no,nr,nqx,nqy,nqz,hop,rvec,klist,en,uni,eq)
        implicit none
        integer(4),intent(in):: no,nr,nqx,nqy,nqz,klist(3)
        real(8),intent(in) rvec(nr,3)
        complex(8),intent(in) hop(nr,no,no)
        real(8),intent(out):: eq(no)
        complex(8),intent(out):: uni(no,no),en(no,no)

        integer(4) k,l,m
        real(8),Parameter::pi=3.141592653589793238462643383279d0
        real(8) eq1(no),phase
        complex(8) en2(no,no)

        ! Lapack blas
        integer(4) info
        real(8) rwork(3*no-2),tmp
        complex(8) work(2*no-1)

        !$omp parallel
        !$omp workshare
        en=0.0d0
        uni=0.0d0
        eq=0.0d0
        !$omp end workshare
        do m=1,no
           do l=m,no
              !$omp do private(phase),reduction(+: en)
              do k=1,nr
                 phase=2*pi*klist[1]/nqx*rvec(k,1)+2*pi*klist[2]/nqy*rvec(k,2)&
                      +2*pi*klist[3]/nqz*rvec(k,3) !k@
                 en(l,m)=en(l,m)+hop(k,l,m)*cmplx(cos(phase),-sin(phase))
              end do
              !$omp end do
              en(m,l)=conjg(en(l,m))
           end do
        end do
        !$omp end parallel

        en2=en
        call zheev('V','U',no,en2,no,eq1,work,2*no-1,rwork,info)
        uni=en2
        eq=eq1

        return
      end subroutine mkhm
