require './BOAST.rb'
require 'rubygems'
require 'narray'
module ConvolutionGenerator

  def ConvolutionGenerator::kinetic_per_ref
    lang = ConvolutionGenerator::get_lang
    ConvolutionGenerator::set_lang(ConvolutionGenerator::FORTRAN)
    kernel = CKernel::new
    kernel.lang = ConvolutionGenerator::FORTRAN
    function_name = "kinetic_per_ref"
    n1 = Variable::new("n1",Int,{:direction => :in, :signed => false})
    n2 = Variable::new("n2",Int,{:direction => :in, :signed => false})
    n3 = Variable::new("n3",Int,{:direction => :in, :signed => false})
    hgrid = Variable::new("hgrid",Real,{:direction => :in, :dimension => [3]})
    c = Variable::new("c",Real,{:direction => :in})
    x = Variable::new("x",Real,{:direction => :in, :dimension => [ Dimension::new(0, n1-1),  Dimension::new(0, n2-1), Dimension::new(0, n3-1)] })
    y = Variable::new("y",Real,{:direction => :out, :dimension => [ Dimension::new(0, n1-1),  Dimension::new(0, n2-1), Dimension::new(0, n3-1)] })
    
    p = Procedure::new(function_name, [n1,n2,n3,hgrid,x,y,c])
    kernel.code.print <<EOF
subroutine convolut_kinetic_per_c(n1,n2,n3,hgrid,x,y,c)
!   applies the kinetic energy operator onto x to get y. Works for periodic BC
  implicit none
  integer, parameter :: wp=kind(1.0d0), gp=kind(1.0d0)
  integer, intent(in) :: n1,n2,n3
  real(gp),intent(in)::c
  real(gp), dimension(3), intent(in) :: hgrid
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: x
  real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: y
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  integer :: i1,i2,i3,i,l,j
  real(wp) :: tt
  real(wp), dimension(3) :: scale
  real(wp), dimension(lowfil:lupfil,3) :: fil
  
  scale(:)=real(-.5_gp/hgrid(:)**2,wp)

  ! second derivative filters for Daubechies 16
  fil(0,:)=   -3.5536922899131901941296809374e0_wp*scale(:)
  fil(1,:)=    2.2191465938911163898794546405e0_wp*scale(:)
  fil(2,:)=   -0.6156141465570069496314853949e0_wp*scale(:)
  fil(3,:)=    0.2371780582153805636239247476e0_wp*scale(:)
  fil(4,:)=   -0.0822663999742123340987663521e0_wp*scale(:)
  fil(5,:)=    0.02207029188482255523789911295638968409e0_wp*scale(:)
  fil(6,:)=   -0.409765689342633823899327051188315485e-2_wp*scale(:)
  fil(7,:)=    0.45167920287502235349480037639758496e-3_wp*scale(:)
  fil(8,:)=   -0.2398228524507599670405555359023135e-4_wp*scale(:)
  fil(9,:)=    2.0904234952920365957922889447361e-6_wp*scale(:)
  fil(10,:)=  -3.7230763047369275848791496973044e-7_wp*scale(:)
  fil(11,:)=  -1.05857055496741470373494132287e-8_wp*scale(:)
  fil(12,:)=  -5.813879830282540547959250667e-11_wp*scale(:)
  fil(13,:)=   2.70800493626319438269856689037647576e-13_wp*scale(:)
  fil(14,:)=  -6.924474940639200152025730585882e-18_wp*scale(:)

  do i=1,14
     fil(-i,:)=fil(i,:)
  enddo
  
  do i3=0,n3
     ! (1/2) d^2/dx^2
     do i2=0,n2
        do i1=0,n1
           tt=x(i1,i2,i3)*c
           do l=lowfil,lupfil
              j=modulo(i1+l,n1+1)
              tt=tt+x(j   ,i2,i3)*fil(l,1)
           enddo
           y(i1,i2,i3)=tt
        enddo
     enddo
     
     ! + (1/2) d^2/dy^2
     do i1=0,n1
        do i2=0,n2
           tt=0.e0_wp
           do l=lowfil,lupfil
              j=modulo(i2+l,n2+1)
              tt=tt+x(i1,j   ,i3)*fil(l,2)
           enddo
           y(i1,i2,i3)=y(i1,i2,i3)+tt
        enddo
     enddo
     
  enddo

  ! + (1/2) d^2/dz^2
  do i2=0,n2
     do i1=0,n1
        do i3=0,n3
           tt=0.e0_wp
           do l=lowfil,lupfil
              j=modulo(i3+l,n3+1)
              tt=tt+x(i1,i2,   j)*fil(l,3)
           enddo
           y(i1,i2,i3)=y(i1,i2,i3)+tt
        enddo
     enddo
  enddo
  
END SUBROUTINE convolut_kinetic_per_c
EOF
    kernel.procedure = p
    ConvolutionGenerator::set_lang(lang)
    return kernel
  end

  def ConvolutionGenerator::kinetic(filt, center, unroll, free=[false,false,false] )
    kernel = CKernel::new
    ConvolutionGenerator::set_output( kernel.code )
    kernel.lang = ConvolutionGenerator::get_lang
    function_name = "kinetic"
    fr = free.collect { |e|
      if e
      then "f"
      else "p"
      end
    }
    function_name += "_" + fr.join("")
    if unroll>0 then
      function_name += "_u#{unroll}"
    end

    n1 = Variable::new("n1",Int,{:direction => :in, :signed => false})
    n2 = Variable::new("n2",Int,{:direction => :in, :signed => false})
    n3 = Variable::new("n3",Int,{:direction => :in, :signed => false})
    hgrid = Variable::new("hgrid",Real,{:direction => :in, :dimension => [Dimension::new(0,2)]})
    c = Variable::new("c",Real,{:direction => :in})
    x = Variable::new("x",Real,{:direction => :in, :dimension => [ Dimension::new(0, n1-1),  Dimension::new(0, n2-1), Dimension::new(0, n3-1)] })
    y = Variable::new("y",Real,{:direction => :out, :dimension => [ Dimension::new(0, n1-1),  Dimension::new(0, n2-1), Dimension::new(0, n3-1)] })
    if ConvolutionGenerator::get_lang == C then
      $output.print "inline #{Int::new.decl} modulo( #{Int::new.decl} a, #{Int::new.decl} b) { return (a+b)%b;}\n"
      $output.print "inline #{Int::new.decl} min( #{Int::new.decl} a, #{Int::new.decl} b) { return a < b ? a : b;}\n"
      $output.print "inline #{Int::new.decl} max( #{Int::new.decl} a, #{Int::new.decl} b) { return a > b ? a : b;}\n"
    end

    lowfil = Variable::new("lowfil",Int,{:constant => 1-center})
    upfil = Variable::new("upfil",Int,{:constant => filt.length - center})
    arr = ConstArray::new(filt,Real)
    zerocinq = Variable::new("zerocinq",Real,{:constant => "0.5"})
    fil = Variable::new("fil",Real,{:constant => arr,:dimension => [ Dimension::new((1-center),(filt.length - center)) ]})
    scal = Variable::new("scal",Real,{:dimension => [Dimension::new(0,2)]})
    i1 = Variable::new("i1",Int)
    i2 = Variable::new("i2",Int)
    i3 = Variable::new("i3",Int)
    l = Variable::new("l",Int)
    k = Variable::new("k",Int)
    tt = [Variable::new("tt1",Real)]
    2.upto(unroll) { |index|
      tt.push(Variable::new("tt#{index}",Real))
    }
    dims = [n1,n2,n3]
    iters = [i1,i2,i3]


    for_BC = lambda { |t,dim,unro,init|
      t.each { |ttt|
       (ttt === 0.0).print
      }
      if( free[dim[2]] ) then
        For::new( l, FuncCall::new( "max", -iters[dim[2]], lowfil), FuncCall::new( "min", upfil, dims[dim[2]] - 1 - iters[dim[2]] ) ) {
          (k === iters[dim[2]]+l).print
          t.each_index { |ind|
            j1 = (dim[2] == 0 ? k : (unro == 0 ? i1 + ind : i1))
            j2 = (dim[2] == 1 ? k : (unro == 1 ? i2 + ind : i2))
            j3 = (dim[2] == 2 ? k : (unro == 2 ? i3 + ind : i3))
            (t[ind] === t[ind] + x[j1,j2,j3]*fil[l] ).print
          }
        }.print
      else
        For::new( l, lowfil, upfil) {
          (k === FuncCall::new( "modulo", iters[dim[2]]+(l), dims[dim[2]] )).print
          t.each_index { |ind|
            j1 = (dim[2] == 0 ? k : (unro == 0 ? i1 + ind : i1))
            j2 = (dim[2] == 1 ? k : (unro == 1 ? i2 + ind : i2))
            j3 = (dim[2] == 2 ? k : (unro == 2 ? i3 + ind : i3))
            (t[ind] === t[ind] + x[j1,j2,j3]*fil[l] ).print
          }
        }.print
      end
      t.each_index { |ind|
          j1 = (unro == 0 ? i1 + ind : i1)
          j2 = (unro == 1 ? i2 + ind : i2)
          j3 = (unro == 2 ? i3 + ind : i3)
          (t[ind] === t[ind] * scal[dim[2]]).print
          (y[j1,j2,j3] === ( init ? t[ind] : y[j1,j2,j3] + t[ind] )).print
      }
    }
    for_noBC = lambda { |t,dim,unro,init|
      t.each { |ttt|
        (ttt === 0.0).print
      }
      For::new( l, lowfil, upfil) {
        t.each_index { |ind|
          j1 = (dim[2] == 0 ? i1 + l : (unro == 0 ? i1 + ind : i1))
          j2 = (dim[2] == 1 ? i2 + l : (unro == 1 ? i2 + ind : i2))
          j3 = (dim[2] == 2 ? i3 + l : (unro == 2 ? i3 + ind : i3))
          (t[ind] === t[ind] + x[j1,j2,j3]*fil[l] ).print
        }
      }.print
      t.each_index { |ind|
          j1 = (unro == 0 ? i1 + ind : i1)
          j2 = (unro == 1 ? i2 + ind : i2)
          j3 = (unro == 2 ? i3 + ind : i3)
          (t[ind] === t[ind] * scal[dim[2]]).print
          (y[j1,j2,j3] === ( init ? t[ind] : y[j1,j2,j3] + t[ind] )).print
      }
    }


    kinetic_d = lambda { |t,dim,unro,init|
      For::new(iters[dim[0]], 0, dims[dim[0]]-1, unro == dim[0] ? unroll : 1 ) {
        For::new(iters[dim[1]], 0, dims[dim[1]]-1, unro == dim[1] ? unroll : 1) {
          For::new(iters[dim[2]], 0, -lowfil) {
            for_BC.call(t, dim, unro, init)
          }.print
          For::new(iters[dim[2]], -lowfil+1, dims[dim[2]] - 1 - upfil) {
            for_noBC.call(t, dim, unro, init)
          }.print
          For::new(iters[dim[2]], dims[dim[2]] - upfil, dims[dim[2]]) {
            for_BC.call(t, dim, unro, init)
          }.print
        }.print
      }.print
    }

    p = Procedure::new(function_name, [n1,n2,n3,hgrid,x,y,c],[lowfil,upfil,zerocinq]){
      i1.decl
      i2.decl
      i3.decl
      l.decl
      k.decl
      scal.decl
      fil.decl
      tt.each { |t| t.decl }
      (0..2).each { |ind|
        (scal[ind] === -zerocinq / (hgrid[ind] * hgrid[ind])).print
      }
      kinetic_d.call(tt, [2,1,0], 1, true)
      kinetic_d.call(tt, [2,0,1], 0, false)
      kinetic_d.call(tt, [1,0,2], 0, false) 
    }
    p.print
    kernel.procedure = p
    return kernel
  end

end
