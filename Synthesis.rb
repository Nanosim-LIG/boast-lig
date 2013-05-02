require './BOAST.rb'
require 'rubygems'
require 'narray'
module ConvolutionGenerator

  def ConvolutionGenerator::synthesis_free_ref
    lang = ConvolutionGenerator::get_lang
    ConvolutionGenerator::set_lang(ConvolutionGenerator::FORTRAN)
    kernel = CKernel::new
    kernel.lang = ConvolutionGenerator::FORTRAN
    function_name = "synthesis_free_ref"
    n = Variable::new("n",Int,{:direction => :in, :signed => false})
    ndat = Variable::new("ndat",Int,{:direction => :in, :signed => false})
    dim_in_min = 0
    dim_in_max = n*2-1
    dim_out_min = -7
    dim_out_max = n*2+6
#    dim_in_min = 0
#    dim_in_max = n*2-1
#    dim_out_min = -7
#    dim_out_max = n*2+6
    x = Variable::new("x",Real,{:direction => :in, :dimension => [ Dimension::new(dim_in_min, dim_in_max), Dimension::new(ndat) ] })
    y = Variable::new("y",Real,{:direction => :out, :dimension => [ Dimension::new(ndat), Dimension::new(dim_out_min, dim_out_max) ] })
    p = Procedure::new(function_name, [n,ndat,x,y])
    kernel.code.print <<EOF
subroutine synthesis_free_ref(n,ndat,x,y)
  implicit real(kind=8) (a-h,o-z)
  dimension x(0:2*n-1,ndat),y(ndat,-7:2*n+6)
  real(kind=8) ch(-8:9) ,cg(-8:9)
  !       Daubechy S16
  data ch  /  0.d0 , -0.0033824159510050025955D0, & 
       -0.00054213233180001068935D0, 0.031695087811525991431D0, & 
       0.0076074873249766081919D0, -0.14329423835127266284D0, & 
       -0.061273359067811077843D0, 0.48135965125905339159D0,  & 
       0.77718575169962802862D0,0.36444189483617893676D0, &
       -0.051945838107881800736D0,-0.027219029917103486322D0, &
       0.049137179673730286787D0,0.0038087520138944894631D0, &
       -0.014952258337062199118D0,-0.00030292051472413308126D0, &
       0.0018899503327676891843D0 , 0.d0 /
  data cg  / 0.d0 , -0.0018899503327676891843D0, &
       -0.00030292051472413308126D0, 0.014952258337062199118D0, &
       0.0038087520138944894631D0, -0.049137179673730286787D0, &
       -0.027219029917103486322D0, 0.051945838107881800736D0, &
       0.36444189483617893676D0, -0.77718575169962802862D0, &
       0.48135965125905339159D0, 0.061273359067811077843D0, &
       -0.14329423835127266284D0, -0.0076074873249766081919D0, &
       0.031695087811525991431D0, 0.00054213233180001068935D0, &
       -0.0033824159510050025955D0 , 0.d0 /

  do j=1,ndat

     i=-4
     so=0.d0
     do l=max(i-n+1,-4),min(i,4)
        so=so+ch(2*l+1)*x(i-l,j)+cg(2*l+1)*x(n+i-l,j)
     enddo
     y(j,2*i+1)=so

     do i=-3,n+2
        se=0.d0
        so=0.d0
        do l=max(i-n+1,-4),min(i,4)
           se=se+ch(2*l  )*x(i-l,j)+cg(2*l  )*x(n+i-l,j)
           so=so+ch(2*l+1)*x(i-l,j)+cg(2*l+1)*x(n+i-l,j)
        enddo
        y(j,2*i  )=se
        y(j,2*i+1)=so
     enddo

     i=n+3
     se=0.d0
     do l=max(i-n+1,-4),min(i,4)
        se=se+ch(2*l  )*x(i-l,j)+cg(2*l  )*x(n+i-l,j)
     enddo
     y(j,2*i  )=se

  enddo

  return
END SUBROUTINE synthesis_free_ref
EOF
    kernel.procedure = p
    ConvolutionGenerator::set_lang(lang)
    return kernel
  end
  def ConvolutionGenerator::synthesis_per_ref
    lang = ConvolutionGenerator::get_lang
    ConvolutionGenerator::set_lang(ConvolutionGenerator::FORTRAN)
    kernel = CKernel::new
    kernel.lang = ConvolutionGenerator::FORTRAN
    function_name = "synthesis_per_ref"
    n = Variable::new("n",Int,{:direction => :in, :signed => false})
    ndat = Variable::new("ndat",Int,{:direction => :in, :signed => false})
    dim_in_min = 0
    dim_in_max = n*2-1
    dim_out_min = 0
    dim_out_max = n*2-1
    x = Variable::new("x",Real,{:direction => :in, :dimension => [ Dimension::new(dim_in_min, dim_in_max), Dimension::new(ndat) ] })
    y = Variable::new("y",Real,{:direction => :out, :dimension => [ Dimension::new(ndat), Dimension::new(dim_out_min, dim_out_max) ] })
    p = Procedure::new(function_name, [n,ndat,x,y])
    kernel.code.print <<EOF
subroutine synthesis_per_ref(n,ndat,x,y)
  implicit none
  integer, parameter :: wp=kind(1.0d0)
  integer, intent(in) :: n,ndat
  real(wp), dimension(0:2*n-1,ndat), intent(in) :: x
  real(wp), dimension(ndat,0:2*n-1), intent(out) :: y
  !local variables
  integer :: i,j,k,l
  real(wp) :: so,se
  real(wp), dimension(-8:9) :: ch,cg
  !       Daubechy S16
  data ch  /  0.e0_wp , -0.0033824159510050025955_wp, & 
       -0.00054213233180001068935_wp, 0.031695087811525991431_wp, & 
       0.0076074873249766081919_wp, -0.14329423835127266284_wp, & 
       -0.061273359067811077843_wp, 0.48135965125905339159_wp,  & 
       0.77718575169962802862_wp,0.36444189483617893676_wp, &
       -0.051945838107881800736_wp,-0.027219029917103486322_wp, &
       0.049137179673730286787_wp,0.0038087520138944894631_wp, &
       -0.014952258337062199118_wp,-0.00030292051472413308126_wp, &
       0.0018899503327676891843_wp , 0.e0_wp /
  data cg  / 0.e0_wp , -0.0018899503327676891843_wp, &
       -0.00030292051472413308126_wp, 0.014952258337062199118_wp, &
       0.0038087520138944894631_wp, -0.049137179673730286787_wp, &
       -0.027219029917103486322_wp, 0.051945838107881800736_wp, &
       0.36444189483617893676_wp, -0.77718575169962802862_wp, &
       0.48135965125905339159_wp, 0.061273359067811077843_wp, &
       -0.14329423835127266284_wp, -0.0076074873249766081919_wp, &
       0.031695087811525991431_wp, 0.00054213233180001068935_wp, &
       -0.0033824159510050025955_wp , 0.e0_wp /

  do j=1,ndat

     do i=0,n-1
        se=0.e0_wp
        so=0.e0_wp
        do l=-4,4
           k=modulo(i-l,n)
           se=se+ch(2*l  )*x(  k,j)+cg(2*l  )*x(n+k  ,j)
           so=so+ch(2*l+1)*x(  k,j)+cg(2*l+1)*x(n+k  ,j)
        enddo
        y(j,2*i  )=se
        y(j,2*i+1)=so
     enddo

  enddo

END SUBROUTINE synthesis_per_ref
EOF
    kernel.procedure = p
    ConvolutionGenerator::set_lang(lang)
    return kernel
  end

  def ConvolutionGenerator::Synthesis(filt, center, unroll, free=false )
    kernel = CKernel::new
    ConvolutionGenerator::set_output( kernel.code )
    kernel.lang = ConvolutionGenerator::get_lang
    function_name = "synthesis"
    if free then
      function_name += "_free"
    else 
      function_name += "_per"
    end
    if unroll>0 then
      function_name += "_u#{unroll}"
    end

    n = Variable::new("n",Int,{:direction => :in, :signed => false})
    ndat = Variable::new("ndat",Int,{:direction => :in, :signed => false})
    lowfil = Variable::new("lowfil",Int,{:constant => (1-center) / 2})
    upfil = Variable::new("upfil",Int,{:constant => (filt.length - center) / 2})

    dim_in_min = 0
    dim_in_max = n*2-1
    if free then
      dim_out_min = lowfil - upfil #-7 in free BC
      dim_out_max = n*2 -1 - lowfil + upfil #2*n+6 in free BC
    else
      dim_out_min = 0
      dim_out_max = n*2-1
    end

    if free then
      lowlimit=lowfil-1 #(dim_out_min-1)/2 #0 in periodic, -4 in free BC
      uplimit=n-2+upfil #(dim_out_max-1)/2 #n-1 in periodic, n+2 in free BC
    else
      lowlimit=0
      uplimit=n-2
    end

    x = Variable::new("x",Real,{:direction => :in, :dimension => [ Dimension::new(dim_in_min, dim_in_max), Dimension::new(ndat) ] })
    y = Variable::new("y",Real,{:direction => :out, :dimension => [ Dimension::new(ndat), Dimension::new(dim_out_min, dim_out_max) ] })
    i = Variable::new("i",Int)
    j = Variable::new("j",Int)
    k = Variable::new("k",Int)
    l = Variable::new("l",Int)
    so = [Variable::new("so1",Real)]
    2.upto(unroll) { |index|
      so.push(Variable::new("so#{index}",Real))
    }
    se = [Variable::new("se1",Real)]
    2.upto(unroll) { |index|
      se.push(Variable::new("se#{index}",Real))
    }
    arr = ConstArray::new(filt,Real)

    fil = Variable::new("fil",Real,{:constant => arr,:dimension => [ Dimension::new((1-center),(filt.length - center)) ]})
    

    if ConvolutionGenerator::get_lang == C then
      $output.print "inline #{Int::new.decl} modulo( #{Int::new.decl} a, #{Int::new.decl} b) { return (a+b)%b;}\n"
      $output.print "inline #{Int::new.decl} min( #{Int::new.decl} a, #{Int::new.decl} b) { return a < b ? a : b;}\n"
      $output.print "inline #{Int::new.decl} max( #{Int::new.decl} a, #{Int::new.decl} b) { return a > b ? a : b;}\n"
    end

    p = Procedure::new(function_name, [n,ndat,x,y], [lowfil,upfil]) {
      
      i.decl
      j.decl
      k.decl
      l.decl
      fil.decl
      so.each{ |s| s.decl }
      se.each{ |s| s.decl }
      if ConvolutionGenerator::get_lang == ConvolutionGenerator::FORTRAN
        $output.print("!$omp parallel default(shared) shared(x,y,ndat,n)&\n")
        $output.print("!$omp private(i,j,k,l)&\n")
        $output.print("!$omp private(")
        $output.print(se.join(","))
        $output.print(")&\n")
        $output.print("!$omp private(")
        $output.print(so.join(","))
        $output.print(")\n")
      end
      if ConvolutionGenerator::get_lang == ConvolutionGenerator::C then
        $output.print("#pragma omp parallel private(")
        $output.print(se.join(","))
        $output.print(",")
        $output.print(so.join(","))
        $output.print(",i,j,k,l) shared(x,y)\n") 
      end
      $output.print("{\n") if ConvolutionGenerator::get_lang == ConvolutionGenerator::C
      $output.print("!$omp do\n") if ConvolutionGenerator::get_lang == ConvolutionGenerator::FORTRAN
      $output.print("#pragma omp for\n") if ConvolutionGenerator::get_lang == ConvolutionGenerator::C

      #internal loop taking care of BCs
      if free then
        forBC = For::new(l, FuncCall::new("max", -i,lowfil), FuncCall::new("min", upfil, n-1-i) ){ |e,o|
          (k === i+l).print
          e.each_index{ |ind|
            (e[ind] === e[ind] + fil[l*2]*x[k,j+ind] + fil[l*-2+3]*x[n+k,j+ind]).print
            (o[ind] === o[ind] + fil[l*2+1]*x[k,j+ind] - fil[l*-2+2]*x[n+k,j+ind]).print
          }
        }
      else
        forBC = For::new(l,lowfil,upfil){ |e,o|
          (k === FuncCall::new( "modulo", i+(l), n)).print
          e.each_index{ |ind|
            (e[ind] === e[ind] + fil[l*2]*x[k,j+ind] + fil[l*-2+3]*x[n+k,j+ind]).print
            (o[ind] === o[ind] + fil[l*2+1]*x[k,j+ind] - fil[l*-2+2]*x[n+k,j+ind]).print
          }
        }
      end

      for_noBC = For::new(l,lowfil,upfil){ |e,o|
        (k === i + l).print
        e.each_index{ |ind|
          (e[ind] === e[ind] + fil[l*2]*x[k,j+ind] + fil[l*-2+3]*x[n+k,j+ind]).print
          (o[ind] === o[ind] + fil[l*2+1]*x[k,j+ind] - fil[l*-2+2]*x[n+k,j+ind]).print
        }
      }
      
      inner_block_uniq = lambda { |e, o, inner_for|
        o.each{ |s| (s === 0.0).print }
        e.each{ |s| (s === 0.0).print }
        inner_for.unroll(e, o)
        o.each_index { |ind|
          (y[j+ind,n*2-1] === o[ind]).print
        }
        e.each_index { |ind|
          (y[j+ind,0] === e[ind]).print
        }
      }

      inner_block = lambda { |e, o, inner_for|
        o.each{ |s| (s === 0.0).print }
        e.each{ |s| (s === 0.0).print }
        inner_for.unroll(e, o)
        o.each_index { |ind|
          (y[j+ind,i*2+1] === o[ind]).print
        }
        e.each_index { |ind|
          (y[j+ind,i*2+2] === e[ind]).print
        }
      }


      if unroll > 0 then
        for1 = For::new(j,1,ndat-(unroll-1), unroll)
      else
        for1 = For::new(j,1,ndat)
      end
      for1.print
      if not free then
        (i === n-1).print
        inner_block_uniq.call(se, so, forBC)
      end

      For::new(i,lowlimit,-lowfil-1) {
        inner_block.call(se, so, forBC)
      }.print
      
      For::new(i,-lowfil,n-1-upfil) {
        inner_block.call(se, so, for_noBC)
      }.print

      For::new(i,n-upfil,uplimit) {
        inner_block.call(se, so, forBC)
      }.print
      for1.close
      $output.print("!$omp end do\n") if ConvolutionGenerator::get_lang == ConvolutionGenerator::FORTRAN


      if unroll>1 then
        $output.print("#pragma omp do\n") if ConvolutionGenerator::get_lang == ConvolutionGenerator::C
        $output.print("!$omp do\n") if ConvolutionGenerator::get_lang == ConvolutionGenerator::FORTRAN
        For::new(j,ndat-FuncCall::new("modulo",ndat,unroll)+1,ndat) {
          ind=0
          if not free then
            (i === n-1).print
            inner_block_uniq.call([se[ind]], [so[ind]], forBC)
          end
          
          For::new(i,lowlimit,uplimit) {
            inner_block.call([se[ind]], [so[ind]], forBC)
          }.print
          
        }.print
        $output.print("!$omp end do\n") if ConvolutionGenerator::get_lang == ConvolutionGenerator::FORTRAN
      end
      $output.print("!$omp end parallel\n") if ConvolutionGenerator::get_lang == ConvolutionGenerator::FORTRAN
      $output.print("}\n") if ConvolutionGenerator::get_lang == ConvolutionGenerator::C
    }
    p.print
    kernel.procedure = p
    return kernel
  end
end


