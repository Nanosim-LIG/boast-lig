require './BOAST.rb'
require 'rubygems'
require 'narray'
module ConvolutionGenerator

  def ConvolutionGenerator::analysis_free_ref
    lang = ConvolutionGenerator::get_lang
    ConvolutionGenerator::set_lang(ConvolutionGenerator::FORTRAN)
    kernel = CKernel::new
    kernel.lang = ConvolutionGenerator::FORTRAN
    function_name = "analysis_free_ref"
    n = Variable::new("n",Int,{:direction => :in, :signed => false})
    ndat = Variable::new("ndat",Int,{:direction => :in, :signed => false})
    dim_in_min = -7
    dim_in_max = n*2+6
    dim_out_min = 0
    dim_out_max = n*2-1
    x = Variable::new("x",Real,{:direction => :in, :dimension => [ Dimension::new(dim_in_min, dim_in_max), Dimension::new(ndat) ] })
    y = Variable::new("y",Real,{:direction => :out, :dimension => [ Dimension::new(ndat), Dimension::new(dim_out_min, dim_out_max) ] })
    p = Procedure::new(function_name, [n,ndat,x,y])
    kernel.code.print <<EOF
subroutine analysis_free_ref(n,ndat,x,y)
  implicit none
  integer, intent(in) :: n,ndat
  real(kind=8), dimension(-7:2*n+6,ndat), intent(in) :: x
  real(kind=8), dimension(ndat,0:2*n-1), intent(out) :: y
  !local variables
  integer :: i,j,l
  real(kind=8) :: ci,di
  real(kind=8), dimension(-7:8) :: ch,cg
  !       Daubechy S16
  data ch  /  -0.0033824159510050025955d0, & 
       -0.00054213233180001068935d0, 0.031695087811525991431d0, & 
       0.0076074873249766081919d0, -0.14329423835127266284d0, & 
       -0.061273359067811077843d0, 0.48135965125905339159d0,  & 
       0.77718575169962802862d0,0.36444189483617893676d0, &
       -0.051945838107881800736d0,-0.027219029917103486322d0, &
       0.049137179673730286787d0,0.0038087520138944894631d0, &
       -0.014952258337062199118d0,-0.00030292051472413308126d0, &
       0.0018899503327676891843d0 /
  data cg  / -0.0018899503327676891843d0, &
       -0.00030292051472413308126d0, 0.014952258337062199118d0, &
       0.0038087520138944894631d0, -0.049137179673730286787d0, &
       -0.027219029917103486322d0, 0.051945838107881800736d0, &
       0.36444189483617893676d0, -0.77718575169962802862d0, &
       0.48135965125905339159d0, 0.061273359067811077843d0, &
       -0.14329423835127266284d0, -0.0076074873249766081919d0, &
       0.031695087811525991431d0, 0.00054213233180001068935d0, &
       -0.0033824159510050025955d0  /

  do j=1,ndat
     do i=0,n-1
        ci=0.d0
        di=0.d0
        do l=-7,8
           ci=ci+ch(l)*x(l+2*i,j)
           di=di+cg(l)*x(l+2*i,j)
        enddo
        y(j,i)=ci
        y(j,n+i)=di
     enddo
  enddo

END SUBROUTINE analysis_free_ref
EOF
    kernel.procedure = p
    ConvolutionGenerator::set_lang(lang)
    return kernel
  end
  def ConvolutionGenerator::analysis_per_ref
    lang = ConvolutionGenerator::get_lang
    ConvolutionGenerator::set_lang(ConvolutionGenerator::FORTRAN)
    kernel = CKernel::new
    kernel.lang = ConvolutionGenerator::FORTRAN
    function_name = "analysis_per_ref"
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
subroutine analysis_per_ref(n,ndat,x,y)
  use module_base
  implicit none
  integer, intent(in) :: n,ndat
  real(wp), dimension(0:2*n-1,ndat), intent(in) :: x
  real(wp), dimension(ndat,0:2*n-1), intent(out) :: y
  !local variables
  integer :: i,j,k,l
  real(wp) :: ci,di
  real(wp), dimension(-7:8) :: ch,cg
  !       Daubechy S16
  data ch  /  -0.0033824159510050025955_wp, & 
       -0.00054213233180001068935_wp, 0.031695087811525991431_wp, & 
       0.0076074873249766081919_wp, -0.14329423835127266284_wp, & 
       -0.061273359067811077843_wp, 0.48135965125905339159_wp,  & 
       0.77718575169962802862_wp,0.36444189483617893676_wp, &
       -0.051945838107881800736_wp,-0.027219029917103486322_wp, &
       0.049137179673730286787_wp,0.0038087520138944894631_wp, &
       -0.014952258337062199118_wp,-0.00030292051472413308126_wp, &
       0.0018899503327676891843_wp /

!cg(l)=(-1)**l ch(1-l)
  data cg  / -0.0018899503327676891843_wp, &
       -0.00030292051472413308126_wp, 0.014952258337062199118_wp, &
       0.0038087520138944894631_wp, -0.049137179673730286787_wp, &
       -0.027219029917103486322_wp, 0.051945838107881800736_wp, &
        0.36444189483617893676_wp, -0.77718575169962802862_wp, &
       0.48135965125905339159_wp, 0.061273359067811077843_wp, &
       -0.14329423835127266284_wp, -0.0076074873249766081919_wp, &
       0.031695087811525991431_wp, 0.00054213233180001068935_wp, &
       -0.0033824159510050025955_wp  /

  do j=1,ndat

     do i=0,n-1
        ci=0.e0_wp
        di=0.e0_wp
        do l=-7,8
           k=modulo(l+2*i,2*n)
            ci=ci+ch(l)*x(k    ,j)
            di=di+cg(l)*x(k    ,j)
        enddo
        y(j,i)=ci
        y(j,n+i)=di
     enddo

  enddo
END SUBROUTINE analysis_per_ref
EOF
    kernel.procedure = p
    ConvolutionGenerator::set_lang(lang)
    return kernel
  end

  def ConvolutionGenerator::Analysis(filt, center, unroll, free=false )
    kernel = CKernel::new
    ConvolutionGenerator::set_output( kernel.code )
    kernel.lang = ConvolutionGenerator::get_lang
    function_name = "analysis"
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
    lowfil = Variable::new("lowfil",Int,{:constant => -center}) #-7 for daub16
    upfil = Variable::new("upfil",Int,{:constant => filt.length - center - 1}) # 8 for daub16

    if free then
      dim_in_min = lowfil #-7 in free BC
      dim_in_max = n*2 -2 + upfil #2*n+6 in free BC
    else
      dim_in_min = 0
      dim_in_max = n*2-1
    end
    dim_out_min = 0
    dim_out_max = n*2-1

    #potentially not needed
    #    if free then
    #      lowlimit=lowfil-1 #(dim_out_min-1)/2 #0 in periodic, -4 in free BC
    #      uplimit=n-2+upfil #(dim_out_max-1)/2 #n-1 in periodic, n+2 in free BC
    #    else
    #      lowlimit=0
    #      uplimit=n-2
    #    end

    x = Variable::new("x",Real,{:direction => :in, :dimension => [ Dimension::new(dim_in_min, dim_in_max), Dimension::new(ndat) ] })
    y = Variable::new("y",Real,{:direction => :out, :dimension => [ Dimension::new(ndat), Dimension::new(dim_out_min, dim_out_max) ] })
    i = Variable::new("i",Int)
    j = Variable::new("j",Int)
    k = Variable::new("k",Int)
    l = Variable::new("l",Int)
    ci = [Variable::new("ci1",Real)]
    2.upto(unroll) { |index|
      ci.push(Variable::new("ci#{index}",Real))
    }
    di = [Variable::new("di1",Real)]
    2.upto(unroll) { |index|
      di.push(Variable::new("di#{index}",Real))
    }
    arr = ConstArray::new(filt.reverse,Real)

    fil = Variable::new("fil",Real,{:constant => arr,:dimension => [ Dimension::new(-center, -center -1 + filt.length) ]})


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
      ci.each{ |s| s.decl }
      di.each{ |s| s.decl }

      $output.print("!$omp parallel default (private) shared(x,y,fil,ndat,n)\n") if ConvolutionGenerator::get_lang == ConvolutionGenerator::FORTRAN
      $output.print("!$omp do\n") if ConvolutionGenerator::get_lang == ConvolutionGenerator::FORTRAN

      #external loop, with unrolling. the rest is excluded
      if unroll > 0 then
        forJ = For::new(j,1,ndat-(unroll-1), unroll)
      else
        forJ = For::new(j,1,ndat)
      end
      #print the do part
      forJ.print
      if free then
        For::new(i,0,n-1) {
          ci.each{ |s| (s === 0.0).print }
          di.each{ |s| (s === 0.0).print }
          For::new(l,lowfil,upfil,2) {
            ci.each_index{ |ind|
              (ci[ind] === ci[ind] + fil[l]*x[l + i*2,j+ind]).print
              (di[ind] === di[ind] - fil[-l +1 ]*x[l + i*2,j+ind]).print

              (ci[ind] === ci[ind] + fil[l+1]*x[l+1 + i*2,j+ind]).print
              (di[ind] === di[ind] + fil[-l]*x[l+1 + i*2,j+ind]).print
            }
          }.unroll
          ci.each_index { |ind|
            (y[j+ind,i] === ci[ind]).print
          }
          di.each_index { |ind|
            (y[j+ind,n+i] === di[ind]).print
          }
        }.print
      end      
      #end do for the external loop
      forJ.close
      $output.print("!$omp end do\n") if ConvolutionGenerator::get_lang == ConvolutionGenerator::FORTRAN
      
      #remaining part after unrolling
      if unroll>1 then
        $output.print("!$omp do\n") if ConvolutionGenerator::get_lang == ConvolutionGenerator::FORTRAN
        for1 = For::new(j,ndat-FuncCall::new("modulo",ndat,unroll)+1,ndat) {
          if free then
            For::new(i,0,n-1) {
              ind=0
              (ci[ind] === 0.0).print 
              (di[ind] === 0.0).print
              For::new(l,lowfil,upfil,2) {
                (ci[ind] === ci[ind] + fil[l]*x[l + i*2,j+ind]).print
                (di[ind] === di[ind] - fil[-l+1]*x[l + i*2,j+ind]).print
                
                (ci[ind] === ci[ind] + fil[l+1]*x[l+1 + i*2,j+ind]).print
                (di[ind] === di[ind] + fil[-l]*x[l+1 + i*2,j+ind]).print
                
              }.unroll
              (y[j+ind,i] === ci[ind]).print
              (y[j+ind,n+i] === di[ind]).print
            }.print
          end      
        }.print
        $output.print("!$omp end do\n") if ConvolutionGenerator::get_lang == ConvolutionGenerator::FORTRAN
      end
      $output.print("!$omp end parallel\n") if ConvolutionGenerator::get_lang == ConvolutionGenerator::FORTRAN
    }
    p.print
    kernel.procedure = p
    return kernel
  end
end

FILTER = ["0.0018899503327676891843",
          "-0.00030292051472413308126",
          "-0.014952258337062199118",
          "0.0038087520138944894631",
          "0.049137179673730286787",
          "-0.027219029917103486322",
          "-0.051945838107881800736",
          "0.36444189483617893676",
          "0.77718575169962802862",
          "0.48135965125905339159",
          "-0.061273359067811077843",
          "-0.14329423835127266284",
          "0.0076074873249766081919",
          "0.031695087811525991431",
          "-0.00054213233180001068935",
          "-0.0033824159510050025955"]

n1 = 124
n2 = 132
n3 = 130
input = NArray.float(n1+14,n2,n3).random
output_ref = NArray.float(n2,n3,n1)
output = NArray.float(n2,n3,n1)
epsilon = 10e-15
ConvolutionGenerator::set_lang( ConvolutionGenerator::FORTRAN )
k = ConvolutionGenerator::analysis_free_ref
stats = k.run(n1/2, n2*n3, input, output_ref)
puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"

(0..8).each{ |unroll|
  k = ConvolutionGenerator::Analysis(FILTER,7,unroll,true)
#  k.print
  k.build({:FC => 'ifort',:CC => 'icc',:FCFLAGS => "-O2 -openmp",:LDFLAGS => "-openmp"})
  stats = k.run(n1/2, n2*n3, input, output)
  stats = k.run(n1/2, n2*n3, input, output)
  diff = (output_ref - output).abs
  diff.each { |elem|
    puts "Warning: residue too big: #{elem}" if elem > epsilon
  }
  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
}
ConvolutionGenerator::set_lang( ConvolutionGenerator::C )
(0..8).each{ |unroll|
  k = ConvolutionGenerator::Analysis(FILTER,7,unroll,true)

  stats = k.run(n1/2, n2*n3, input, output)
  stats = k.run(n1/2, n2*n3, input, output)
  diff = (output_ref - output).abs
  diff.each { |elem|
    puts "Warning: residue too big: #{elem}" if elem > epsilon
  }
  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
}

