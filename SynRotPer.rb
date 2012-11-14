require './Algorithm.rb'
module ConvolutionGenerator
#subroutine syn_rot_per_temp(n,ndat,x,y)
#  use module_base
#  implicit none
#  integer, intent(in) :: n,ndat
#  real(wp), dimension(0:2*n-1,ndat), intent(in) :: x
#  real(wp), dimension(ndat,0:2*n-1), intent(out) :: y
#  !local variables
#  integer :: i,j,k,l
#  real(wp) :: so,se
#  !       Daubechy S16
#  real(wp), dimension(-6:9), parameter :: ch=(/&
#                                    0.0018899503327676891843_wp, &
#       -0.00030292051472413308126_wp,-0.014952258337062199118_wp, &
#       0.0038087520138944894631_wp, 0.049137179673730286787_wp, &
#       -0.027219029917103486322_wp, -0.051945838107881800736_wp, &
#       0.36444189483617893676_wp, 0.77718575169962802862_wp, &
#       0.48135965125905339159_wp, -0.061273359067811077843_wp, &
#       -0.14329423835127266284_wp, 0.0076074873249766081919_wp, &
#       0.031695087811525991431_wp, -0.00054213233180001068935_wp, &
#       -0.0033824159510050025955_wp/)
#  real(wp), dimension(-6:9), parameter :: cg=(/&
#                                    -0.0033824159510050025955_wp, & 
#       0.00054213233180001068935_wp, 0.031695087811525991431_wp, & 
#       -0.0076074873249766081919_wp, -0.14329423835127266284_wp, & 
#       0.061273359067811077843_wp, 0.48135965125905339159_wp,  & 
#       -0.77718575169962802862_wp,0.36444189483617893676_wp, &
#       0.051945838107881800736_wp,-0.027219029917103486322_wp, &
#       -0.049137179673730286787_wp,0.0038087520138944894631_wp, &
#       0.014952258337062199118_wp,-0.00030292051472413308126_wp, &
#       -0.0018899503327676891843_wp  /)
#  
#  do j=1,ndat
#
#     !fifth changement
#     so=0.0_wp
#     se=0.0_wp
#     do l=-3,4
#        k=modulo(n-1+l,n)
#        se=se+ch(2*l)*x(k,j)+cg(2*l)*x(n+k,j)
#        so=so+ch(2*l+1)*x(k,j)+cg(2*l+1)*x(n+k  ,j)
#     end do
#     y(j,2*n-1)=so
#     y(j,0  )=se
#
#     do i=0,n-2
#        so=0.0_wp
#        se=0.0_wp
#        do l=-3,4
#           k=modulo(i+l,n)
#           se=se+ch(2*l)*x(k,j)+cg(2*l)*x(n+k,j)
#           so=so+ch(2*l+1)*x(k,j)+cg(2*l+1)*x(n+k ,j)
#        end do
#        y(j,2*i+1)=so
#        y(j,2*i+2)=se
#     end do
#  end do
#
#END SUBROUTINE syn_rot_per_temp

  def ConvolutionGenerator::Synthesis(filt, center, unroll, free=false )
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
    
    p = Procedure::new(function_name, [n,ndat,x,y], [lowfil,upfil]) {
      
      i.decl
      j.decl
      k.decl
      l.decl
      fil.decl
      so.each{ |s| s.decl }
      se.each{ |s| s.decl }

      $output.print("!$omp parallel default (private) shared(x,y,fil,ndat,n)\n")
      $output.print("!$omp do\n")

      #internal loop taking care of BCs
      if free then
        forBC = For::new(l, FuncCall::new("max", -i,lowfil), FuncCall::new("min", upfil, n-1-i) ){
          (k === i+l).print
          se.each_index{ |ind|
            (se[ind] === se[ind] + fil[l*2]*x[k,j+ind] + fil[l*-2+3]*x[n+k,j+ind]).print
            (so[ind] === so[ind] + fil[l*2+1]*x[k,j+ind] - fil[l*-2+2]*x[n+k,j+ind]).print
          }
        }
      else
        forBC = For::new(l,lowfil,upfil){
          (k === FuncCall::new( "modulo", i+l, n)).print
          se.each_index{ |ind|
            (se[ind] === se[ind] + fil[l*2]*x[k,j+ind] + fil[l*-2+3]*x[n+k,j+ind]).print
            (so[ind] === so[ind] + fil[l*2+1]*x[k,j+ind] - fil[l*-2+2]*x[n+k,j+ind]).print
          }
        }
      end


      if unroll > 0 then
        for1 = For::new(j,1,ndat-(unroll-1), unroll)
      else
        for1 = For::new(j,1,ndat)
      end
      for1.print
      if not free then
        so.each{ |s| (s === 0.0).print }
        se.each{ |s| (s === 0.0).print }
        (i === n-1).print
        forBC.unroll
        so.each_index { |ind|
          (y[j+ind,n*2-1] === so[ind]).print
        }
        se.each_index { |ind|
          (y[j+ind,0] === se[ind]).print
        }
      end

      For::new(i,lowlimit,-lowfil-1) {
        so.each{ |s| (s === 0.0).print }
        se.each{ |s| (s === 0.0).print }
        if free then
          forBC.print
        else
          forBC.unroll
        end

        so.each_index { |ind|
          (y[j+ind,i*2+1] === so[ind]).print
        }
        se.each_index { |ind|
          (y[j+ind,i*2+2] === se[ind]).print
        }
      }.print
      
      For::new(i,-lowfil,n-1-upfil) {
        so.each{ |s| (s === 0.0).print }
        se.each{ |s| (s === 0.0).print }
        For::new(l,lowfil,upfil){
          (k === i + l).print
          se.each_index{ |ind|
            (se[ind] === se[ind] + fil[l*2]*x[k,j+ind] + fil[l*-2+3]*x[n+k,j+ind]).print
            (so[ind] === so[ind] + fil[l*2+1]*x[k,j+ind] - fil[l*-2+2]*x[n+k,j+ind]).print
          }
        }.unroll
        so.each_index { |ind|
          (y[j+ind,i*2+1] === so[ind]).print
        }
        se.each_index { |ind|
          (y[j+ind,i*2+2] === se[ind]).print
        }
      }.print

      For::new(i,n-upfil,uplimit) {
        so.each{ |s| (s === 0.0).print }
        se.each{ |s| (s === 0.0).print }
        if free then
          forBC.print
        else
          forBC.unroll
        end
        so.each_index { |ind|
          (y[j+ind,i*2+1] === so[ind]).print
        }
        se.each_index { |ind|
          (y[j+ind,i*2+2] === se[ind]).print
        }
      }.print
      for1.close
      $output.print("!$omp end do\n")


      if unroll>1 then
        $output.print("!$omp do\n")
        for1 = For::new(j,ndat-FuncCall::new("modulo",ndat,unroll)+1,ndat) {
          ind=0
          if not free then
            (so[ind] === 0.0).print
            (se[ind] === 0.0).print
            
            (i === n-1).print
            for3 = For::new(l,lowfil,upfil){
              (k === FuncCall::new( "modulo", i+l, n)).print
              (se[ind] === se[ind] + fil[l*2]*x[k,j+ind] + fil[l*-2+3]*x[n+k,j+ind]).print
              (so[ind] === so[ind] + fil[l*2+1]*x[k,j+ind] - fil[l*-2+2]*x[n+k,j+ind]).print
            }
            for3.unroll
            (y[j+ind,n*2-1] === so[ind]).print
            (y[j+ind,0] === se[ind]).print
          end
          
          for2 = For::new(i,lowlimit,uplimit) {
            (so[ind] === 0.0).print
            (se[ind] === 0.0).print
            if free then
              For::new(l, FuncCall::new("max", -i,lowfil), FuncCall::new("min", upfil, n-1-i) ){
                (k === i+l).print
                (se[ind] === se[ind] + fil[l*2]*x[k,j+ind] + fil[l*-2+3]*x[n+k,j+ind]).print
                (so[ind] === so[ind] + fil[l*2+1]*x[k,j+ind] - fil[l*-2+2]*x[n+k,j+ind]).print
              }.print
            else
              for3.unroll
            end
            (y[j+ind,i*2+1] === so[ind]).print
            (y[j+ind,i*2+2] === se[ind]).print
          }.print
          
        }.print
        $output.print("!$omp end do\n")
      end
      $output.print("!$omp end parallel\n")
    }.print
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

ConvolutionGenerator::set_lang( ConvolutionGenerator::FORTRAN )
(0..8).each{ |unroll|
  ConvolutionGenerator::Synthesis(FILTER,7,unroll,true)
}
(0..8).each{ |unroll|
  ConvolutionGenerator::Synthesis(FILTER,7,unroll,false)
}
