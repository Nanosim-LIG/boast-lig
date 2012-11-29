require "./BOAST.rb"
require 'narray'
module ConvolutionGenerator
  def ConvolutionGenerator::magicfilter_per_ref
    lang = ConvolutionGenerator::get_lang
    ConvolutionGenerator::set_lang(ConvolutionGenerator::FORTRAN)
    kernel = CKernel::new
    kernel.lang = ConvolutionGenerator::FORTRAN
    function_name = "magicfilter_per_ref"
    n = Variable::new("n",Int,{:direction => :in, :signed => false})
    ndat = Variable::new("ndat",Int,{:direction => :in, :signed => false})
    dim_in_min = 0
    dim_in_max = n-1
    dim_out_min = 0
    dim_out_max = n-1
    x = Variable::new("x",Real,{:direction => :in, :dimension => [ Dimension::new(dim_in_min, dim_in_max), Dimension::new(ndat) ] })
    y = Variable::new("y",Real,{:direction => :out, :dimension => [ Dimension::new(ndat), Dimension::new(dim_out_min, dim_out_max) ] })
    p = Procedure::new(function_name, [n,ndat,x,y])
    kernel.code.print <<EOF
subroutine magicfilter_per_ref(n1,ndat,x,y)
  implicit none
  integer, parameter :: wp=kind(1.0d0)
  integer, intent(in) :: n1,ndat
  real(kind=8), dimension(0:(n1-1),ndat), intent(in) :: x
  real(kind=8), dimension(ndat,0:(n1-1)), intent(out) :: y
  !local variables
  integer, parameter :: lowfil=-8,lupfil=7
  integer :: i,j,k,l
  real(kind=8) :: tt

  !          THE MAGIC FILTER FOR DAUBECHIES-16
  real(kind=8) fil(lowfil:lupfil)
  DATA fil / &
       8.4334247333529341094733325815816e-7_wp,&
       -0.1290557201342060969516786758559028e-4_wp,&
       0.8762984476210559564689161894116397e-4_wp,&
       -0.30158038132690463167163703826169879e-3_wp,&
       0.174723713672993903449447812749852942e-2_wp,&
       -0.942047030201080385922711540948195075e-2_wp,&
       0.2373821463724942397566389712597274535e-1_wp,&
       0.612625895831207982195380597e-1_wp,&
       0.9940415697834003993178616713_wp,&
       -0.604895289196983516002834636e-1_wp, &
       -0.2103025160930381434955489412839065067e-1_wp,&
       0.1337263414854794752733423467013220997e-1_wp,&
       -0.344128144493493857280881509686821861e-2_wp,&
       0.49443227688689919192282259476750972e-3_wp,&
       -0.5185986881173432922848639136911487e-4_wp,&
       2.72734492911979659657715313017228e-6_wp /
  do j=1,ndat
     do i=0,n1-1
        tt=0.e0_wp
        do l=lowfil,lupfil
           k=modulo(i+l,n1)  
           tt=tt+x(  k,j)*fil(l)
        enddo
        y(j,i)=tt
     enddo
  enddo

END SUBROUTINE magicfilter_per_ref
EOF
    kernel.procedure = p
    ConvolutionGenerator::set_lang(lang)
    return kernel
  end

  def ConvolutionGenerator::MagicFilter(filt, center, unroll, invert, free=false )
    kernel = CKernel::new
    ConvolutionGenerator::set_output( kernel.code )
    kernel.lang = ConvolutionGenerator::get_lang
    function_name = "magicfilter"
    if free then
      function_name += "_free"
    else 
      function_name += "_per"
    end
    if invert then
      function_name += "_inv"
    end
    if unroll>0 then
      function_name += "_u#{unroll}"
    end

    n = Variable::new("n",Int,{:direction => :in, :signed => false})
    ndat = Variable::new("ndat",Int,{:direction => :in, :signed => false})
    if invert then
      lowfil = Variable::new("lowfil",Int,{:constant => 1-center})
      upfil = Variable::new("upfil",Int,{:constant => filt.length-center})
    else
      lowfil = Variable::new("lowfil",Int,{:constant => center-filt.length})
      upfil = Variable::new("upfil",Int,{:constant => center-1})
    end

    dim_in_min = 0
    dim_in_max = n-1
    dim_out_min = 0
    dim_out_max = n-1
    if free then
      if invert then
        dim_in_min = lowfil
        dim_in_max = n-1+upfil
      else
        dim_out_min = -upfil
        dim_out_max = n-1-lowfil
      end
    end

    x = Variable::new("x",Real,{:direction => :in, :dimension => [ Dimension::new(dim_in_min, dim_in_max), Dimension::new(ndat) ] })
    y = Variable::new("y",Real,{:direction => :out, :dimension => [ Dimension::new(ndat), Dimension::new(dim_out_min, dim_out_max) ] })
    i = Variable::new("i",Int)
    j = Variable::new("j",Int)
    k = Variable::new("k",Int)
    l = Variable::new("l",Int)
    tt = [Variable::new("tt1",Real)]
    2.upto(unroll) { |index|
      tt.push(Variable::new("tt#{index}",Real))
    }
    if invert then filt = filt.reverse end 
    arr = ConstArray::new(filt,Real)
    fil = Variable::new("fil",Real,{:constant => arr,:dimension => [ Dimension::new(lowfil,upfil) ]})

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
      tt.each { |t|
        t.decl
      }
  
      fil.decl
  
      if unroll>1 then 
        for1 = For::new(j,1,ndat-(unroll-1), unroll)
      else
        for1 = For::new(j,1,ndat)
      end
      for1.print
      for2 = For::new(i,dim_out_min, dim_out_max) {
        tt.each{ |t|
          (t === 0.0).print
        }
        if free and not invert
          for3 = For::new(l, FuncCall::new("max", -i, lowfil), FuncCall::new("min", upfil, n-1-i) )
        else
          for3 = For::new(l, lowfil, upfil)
        end
        for3.print
        if not free then
          (k === FuncCall::new( "modulo", i+l, n)).print
        else
          (k === i+l).print
        end
        tt.each_index{ |index|
          (tt[index] === tt[index] + x[k,j+index]*fil[l]).print
        }
        for3.close
        tt.each_index{ |index|
          (y[j+index,i] === tt[index]).print
        }
      }.print
      for1.close
  
      if unroll>1 then
        for1 = For::new(j,ndat-FuncCall::new("modulo",ndat,unroll)+1,ndat) {
          for2 = For::new(i,dim_out_min,dim_out_max) {
            (tt[0] === 0.0).print
            if free and not invert then
              for3 = For::new(l, FuncCall::new("max", -i, lowfil), FuncCall::new("min", upfil, n-1-i) ) {
                (k === i+l).print
                (tt[0] === tt[0] + x[k,j] * fil[l] ).print
              }.print
            else
              for3 = For::new(l, lowfil, upfil) {
                if not free then
                  (k === FuncCall::new( "modulo", i+l, n)).print
                else
                  (k === i+l).print
                end
                (tt[0] === tt[0] + x[k,j] * fil[l] ).print
              }.print
            end
            (y[j,i] === tt[0]).print
          }.print
        }.print
      end
    }
    p.print
    kernel.procedure = p
    return kernel
  end
end

FILTER = [ "8.4334247333529341094733325815816e-7",
       "-0.1290557201342060969516786758559028e-4",
       "0.8762984476210559564689161894116397e-4",
       "-0.30158038132690463167163703826169879e-3",
       "0.174723713672993903449447812749852942e-2",
       "-0.942047030201080385922711540948195075e-2",
       "0.2373821463724942397566389712597274535e-1",
       "0.612625895831207982195380597e-1",
       "0.9940415697834003993178616713",
       "-0.604895289196983516002834636e-1",
       "-0.2103025160930381434955489412839065067e-1",
       "0.1337263414854794752733423467013220997e-1",
       "-0.344128144493493857280881509686821861e-2",
       "0.49443227688689919192282259476750972e-3",
       "-0.5185986881173432922848639136911487e-4",
       "2.72734492911979659657715313017228e-6" ]



n1 = 126
n2 = 128
n3 = 130
input = NArray.float(n1,n2,n3).random
output_ref = NArray.float(n1,n2,n3)
output = NArray.float(n1,n2,n3)
epsilon = 10e-15



ConvolutionGenerator::set_lang( ConvolutionGenerator::FORTRAN )
ConvolutionGenerator::magicfilter_per_ref.run(n1, n2*n3, input, output_ref)
#ConvolutionGenerator::MagicFilter(FILTER,8,0,false).run(32, 32, " "*32*32*8, " "*32*32*8)
1.upto(15) { |i|
  ConvolutionGenerator::MagicFilter(FILTER,8,i,false).run(n1, n2*n3, input, output)
  diff = (output_ref - output).abs
  diff.each { |elem|
    puts "Warning: residue too big: #{elem}" if elem > epsilon
  }
}

ConvolutionGenerator::MagicFilter(FILTER,8,5,true).build
ConvolutionGenerator::MagicFilter(FILTER,8,3,false,true).build
ConvolutionGenerator::MagicFilter(FILTER,8,4,true,true).build
ConvolutionGenerator::set_lang( ConvolutionGenerator::C )
ConvolutionGenerator::MagicFilter(FILTER,8,0,false).run(32, 32, NArray.float(32*32), NArray.float(32*32))
ConvolutionGenerator::MagicFilter(FILTER,8,0,false).run(32, 32, " "*32*32*8, " "*32*32*8)
ConvolutionGenerator::MagicFilter(FILTER,8,0,false).build
ConvolutionGenerator::MagicFilter(FILTER,8,5,true).build
ConvolutionGenerator::MagicFilter(FILTER,8,3,false,true).build
ConvolutionGenerator::MagicFilter(FILTER,8,4,true,true).build

