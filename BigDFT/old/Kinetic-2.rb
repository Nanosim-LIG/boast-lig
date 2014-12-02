require 'BOAST'
require 'rubygems'
require 'narray'
require "./Kinetic.rb"
module ConvolutionGenerator

  def ConvolutionGenerator::kinetic_1d(filt, center, unroll, ekin = false, free=false,mod_arr=false)
    kernel = CKernel::new
    ConvolutionGenerator::set_output( kernel.code )
    kernel.lang = ConvolutionGenerator::get_lang
    function_name = "kinetic_1d"
    if free
    then function_name +="_f"
    else function_name +="_p"
    end
    if unroll > 0 then
      function_name += "_#{unroll}"
    end

    n = Variable::new("n",Int,{:direction => :in, :signed => false})
    ndat = Variable::new("ndat",Int,{:direction => :in, :signed => false})
    scal = Variable::new("scal",Real,{:direction => :in})
    c = Variable::new("c",Real,{:direction => :in})
    ek = Variable::new("kstrten",Real,{:direction => :out}) if ekin
    x_in = Variable::new("x_in",Real,{:direction => :in, :dimension => [ Dimension::new(0, n-1),  Dimension::new(ndat)] })
    y_in = Variable::new("y_in",Real,{:direction => :in, :dimension => [ Dimension::new(0, n-1), Dimension::new(ndat)] })
    x_out = Variable::new("x_out",Real,{:direction => :out, :dimension => [ Dimension::new(ndat), Dimension::new(0, n-1)] })
    #y_out = Variable::new("y_out",Real,{:direction => :out, :dimension => [ Dimension::new(0, n-1),  Dimension::new(ndat)] })
    y_out = Variable::new("y_out",Real,{:direction => :out, :dimension => [ Dimension::new(ndat), Dimension::new(0, n-1)] })

    if ConvolutionGenerator::get_lang == C then
      @@output.print "inline #{Int::new.decl} modulo( #{Int::new.decl} a, #{Int::new.decl} b) { return (a+b)%b;}\n"
      @@output.print "inline #{Int::new.decl} min( #{Int::new.decl} a, #{Int::new.decl} b) { return a < b ? a : b;}\n"
      @@output.print "inline #{Int::new.decl} max( #{Int::new.decl} a, #{Int::new.decl} b) { return a > b ? a : b;}\n"
    end

    lowfil = Variable::new("lowfil",Int,{:constant => -center})
    upfil = Variable::new("upfil",Int,{:constant => filt.length - center -1})
    arr = ConstArray::new(filt,Real)
    fil = Variable::new("fil",Real,{:constant => arr,:dimension => [ Dimension::new((-center),(filt.length - center -1)) ]})
    i = Variable::new("i",Int)
    j = Variable::new("j",Int)
    l = Variable::new("l",Int)
    k = Variable::new("k",Int)

    #to replace modulo operation
    mods=nil
    mods=Variable::new("mod_arr",Real,{:dimension => [ Dimension::new((-center),n+(filt.length - center -1)) ]}) if mod_arr and not free     

    tt = [Variable::new("tt1",Real)]
    2.upto(unroll) { |index|
      tt.push(Variable::new("tt#{index}",Real))
    }

    #define the internal convolution. Fundamental operation which is associated to the kernel
    convolution_op = lambda { |t,border|
      if border then
        (k ===  i+l - ((i+ l  +  n * 2 )/n - 2) *n  ).print unless mods
        ii=(mods ? mods[i+l] : k )
      else
        ii=i+l
      end
      t.each_index { |ind|
          jj=j+ind
        (t[ind] === t[ind] + x_in[ ii , jj ]*fil[l] ).print
        }
    }

    #this is the sum over all the filters 
    for_uniq = lambda { |t,border|
      t.each_index { |ind|
        (t[ind] === 0.0).print
      }

      if (border) then
        if( free ) then
          For::new( l, FuncCall::new( "max", -i, lowfil), FuncCall::new( "min", upfil, n - 1 - i ) ) {
            convolution_op.call(t,false) #not free 
          }.unroll#.print#
        else
          For::new( l, lowfil, upfil) {
            convolution_op.call(t,true) #not free
          }.unroll#.print#
        end
      else
        For::new( l, lowfil, upfil) {
          convolution_op.call(t,false) 
        }.unroll#.print#
      end
      t.each_index { |ind|
        (t[ind] === t[ind] * scal).print
        (eks === eks + t[ind] * x_in[i,j+ind]).print if ekin
        (y_out[j+ind,i] ===  y_in[i,j+ind] + t[ind]).print
        #(y_out[i,j+ind] ===  y_in[i,j+ind] + t[ind]).print
        (x_out[j+ind,i] ===  x_in[i,j+ind]).print
      }
    }


    kinetic_1d = lambda { |t,js|
      unro=t.length
      @@output.print("!$omp do\n") if ConvolutionGenerator::get_lang == ConvolutionGenerator::FORTRAN
      For::new(j,js,(ndat/unro)*unro, step: unro) {
        For::new(i, 0, -lowfil) {
          for_uniq.call(t,true)
        }.print
        For::new(i, -lowfil+1, n-1 - upfil) {
          for_uniq.call(t,false)
        }.print
        For::new(i, n - upfil, n-1) {
          for_uniq.call(t, true)
        }.print
      }.print
      @@output.print("!$omp end do\n") if ConvolutionGenerator::get_lang == ConvolutionGenerator::FORTRAN
    }

    vars = [n,ndat,scal,x_in,y_in,x_out,y_out]
    vars += [ek] if ekin
    p = Procedure::new(function_name, vars, [lowfil,upfil]){
      i.decl
      j.decl
      l.decl
      k.decl unless mods
      fil.decl
      ek.decl if ekin
      tt.each { |t| t.decl }
      mods.decl if mod_arr

      (ek === 0.0).print if ekin

      if mod_arr then
          For::new(j,(-center),n+(filt.length - center -1)) {
            (mods[j] === FuncCall::new( "modulo", j, n)).print
          }.print
      end

      
      if ConvolutionGenerator::get_lang == ConvolutionGenerator::FORTRAN then
        @@output.print("!$omp parallel default(shared)&\n")
        @@output.print("!$omp reduction(+:kstrten)&\n") if ekin
        @@output.print("!$omp private(i,j,k,l)&\n")
        @@output.print("!$omp private(")
        @@output.print(tt.join(","))
        @@output.print(")\n")
      end
      kinetic_1d.call(tt,1) 
      kinetic_1d.call([tt[0]],(ndat/unroll)*unroll+1) if unroll > 1

      @@output.print("!$omp  end parallel\n") if ConvolutionGenerator::get_lang == ConvolutionGenerator::FORTRAN
    }
    p.print
    kernel.procedure = p
    return kernel
  end

end




