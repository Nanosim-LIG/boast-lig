module BOAST

  require './compute_strain_product_helper.rb'

  def BOAST::compute_element_strain_undo_att(n_gllx = 5, n_gll2 = 25, n_gll3 = 125, n_gll3_padded = 128)

    function_name = "compute_element_strain_undo_att"

    ngllx = Int("NGLLX", :const => n_gllx)
    ngll2 = Int("NGLL2", :const => n_gll2)
    ngll3 = Int("NGLL3", :const => n_gll3)
    ngll3_padded = Int("NGLL3_PADDED", :const => n_gll3_padded)

    v = []
    v.push ispec = Int("ispec", :dir => :in)
    v.push ijk_ispec = Int("ijk_ispec", :dir => :in)
    v.push d_ibool   = Int( "d_ibool",  :dir => :in, :dim => [Dim()] ) #unused
    v.push *s_dummy_loc = ["x", "y", "z"].collect { |a|
      Real("s_dummy#{a}_loc", :dir => :in, :dim => [Dim(ngll3)], :local => true )
    }
    v.push *d_xi =    [d_xix      = Real("d_xix",    :dir => :in, :dim => [Dim()] ), d_xiy    = Real("d_xiy",   :dir => :in, :dim => [Dim()] ), d_xiz    = Real("d_xiz",   :dir => :in, :dim => [Dim()] ) ]
    v.push *d_eta =   [d_etax     = Real("d_etax",                  :dir => :in, :dim => [Dim()] ), d_etay = Real("d_etay",:dir => :in, :dim => [Dim()] ), d_etaz = Real("d_etaz",:dir => :in, :dim => [Dim()] ) ]
    v.push *d_gamma = [d_gammax   = Real("d_gammax",                :dir => :in, :dim => [Dim()] ), d_gammay = Real("d_gammay",:dir => :in, :dim => [Dim()] ), d_gammaz = Real("d_gammaz",:dir => :in, :dim => [Dim()] ) ]
    v.push sh_hprime_xx = Real("sh_hprime_xx", :dir => :in, :dim => [Dim(ngll2)], :local => true )
    v.push epsilondev_loc     = Real("epsilondev_loc",          :dir => :out, :dim => [Dim(5)], :register => true )
    v.push epsilon_trace_over_3  = Real("epsilon_trace_over_3",          :dir => :out, :dim => [Dim(1)], :register => true )


    sub = Procedure(function_name, v, [], :local => true) {
      decl tx = Int("tx")
      decl k = Int("K") 
      decl j = Int("J") 
      decl i = Int("I")
      decl l = Int("l")
      decl offset = Int("offset")
      tempanl = ["x", "y", "z"].collect { |a|
        [ 1, 2, 3 ].collect { |n|
          Real("temp#{a}#{n}l")
        }
      }
      decl *(tempanl.flatten)
      decl *xil   = [ Real("xixl"),   Real("xiyl"),   Real("xizl")   ]
      decl *etal  = [ Real("etaxl"),  Real("etayl"),  Real("etazl")  ]
      decl *gammal= [ Real("gammaxl"),Real("gammayl"),Real("gammazl")]
      dudl = ["x", "y", "z"].collect { |a1|
        ["x", "y", "z"].collect { |a2|
          Real("du#{a1}d#{a2}l")
        }
      }
      decl *(dudl.flatten)
      decl templ= Real("templ")
      decl *fac = (1..3).collect { |n| Real("fac#{n}") }

      print tx === get_local_id(0)
      print k === tx/ngll2
      print j === (tx-k*ngll2)/ngllx
      print i === tx - k*ngll2 - j*ngllx

      tempanl.flatten.each { |t|
        print t === 0.0
      }
      print For(l, 0, ngllx - 1) {
        print fac[0] === sh_hprime_xx[l*ngllx + i]
        (0..2).each { |indx|
          print tempanl[indx][0] === tempanl[indx][0] + s_dummy_loc[indx][k*ngll2 + j*ngllx + l]*fac[0]
        }
        print fac[1] === sh_hprime_xx[l*ngllx + j]
        (0..2).each { |indx|
          print tempanl[indx][1] === tempanl[indx][1] + s_dummy_loc[indx][k*ngll2 + l*ngllx + i]*fac[1]
        }
        print fac[2] === sh_hprime_xx[l*ngllx + k]
        (0..2).each { |indx|
          print tempanl[indx][2] === tempanl[indx][2] + s_dummy_loc[indx][l*ngll2 + j*ngllx + i]*fac[2]
        }
      }
      print offset === ispec*ngll3_padded + tx
      (0..2).each { |indx|
        print xil[indx] === d_xi[indx][offset]
        print etal[indx] === d_eta[indx][offset]
        print gammal[indx] === d_gamma[indx][offset]
      }
      (0..2).each { |indx1|
        (0..2).each { |indx2|
          print dudl[indx1][indx2] === xil[indx2]*tempanl[indx1][0] + etal[indx2]*tempanl[indx1][1] + gammal[indx2]*tempanl[indx1][2]
        }
      }
      print templ === (dudl[0][0] + dudl[1][1] + dudl[2][2]) * 0.33333333333333333333
      print epsilondev_loc[0] === dudl[0][0] - templ;
      print epsilondev_loc[1] === dudl[1][1] - templ;
      print epsilondev_loc[2] === ( dudl[0][1] + dudl[1][0] ) * 0.5;
      print epsilondev_loc[3] === ( dudl[2][0] + dudl[0][2] ) * 0.5;
      print epsilondev_loc[4] === ( dudl[2][1] + dudl[1][2] ) * 0.5;
      print epsilon_trace_over_3.dereference === templ
    }
    return sub
  end

  def BOAST::compute_ani_undo_att_kernel(ref = true, n_gllx = 5, n_gll2 = 25,  n_gll3 = 125)
    push_env( :array_start => 0 )
    kernel = CKernel::new
    function_name = "compute_ani_undo_att_kernel"
    v = []
    v.push epsilondev_xx          = Real("epsilondev_xx",          :dir => :in, :dim => [Dim()] )
    v.push epsilondev_yy          = Real("epsilondev_yy",          :dir => :in, :dim => [Dim()] )
    v.push epsilondev_xy          = Real("epsilondev_xy",          :dir => :in, :dim => [Dim()] )
    v.push epsilondev_xz          = Real("epsilondev_xz",          :dir => :in, :dim => [Dim()] )
    v.push epsilondev_yz          = Real("epsilondev_yz",          :dir => :in, :dim => [Dim()] )
    v.push epsilon_trace_over_3   = Real("epsilon_trace_over_3",   :dir => :in, :dim => [Dim()] )
    v.push cijkl_kl               = Real("cijkl_kl",               :dir => :inout,:dim => [Dim(21),Dim()] )
    v.push nspec                  = Int( "NSPEC",                  :dir => :in)
    v.push deltat                 = Real("deltat",                 :dir => :in)
    v.push d_ibool                = Int( "d_ibool",                :dir => :in, :dim => [Dim()] )
    v.push d_b_displ              = Real("d_b_displ",              :dir => :in, :dim => [Dim(3), Dim()] )
    v.push *d_xi =    [d_xix      = Real("d_xix",    :dir => :in, :dim => [Dim()] ), d_xiy    = Real("d_xiy",   :dir => :in, :dim => [Dim()] ), d_xiz    = Real("d_xiz",   :dir => :in, :dim => [Dim()] ) ]
    v.push *d_eta =   [d_etax     = Real("d_etax",                  :dir => :in, :dim => [Dim()] ), d_etay = Real("d_etay",:dir => :in, :dim => [Dim()] ), d_etaz = Real("d_etaz",:dir => :in, :dim => [Dim()] ) ]
    v.push *d_gamma = [d_gammax   = Real("d_gammax",                :dir => :in, :dim => [Dim()] ), d_gammay = Real("d_gammay",:dir => :in, :dim => [Dim()] ), d_gammaz = Real("d_gammaz",:dir => :in, :dim => [Dim()] ) ]
    v.push d_hprime_xx             = Real("d_hprime_xx",             :dir => :in, :dim => [Dim()] )

    epsilondev = [ epsilondev_xx, epsilondev_yy, epsilondev_xy, epsilondev_xz, epsilondev_yz ]

    ngllx = Int("NGLLX", :const => n_gllx)
    ngll2 = Int("NGLL2", :const => n_gll2)
    ngll3 = Int("NGLL3", :const => n_gll3)
    ngll3_padded = Int("NGLL3_PADDED", :const => n_gll3_padded)

    p = Procedure(function_name, v)
    if(get_lang == CUDA and ref) then
      @@output.print File::read("references/#{function_name}.cu")
    elsif(get_lang == CL or get_lang == CUDA) then
      make_specfem3d_header( :ngllx => n_gllx, :ngll2 => n_gll2, :ngll3 => n_gll3, :ngll3_padded => n_gll3_padded )
      sub_compute_strain_product =  compute_strain_product()
      print sub_compute_strain_product
      sub_compute_element_strain_undo_att = compute_element_strain_undo_att(n_gllx, n_gll2, n_gll3, n_gll3_padded )
      print sub_compute_element_strain_undo_att
      decl p
        decl i = Int("i")
        decl ispec = Int("ispec")
        decl ijk_ispec = Int("ijk_ispec")
        decl tx = Int("tx")
        decl iglob = Int("iglob")

        decl eps_trace_over_3 = Real("eps_trace_over_3")
        decl b_eps_trace_over_3 = Real("b_eps_trace_over_3")
        decl prod = Real("prod", :dim => [Dim(21)], :allocate => true)
        decl epsdev = Real("epsdev", :dim => [Dim(5)], :allocate => true)
        decl b_epsdev = Real("b_epsdev", :dim => [Dim(5)], :allocate => true)

        decl *s_dummy_loc = ["x", "y", "z"].collect { |a|
          Real("s_dummy#{a}_loc", :local => true, :dim => [Dim(ngll3)] )
        }
        decl sh_hprime_xx = Real("sh_hprime_xx",     :local => true, :dim => [Dim(ngll2)] )
        
        print ispec === get_group_id(0) + get_group_id(1)*get_num_groups(0)
        print ijk_ispec === get_local_id(0) + ngll3*ispec
        print tx === get_local_id(0)

        print If(tx < ngll2) {
          print sh_hprime_xx[tx] === d_hprime_xx[tx]
        }

        print If(ispec < nspec) {
          print iglob === d_ibool[ijk_ispec] - 1
          (0..2).each { |indx|
            print s_dummy_loc[indx][tx] === d_b_displ[iglob*3+indx]
          }
        }
        print barrier(:local)

        print If(ispec < nspec) {
          (0..4).each { |indx|
            print epsdev[indx] === epsilondev[indx][ijk_ispec]
          }

          print eps_trace_over_3 === epsilon_trace_over_3[ijk_ispec]

          print sub_compute_element_strain_undo_att.call(ispec,ijk_ispec,\
                                                         d_ibool,\
                                                         s_dummy_loc[0],s_dummy_loc[1],s_dummy_loc[2],\
                                                         d_xix,d_xiy,d_xiz,d_etax,d_etay,d_etaz,d_gammax,d_gammay,d_gammaz,\
                                                         sh_hprime_xx,\
                                                         b_epsdev,b_eps_trace_over_3.address)
          
          print sub_compute_strain_product.call(prod,eps_trace_over_3,epsdev,b_eps_trace_over_3,b_epsdev)
          print For(i, 0, 21-1) {
            print cijkl_kl[i + ijk_ispec*21] === cijkl_kl[i + ijk_ispec*21] + deltat * prod[i]
          }
        }


      close p
    else
      raise "Unsupported language!"
    end
    pop_env( :array_start )
    kernel.procedure = p
    return kernel
  end
end

