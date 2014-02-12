module BOAST
  def BOAST::compute_element_oc_rotation( n_gll3 = 125)
    function_name = "compute_element_oc_rotation"
    v = []
    v.push tx                    = Int( "tx",                    :dir => :in)
    v.push working_element       = Int( "working_element",       :dir => :in)
    v.push time                  = Real("time",                  :dir => :in)
    v.push two_omega_earth       = Real("two_omega_earth",       :dir => :in)
    v.push deltat                = Real("deltat",                :dir => :in)
    v.push d_A_array_rotation    = Real("d_A_array_rotation",    :dir => :inout, :dim => [Dim()] )
    v.push d_B_array_rotation    = Real("d_B_array_rotation",    :dir => :inout, :dim => [Dim()] )
    v.push dpotentialdxl         = Real("dpotentialdxl",         :dir => :in)
    v.push dpotentialdyl         = Real("dpotentialdyl",         :dir => :in)
    v.push dpotentialdx_with_rot = Real("dpotentialdx_with_rot", :dir => :out, :dim => [Dim(1)], :register => true)
    v.push dpotentialdy_with_rot = Real("dpotentialdy_with_rot", :dir => :out, :dim => [Dim(1)], :register => true)

    ngll3 = Int("NGLL3", :const => n_gll3)

    p = Procedure(function_name, v, [], :local => true) {
      decl two_omega_deltat = Real("two_omega_deltat")
      decl cos_two_omega_t  = Real("cos_two_omega_t")
      decl sin_two_omega_t  = Real("sin_two_omega_t")
      decl a_rotation = Real("A_rotation")
      decl b_rotation = Real("b_rotation")
      decl source_euler_A = Real("source_euler_A")
      decl source_euler_B = Real("source_euler_B")

      print cos_two_omega_t === cos(two_omega_earth*time)
      print sin_two_omega_t === sin(two_omega_earth*time)

      print two_omega_deltat === deltat * two_omega_earth

      print source_euler_A === two_omega_deltat * (cos_two_omega_t * dpotentialdyl + sin_two_omega_t * dpotentialdxl)
      print source_euler_B === two_omega_deltat * (sin_two_omega_t * dpotentialdyl - cos_two_omega_t * dpotentialdxl)

      print a_rotation = d_A_array_rotation[tx + working_element*ngll3]
      print b_rotation = d_B_array_rotation[tx + working_element*ngll3]

      print dpotentialdx_with_rot[0] === dpotentialdxl + ( a_rotation*cos_two_omega_t + b_rotation*sin_two_omega_t)
      print dpotentialdy_with_rot[0] === dpotentialdyl + (-a_rotation*sin_two_omega_t + b_rotation*cos_two_omega_t)

      print d_A_array_rotation[tx + working_element*ngll3] === d_A_array_rotation[tx + working_element*ngll3] + source_euler_A
      print d_B_array_rotation[tx + working_element*ngll3] === d_B_array_rotation[tx + working_element*ngll3] + source_euler_B
    }
    return p
  end

  def BOAST::outer_core_impl_kernel(ref = true, mesh_coloring = false, textures_fields = false, textures_constants = false, unroll_loops = true, n_gllx = 5, n_gll2 = 25, n_gll3 = 125, n_gll3_padded = 128, r_earth_km = 6371.0, coloring_min_nspec_outer_core = 1000)
    push_env( :array_start => 0 )
    kernel = CKernel::new
    v = []
    function_name = "outer_core_impl_kernel"
    v.push nb_blocks_to_compute    = Int("nb_blocks_to_compute",     :dir => :in)
    v.push nglob                   = Int("NGLOB",                    :dir => :in)
    v.push d_ibool                 = Int("d_ibool",                  :dir => :in, :dim => [Dim()] )
    v.push d_phase_ispec_inner     = Int("phase_ispec_inner",        :dir => :in, :dim => [Dim()] )
    v.push num_phase_ispec         = Int("num_phase_ispec",          :dir => :in)
    v.push d_iphase                = Int("d_iphase",                 :dir => :in)
    v.push use_mesh_coloring_gpu   = Int("use_mesh_coloring_gpu",    :dir => :in)
    v.push d_potential             = Real("d_potential",             :dir => :in, :dim => [Dim()] )
    v.push d_potential_dot_dot     = Real("d_potential_dot_dot",     :dir => :inout, :dim => [Dim()] )
    v.push *d_xi = [d_xix          = Real("d_xix",                   :dir => :in, :dim => [Dim()] ), d_xiy = Real("d_xiy",:dir => :in, :dim => [Dim()] ), d_xiz = Real("d_xiz",:dir => :in, :dim => [Dim()] ) ]
    v.push *d_eta = [d_etax        = Real("d_etax",                  :dir => :in, :dim => [Dim()] ), d_etay = Real("d_etay",:dir => :in, :dim => [Dim()] ), d_etaz = Real("d_etaz",:dir => :in, :dim => [Dim()] ) ]
    v.push *d_gamma = [d_gammax    = Real("d_gammax",                :dir => :in, :dim => [Dim()] ), d_gammay = Real("d_gammay",:dir => :in, :dim => [Dim()] ), d_gammaz = Real("d_gammaz",:dir => :in, :dim => [Dim()] ) ]
    v.push d_hprime_xx             = Real("d_hprime_xx",             :dir => :in, :dim => [Dim()] )
    v.push d_hprimewgll_xx         = Real("d_hprimewgll_xx",         :dir => :in, :dim => [Dim()] )
    v.push wgllwgll_xy             = Real("wgllwgll_xy",             :dir => :in, :dim => [Dim()] )
    v.push wgllwgll_xz             = Real("wgllwgll_xz",             :dir => :in, :dim => [Dim()] )
    v.push wgllwgll_yz             = Real("wgllwgll_yz",             :dir => :in, :dim => [Dim()] )
    v.push gravity                 = Int( "GRAVITY",                 :dir => :in)
    v.push *d_store = [d_xstore    = Real("d_xstore",                :dir => :in, :dim => [Dim()] ), d_ystore = Real("d_ystore",:dir => :in, :dim => [Dim()] ), d_zstore = Real("d_zstore",:dir => :in, :dim => [Dim()] ) ]
    v.push d_d_ln_density_dr_table = Real("d_d_ln_density_dr_table", :dir => :in, :dim => [Dim()] )
    v.push d_minus_rho_g_over_kappa_fluid = Real("d_minus_rho_g_over_kappa_fluid", :dir => :in, :dim => [Dim()] )  
    v.push wgll_cube               = Real("wgll_cube",               :dir => :in, :dim => [Dim()] )
    v.push rotation                = Int( "ROTATION",                :dir => :in)
    v.push time                    = Real("time",                    :dir => :in)
    v.push two_omega_earth         = Real("two_omega_earth",         :dir => :in)
    v.push deltat                  = Real("deltat",                  :dir => :in)
    v.push d_A_array_rotation      = Real("d_A_array_rotation",      :dir => :inout, :dim => [Dim()] )
    v.push d_B_array_rotation      = Real("d_B_array_rotation",      :dir => :inout, :dim => [Dim()] )
    v.push nspec_outer_core        = Int( "NSPEC_OUTER_CORE",        :dir => :in)

    ngllx        = Int("NGLLX", :const => n_gllx)
    ngll2        = Int("NGLL2", :const => n_gll2)
    ngll3        = Int("NGLL3", :const => n_gll3)
    ngll3_padded = Int("NGLL3_PADDED", :const => n_gll3_padded)

    constants = [] #ngllx, ngll2, ngll3, ngll3_padded]

    if textures_fields then
      d_displ_oc_tex = Real("d_displ_oc_tex", :texture => true, :dir => :in, :dim => [Dim()] )
      d_accel_oc_tex = Real("d_accel_oc_tex", :texture => true, :dir => :in, :dim => [Dim()] )
      if get_lang == CL then
        v.push(d_displ_oc_tex, d_accel_oc_tex)
        constants.push( d_displ_oc_tex.sampler, d_accel_oc_tex.sampler )
      end
    end
    if textures_constants then
      d_hprime_xx_oc_tex = Real("d_hprime_xx_oc_tex", :texture => true, :dir => :in, :dim => [Dim()] )
      if get_lang == CL then
        v.push(d_hprime_xx_oc_tex)
        constants.push( d_hprime_xx_oc_tex.sampler )
      end
    end

    p = Procedure(function_name, v, constants)
    if(get_lang == CUDA and ref) then
      @@output.print File::read("specfem3D/#{function_name}.cu")
    elsif(get_lang == CL or get_lang == CUDA) then
      make_specfem3d_header(:ngllx => n_gllx, :ngll2 => n_gll2, :ngll3 => n_gll3, :ngll3_padded => n_gll3_padded, :r_earth_km => r_earth_km, :coloring_min_nspec_outer_core => coloring_min_nspec_outer_core)
      if get_lang == CUDA
        if textures_fields then
          decl d_displ_oc_tex
          decl d_accel_oc_tex
        end
        if textures_constants then
          decl d_hprime_xx_oc_tex
        end
      end
      sub_kernel =  compute_element_oc_rotation(n_gll3)
      print sub_kernel
      decl p
      decl bx = Int("bx")
      decl tx = Int("tx")
      decl k  = Int("K"), j = Int("J"), i = Int("I")
      decl active = Int("active"), offset = Int("offset")
      decl iglob  = Int("iglob")
      decl working_element = Int("working_element")
      decl l = Int("l")

      decl *templ = [ Real("temp1l"), Real("temp2l"), Real("temp3l") ]
      decl *xil   = [ Real("xixl"),   Real("xiyl"),   Real("xizl")   ]
      decl *etal  = [ Real("etaxl"),  Real("etayl"),  Real("etazl")  ]
      decl *gammal= [ Real("gammaxl"),Real("gammayl"),Real("gammazl")]
      decl jacobianl= Real("jacobianl")
      decl *dpotentialdl = [ Real("dpotentialdxl"), Real("dpotentialdyl"), Real("dpotentialdzl") ]
      decl dpotentialdx_with_rot = Real("dpotentialdx_with_rot")
      decl dpotentialdy_with_rot = Real("dpotentialdy_with_rot")

      decl sum_terms = Real("sum_terms")
      decl gravity_term = Real("gravity_term")
      decl *gl   = [ Real("gxl"),   Real("gyl"),   Real("gzl")   ]
      
      decl radius = Real("radius"), theta = Real("theta"), phi = Real("phi")
      decl cos_theta = Real("cos_theta"), sin_theta = Real("sin_theta")
      decl cos_phi   = Real("cos_phi"),   sin_phi   = Real("sin_phi")
      decl *grad_ln_rho = [ Real("grad_x_ln_rho"), Real("grad_y_ln_rho"), Real("grad_z_ln_rho") ]

      decl int_radius = Int("int_radius")

      decl s_dummy_loc = Real("s_dummy_loc", :shared => true, :dim => [Dim(ngll3)] )

      decl s_temp1 = Real("s_temp1", :shared => true, :dim => [Dim(ngll3)])
      decl s_temp2 = Real("s_temp2", :shared => true, :dim => [Dim(ngll3)])
      decl s_temp3 = Real("s_temp3", :shared => true, :dim => [Dim(ngll3)])

      decl sh_hprime_xx     = Real("sh_hprime_xx",     :shared => true, :dim => [Dim(ngll2)] )
      decl sh_hprimewgll_xx = Real("sh_hprimewgll_xx", :shared => true, :dim => [Dim(ngll2)] )

      print bx === get_group_id(1)*get_num_groups(0)+get_group_id(0)
      print tx === get_local_id(0)

      print k === tx/ngll2
      print j === (tx-k*ngll2)/ngllx
      print i === tx - k*ngll2 - j*ngllx

      print active === Ternary( Expression("&&", tx < ngll3, bx < nb_blocks_to_compute), 1, 0)
      print If( active ) {
        if mesh_coloring then
          print working_element === bx
        else
          print If( use_mesh_coloring_gpu, lambda {
            print working_element === bx
          }, lambda {
            print working_element === d_phase_ispec_inner[bx + num_phase_ispec*(d_iphase-1)]-1
          })
        end

        print iglob === d_ibool[working_element*ngll3 + tx]-1

        if textures_fields then
          print s_dummy_loc[tx] === d_displ_oc_tex[iglob]
        else
          print s_dummy_loc[tx] === d_potential[iglob]
        end
      }
      print If(tx < ngll2) {
        if textures_constants then
          print sh_hprime_xx[tx] === d_hprime_xx_oc_tex[tx]
        else
          print sh_hprime_xx[tx] === d_hprime_xx[tx]
        end

        print sh_hprimewgll_xx[tx] === d_hprimewgll_xx[tx]
      }
      print barrier(:local)

      print If( active ) {
        (0..2).each { |indx| print templ[indx] === 0.0 }
        for_loop = For(l, 0, ngllx-1) {
           print templ[0] === templ[0] + s_dummy_loc[k*ngll2+j*ngllx+l]*sh_hprime_xx[l*ngllx+i]
           print templ[1] === templ[1] + s_dummy_loc[k*ngll2+l*ngllx+i]*sh_hprime_xx[l*ngllx+j]
           print templ[2] === templ[2] + s_dummy_loc[l*ngll2+j*ngllx+i]*sh_hprime_xx[l*ngllx+k]
        }
        if unroll_loops then
          for_loop.unroll
        else
          print for_loop
        end
        print offset === working_element*ngll3_padded + tx
        (0..2).each { |indx|
          print xil[indx]    === d_xi[indx][offset]
          print etal[indx]   === d_eta[indx][offset]
          print gammal[indx] === d_gamma[indx][offset]
        }

        print jacobianl === Expression("/", 1.0, xil[0]*(etal[1]*gammal[2] - etal[2]*gammal[1]) - xil[1]*(etal[0]*gammal[2] - etal[2]*gammal[0]) + xil[2]*(etal[0]*gammal[1] - etal[1]*gammal[0]))
        (0..2).each { |indx|
          print dpotentialdl[indx] === xil[indx]*templ[0] + etal[indx]*templ[1] + gammal[indx]*templ[2]
        }

        print If( rotation , lambda {
          print sub_kernel.call(tx,working_element,time,two_omega_earth,deltat,d_A_array_rotation,d_B_array_rotation,dpotentialdl[0],dpotentialdl[1],dpotentialdx_with_rot.address, dpotentialdy_with_rot.address)
        }, lambda {
          print dpotentialdx_with_rot === dpotentialdl[0]
          print dpotentialdy_with_rot === dpotentialdl[1]
        })

        print radius === d_store[0][iglob]
        print theta  === d_store[1][iglob]
        print phi    === d_store[2][iglob]
        if(get_lang == CL) then
          print sin_theta === sincos(theta, cos_theta.address)
          print sin_phi   === sincos(phi,   cos_phi.address)
        else
          if(get_default_real_size == 4) then
            print sincosf(theta, sin_theta.address, cos_theta.address)
            print sincosf(phi,   sin_phi.address,   cos_phi.address)
          else
            print cos_theta === cos(theta)
            print sin_theta === sin(theta)
            print cos_phi   === cos(phi)
            print sin_phi   === sin(phi)
          end
        end
        print int_radius === rint(radius * r_earth_km * 10.0) - 1
        If( !gravity , lambda {
          print grad_ln_rho[0] === sin_theta * cos_phi * d_d_ln_density_dr_table[int_radius]
          print grad_ln_rho[1] === sin_theta * sin_phi * d_d_ln_density_dr_table[int_radius]
          print grad_ln_rho[2] ===           cos_theta * d_d_ln_density_dr_table[int_radius]

          print dpotentialdx_with_rot === dpotentialdx_with_rot + s_dummy_loc[tx] * grad_ln_rho[0]
          print dpotentialdy_with_rot === dpotentialdy_with_rot + s_dummy_loc[tx] * grad_ln_rho[1]
          print dpotentialdl[2]       === dpotentialdl[2] +       s_dummy_loc[tx] * grad_ln_rho[2]
        }, lambda {
          print gl[0] === sin_theta*cos_phi
          print gl[1] === sin_theta*sin_phi
          print gl[2] === cos_theta

          print gravity_term === d_minus_rho_g_over_kappa_fluid[int_radius] * jacobianl * wgll_cube[tx] * \
                                 (dpotentialdx_with_rot*gl[0] + dpotentialdy_with_rot*gl[1] + dpotentialdl[2]*gl[2])
        })

        print s_temp1[tx] === jacobianl*(   xil[0]*dpotentialdx_with_rot +    xil[1]*dpotentialdy_with_rot +    xil[2]*dpotentialdl[2])
        print s_temp1[tx] === jacobianl*(  etal[0]*dpotentialdx_with_rot +   etal[1]*dpotentialdy_with_rot +   etal[2]*dpotentialdl[2])
        print s_temp1[tx] === jacobianl*(gammal[0]*dpotentialdx_with_rot + gammal[1]*dpotentialdy_with_rot + gammal[2]*dpotentialdl[2])
      }
      print barrier(:local)
      print If( active ){
        (0..2).each { |indx| print templ[indx] === 0.0 }
        for_loop = For(l, 0, ngllx-1) {
           print templ[0] === templ[0] + s_temp1[k*ngll2+j*ngllx+l]*sh_hprimewgll_xx[i*ngllx+l]
           print templ[1] === templ[1] + s_temp2[k*ngll2+l*ngllx+i]*sh_hprimewgll_xx[j*ngllx+l]
           print templ[2] === templ[2] + s_temp3[l*ngll2+j*ngllx+i]*sh_hprimewgll_xx[k*ngllx+l]
        }
        if unroll_loops then
          for_loop.unroll
        else
          print for_loop
        end
        print sum_terms === -(wgllwgll_yz[k*ngllx+j]*templ[0] + wgllwgll_xz[k*ngllx+i]*templ[1] + wgllwgll_xy[j*ngllx+i]*templ[2])

        print If( gravity ) {
          print  sum_terms === sum_terms + gravity_term
        }
        if mesh_coloring then
          if textures_fields then
            print d_potential_dot_dot[iglob] === d_accel_oc_tex[iglob] + sum_terms
          else
            print d_potential_dot_dot[iglob] === d_potential_dot_dot[iglob] + sum_terms
          end
        else
          print If( use_mesh_coloring_gpu, lambda {
            print If( nspec_outer_core > coloring_min_nspec_outer_core, lambda {
              if textures_fields then
                print d_potential_dot_dot[iglob] === d_accel_oc_tex[iglob] + sum_terms
              else
                print d_potential_dot_dot[iglob] === d_potential_dot_dot[iglob] + sum_terms
              end
            }, lambda{
              print atomicAdd(d_potential_dot_dot+iglob,sum_terms)
            })
          }, lambda {
            print atomicAdd(d_potential_dot_dot+iglob,sum_terms)
          })
        end
      }
      close p
    else
      raise "Unsupported language!"
    end
    pop_env(:array_start)
    kernel.procedure = p
    return kernel
  end
end
