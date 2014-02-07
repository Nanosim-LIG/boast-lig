module BOAST
  def BOAST::compute_element_ic_att_stress( n_gll3 = 125, n_sls = 3 )
    function_name = "compute_element_ic_att_stress"
    v = []
    v.push tx              = Int( "tx",              :dir => :in)
    v.push working_element = Int( "working_element", :dir => :in)
    v.push r_xx            = Real("R_xx",            :dir => :in, :dim => [Dim()] )
    v.push r_yy            = Real("R_yy",            :dir => :in, :dim => [Dim()] )
    v.push r_xy            = Real("R_xy",            :dir => :in, :dim => [Dim()] )
    v.push r_xz            = Real("R_xz",            :dir => :in, :dim => [Dim()] )
    v.push r_yz            = Real("R_yz",            :dir => :in, :dim => [Dim()] )
    v.push sigma_xx        = Real("sigma_xx",        :dir => :inout, :dim => [Dim()] )
    v.push sigma_yy        = Real("sigma_yy",        :dir => :inout, :dim => [Dim()] )
    v.push sigma_zz        = Real("sigma_zz",        :dir => :inout, :dim => [Dim()] )
    v.push sigma_xy        = Real("sigma_xy",        :dir => :inout, :dim => [Dim()] )
    v.push sigma_xz        = Real("sigma_xz",        :dir => :inout, :dim => [Dim()] )
    v.push sigma_yz        = Real("sigma_yz",        :dir => :inout, :dim => [Dim()] )

    ngll3 = Int("NGLL3", :const => n_gll3)
    nsls  = Int("N_SLS", :const => n_sls)

    p = Procedure(function_name, v, [ngll3,nsls]) {
      decl offset = Int("offset")
      decl i_sls  = Int("i_sls")
      decl r_xx_val = Real("R_xx_val")
      decl r_yy_val = Real("R_yy_val")
      print For( i_sls, 0, nsls - 1) {
        print offset === i_sls + nsls*(tx + ngll3*working_element)
        print r_xx_val === r_xx[offset]
        print r_yy_val === r_yy[offset]
        sigma_xx[0] === sigma_xx[0] - r_xx_val
        sigma_yy[0] === sigma_yy[0] - r_yy_val
        sigma_zz[0] === sigma_zz[0] + r_xx_val + r_yy_val
        sigma_xy[0] === sigma_xy[0] - r_xy[offset]
        sigma_xz[0] === sigma_xz[0] - r_xz[offset]
        sigma_yz[0] === sigma_yz[0] - r_yz[offset]
      }
    }
    return p
  end

  def BOAST::compute_element_ic_att_memory(n_gll3 = 125, n_gll3_padded = 128, n_sls = 3 )
    function_name = "compute_element_ic_att_memory"
    v = []
    v.push tx                = Int( "tx",                :dir => :in)
    v.push working_element   = Int( "working_element",   :dir => :in)
    v.push d_muv             = Real("d_muv",             :dir => :in, :dim => [Dim()] )
    v.push factor_common     = Real("factor_common",     :dir => :in, :dim => [Dim()] )
    v.push alphaval          = Real("alphaval",          :dir => :in, :dim => [Dim()] )
    v.push betaval           = Real("betaval",           :dir => :in, :dim => [Dim()] )
    v.push gammaval          = Real("gammaval",          :dir => :in, :dim => [Dim()] )
    v.push r_xx              = Real("R_xx",              :dir => :inout, :dim => [Dim()] )
    v.push r_yy              = Real("R_yy",              :dir => :inout, :dim => [Dim()] )
    v.push r_xy              = Real("R_xy",              :dir => :inout, :dim => [Dim()] )
    v.push r_xz              = Real("R_xz",              :dir => :inout, :dim => [Dim()] )
    v.push r_yz              = Real("R_yz",              :dir => :inout, :dim => [Dim()] )
    v.push epsilondev_xx     = Real("epsilondev_xx",     :dir => :in, :dim => [Dim()] )
    v.push epsilondev_yy     = Real("epsilondev_yy",     :dir => :in, :dim => [Dim()] )
    v.push epsilondev_xy     = Real("epsilondev_xy",     :dir => :in, :dim => [Dim()] )
    v.push epsilondev_xz     = Real("epsilondev_xz",     :dir => :in, :dim => [Dim()] )
    v.push epsilondev_yz     = Real("epsilondev_yz",     :dir => :in, :dim => [Dim()] )
    v.push epsilondev_xx_loc = Real("epsilondev_xx_loc", :dir => :in)
    v.push epsilondev_yy_loc = Real("epsilondev_yy_loc", :dir => :in)
    v.push epsilondev_xy_loc = Real("epsilondev_xy_loc", :dir => :in)
    v.push epsilondev_xz_loc = Real("epsilondev_xz_loc", :dir => :in)
    v.push epsilondev_yz_loc = Real("epsilondev_yz_loc", :dir => :in)
    v.push use_3d_attenuation_arrays = Int( "USE_3D_ATTENUATION_ARRAYS",    :dir => :in)

    ngll3 = Int("NGLL3", :const => n_gll3)
    ngll3_padded = Int("NGLL3_PADDED", :const => n_gll3_padded)
    nsls  = Int("N_SLS", :const => n_sls)

    p = Procedure(function_name, v, [ngll3,ngll3_padded, nsls]) {
      decl offset = Int("offset")
      decl i_sls  = Int("i_sls")
      decl mul = Real("mul")
      decl alphaval_loc = Real("alphaval_loc")
      decl betaval_loc = Real("betaval_loc")
      decl gammaval_loc = Real("gammaval_loc")
      decl factor_loc = Real("factor_loc")
      decl sn = Real("sn")
      decl snp1 = Real("snp1")

      print mul === d_muv[tx + ngll3_padded*working_element]
      print For( i_sls, 0, nsls - 1 ) {
        print offset === i_sls + nsls*(tx + ngll3*working_element)
        print If(use_3d_attenuation_arrays, lambda {
            print factor_loc  === mul * factor_common[offset]
        }, lambda {
            print factor_loc  === mul * factor_common[i_sls + nsls*working_element ]
        })
        print alphaval_loc === alphaval[i_sls]
        print  betaval_loc ===  betaval[i_sls]
        print gammaval_loc === gammaval[i_sls]

        [[r_xx, epsilondev_xx, epsilondev_xx_loc],
         [r_yy, epsilondev_yy, epsilondev_yy_loc],
         [r_xy, epsilondev_xy, epsilondev_xy_loc],
         [r_xz, epsilondev_xz, epsilondev_xz_loc],
         [r_yz, epsilondev_yz, epsilondev_yz_loc]].each { |r, epsilondev, epsilondev_loc|
          print sn   === factor_loc * epsilondev[tx + ngll3 * working_element]
          print snp1 === factor_loc * epsilondev_loc
          print r[offset] === alphaval_loc*r[offset] + betaval_loc*sn + gammaval_loc*snp1
        }
      }
    }
    return p
  end

  def BOAST::compute_element_ic_gravity( n_gll3 = 125, r_earth_km = 6371.0)
    function_name = "compute_element_ic_gravity"
    v = []
    v.push tx                = Int( "tx",                :dir => :in)
    v.push working_element   = Int( "working_element",   :dir => :in)
    v.push d_ibool                 = Int("d_ibool",                  :dir => :in, :dim => [Dim()] )
    v.push *d_store = [d_xstore    = Real("d_xstore",                :dir => :in, :dim => [Dim()] ), d_ystore = Real("d_ystore",:dir => :in, :dim => [Dim()] ), d_zstore = Real("d_zstore",:dir => :in, :dim => [Dim()] ) ]
    v.push d_minus_gravity_table   = Real("d_minus_gravity_table",   :dir => :in, :dim => [Dim()] )
    v.push d_minus_deriv_gravity_table = Real("d_minus_deriv_gravity_table", :dir => :in, :dim => [Dim()] )
    v.push d_density_table         = Real("d_density_table",         :dir => :in, :dim => [Dim()] )
    v.push wgll_cube               = Real("wgll_cube",               :dir => :in, :dim => [Dim()] )
    v.push jacobianl               = Real("jacobianl", :dir => :in)
    v.push *s_dummy_loc = ["x", "y", "z"].collect { |a|
                                     Real("s_dummy#{a}_loc", :dir => :in, :shared => true, :dim => [Dim(ngll3)] )
    }
    sigma = ["x", "y", "z"].collect { |a1|
        ["x", "y", "z"].collect { |a2|
                                     Real("sigma_#{a1}#{a2}", :dir => :inout, :dim => [Dim()] )
        }
    }
    v.push sigma[0][0], sigma[1][1], sigma[2][2],
           sigma[0][1], sigma[1][0], sigma[0][2],\
           sigma[2][0], sigma[1][2], sigma[2][1]
    v.push *rho_s_H = [1,2,3].collect {|n|
                                     Real("rho_s_H#{n}", :dir => :inout, :dim => [Dim()] )
    }

    ngll3 = Int("NGLL3", :const => n_gll3)
    p = Procedure(function_name, v, [ngll3]) {
      decl radius = Real("radius"), theta = Real("theta"), phi = Real("phi")
      decl cos_theta = Real("cos_theta"), sin_theta = Real("sin_theta"), cos_phi = Real("cos_phi"), sin_phi = Real("sin_phi")
      decl cos_theta_sq = Real("cos_theta_sq"), sin_theta_sq = Real("sin_theta_sq"), cos_phi_sq = Real("cos_phi_sq"), sin_phi_sq = Real("sin_phi_sq")
      decl minus_g = Real("minus_g"), minus_dg = Real("minus_dg")
      decl rho = Real("rho")
      decl *gl = [ Real("gxl"), Real("gyl"), Real("gzl") ]
      decl minus_g_over_radius = Real("minus_g_over_radius"), minus_dg_plus_g_over_radius = Real("minus_dg_plus_g_over_radius")
      decl hxxl = Real("Hxxl"), hyyl = Real("Hyyl"), hzzl = Real("Hzzl"), hxyl = Real("Hxyl"), hxzl = Real("Hxzl"), hyzl = Real("Hyzl")
      decl *s_l = [ Real("sx_l"), Real("sy_l"), Real("sz_l") ]
      decl factor = Real("factor")
      decl iglob = Int("iglob")
      decl int_radius = Int("int_radius")

      print iglob === d_ibool[working_element*ngll3 + tx]-1
      print radius === d_xstore[iglob]
      print If(radius < ( 100.0 / (r_earth_km*1000.0))) {
        print radius ===  100.0 / (r_earth_km*1000.0)
      }
      print theta === d_ystore[iglob]
      print phi === d_zstore[iglob]
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
      print int_radius === rint(radius * r_earth_km * 10.0 ) - 1
      print If( int_radius < 0 ) {
        print int_radius === 0
      }
      print minus_g  = d_minus_gravity_table[int_radius]
      print minus_dg = d_minus_deriv_gravity_table[int_radius]
      print rho      = d_density_table[int_radius]

      print gl[0] === minus_g*sin_theta*cos_phi
      print gl[1] === minus_g*sin_theta*sin_phi
      print gl[2] === minus_g*cos_theta

      print minus_g_over_radius === minus_g / radius
      print minus_dg_plus_g_over_radius === minus_dg - minus_g_over_radius

      print cos_theta_sq === cos_theta*cos_theta
      print sin_theta_sq === sin_theta*sin_theta
      print cos_phi_sq   === cos_phi*cos_phi
      print sin_phi_sq   === sin_phi*sin_phi

      print hxxl === minus_g_over_radius*(cos_phi_sq*cos_theta_sq + sin_phi_sq) + cos_phi_sq*minus_dg*sin_theta_sq
      print hyyl === minus_g_over_radius*(cos_phi_sq + cos_theta_sq*sin_phi_sq) + minus_dg*sin_phi_sq*sin_theta_sq
      print hzzl === cos_theta_sq*minus_dg + minus_g_over_radius*sin_theta_sq
      print hxyl === cos_phi*minus_dg_plus_g_over_radius*sin_phi*sin_theta_sq
      print hxzl === cos_phi*cos_theta*minus_dg_plus_g_over_radius*sin_theta
      print hyzl === cos_theta*minus_dg_plus_g_over_radius*sin_phi*sin_theta

      (0..2).each { |indx|
        print s_l[indx] === rho * s_dummy_loc[indx][tx]
      }
      print sigma[0][0].dereference === sigma[0][0].dereference + s_l[1]*gl[1] + s_l[2]*gl[2];
      print sigma[1][1].dereference === sigma[1][1].dereference + s_l[0]*gl[0] + s_l[2]*gl[2];
      print sigma[2][2].dereference === sigma[2][2].dereference + s_l[0]*gl[0] + s_l[1]*gl[1];

      print sigma[0][1].dereference === sigma[0][1].dereference - s_l[0] * gl[1];
      print sigma[1][0].dereference === sigma[1][0].dereference - s_l[1] * gl[0];

      print sigma[0][2].dereference === sigma[0][2].dereference - s_l[0] * gl[2];
      print sigma[2][0].dereference === sigma[2][0].dereference - s_l[2] * gl[0];

      print sigma[1][2].dereference === sigma[1][2].dereference - s_l[1] * gl[2];
      print sigma[2][1].dereference === sigma[2][1].dereference - s_l[2] * gl[1];

      print factor === jacobianl * wgll_cube[tx]
      print rho_s_H[0][0] === factor * (s_l[0]*hxxl + s_l[1]*hxyl + s_l[2]*hxzl)
      print rho_s_H[1][0] === factor * (s_l[0]*hxyl + s_l[1]*hyyl + s_l[2]*hyzl)
      print rho_s_H[2][0] === factor * (s_l[0]*hxzl + s_l[1]*hyzl + s_l[2]*hzzl)
    }
    return p
  end

  def BOAST::inner_core_impl_kernel(ref = true, mesh_coloring = false, textures_fields = false, textures_constants = false, unroll_loops = false, n_gllx = 5, n_gll2 = 25, n_gll3 = 125, n_gll3_padded = 128, r_earth_km = 6371.0, coloring_min_nspec_inner_core = 1000, i_flag_in_fictitious_cube = 11, n_sls = 3)
    push_env( :array_start => 0 )
    kernel = CKernel::new
    v = []
    function_name = "outer_core_impl_kernel"
    v.push nb_blocks_to_compute    = Int("nb_blocks_to_compute",     :dir => :in)
    v.push nglob                   = Int("NGLOB",                    :dir => :in)
    v.push d_ibool                 = Int("d_ibool",                  :dir => :in, :dim => [Dim()] )
    v.push d_idoubling             = Int("d_idoubling",              :dir => :in, :dim => [Dim()] )
    v.push d_phase_ispec_inner     = Int("phase_ispec_inner",        :dir => :in, :dim => [Dim()] )
    v.push num_phase_ispec         = Int("num_phase_ispec",          :dir => :in)
    v.push d_iphase                = Int("d_iphase",                 :dir => :in)
    v.push deltat                  = Real("deltat",                  :dir => :in)
    v.push use_mesh_coloring_gpu   = Int("use_mesh_coloring_gpu",    :dir => :in)
    v.push d_displ                 = Real("d_displ",                 :dir => :in, :dim => [Dim(3), Dim()] )
    v.push d_veloc                 = Real("d_veloc",                 :dir => :in, :dim => [Dim(3), Dim()] ) #unused in original code, to remove
    v.push d_accel                 = Real("d_accel",                 :dir => :inout, :dim => [Dim(3), Dim()] )
    v.push *d_xi = [d_xix          = Real("d_xix",                   :dir => :in, :dim => [Dim()] ), d_xiy = Real("d_xiy",:dir => :in, :dim => [Dim()] ), d_xiz = Real("d_xiz",:dir => :in, :dim => [Dim()] ) ]
    v.push *d_eta = [d_etax        = Real("d_etax",                  :dir => :in, :dim => [Dim()] ), d_etay = Real("d_etay",:dir => :in, :dim => [Dim()] ), d_etaz = Real("d_etaz",:dir => :in, :dim => [Dim()] ) ]
    v.push *d_gamma = [d_gammax    = Real("d_gammax",                :dir => :in, :dim => [Dim()] ), d_gammay = Real("d_gammay",:dir => :in, :dim => [Dim()] ), d_gammaz = Real("d_gammaz",:dir => :in, :dim => [Dim()] ) ]
    v.push d_hprime_xx             = Real("d_hprime_xx",             :dir => :in, :dim => [Dim()] )
    v.push d_hprimewgll_xx         = Real("d_hprimewgll_xx",         :dir => :in, :dim => [Dim()] )
    v.push d_wgllwgll_xy           = Real("d_wgllwgll_xy",           :dir => :in, :dim => [Dim()] )
    v.push d_wgllwgll_xz           = Real("d_wgllwgll_xz",           :dir => :in, :dim => [Dim()] )
    v.push d_wgllwgll_yz           = Real("d_wgllwgll_yz",           :dir => :in, :dim => [Dim()] )
    v.push d_kappav                = Real("d_kappav",                :dir => :in, :dim => [Dim()] )
    v.push d_muv                   = Real("d_muv",                   :dir => :in, :dim => [Dim()] )
    v.push compute_and_store_strain= Int( "COMPUTE_AND_STORE_STRAIN",:dir => :in)
    v.push epsilondev_xx           = Real("epsilondev_xx",           :dir => :in, :dim => [Dim()] )
    v.push epsilondev_yy           = Real("epsilondev_yy",           :dir => :in, :dim => [Dim()] )
    v.push epsilondev_xy           = Real("epsilondev_xy",           :dir => :in, :dim => [Dim()] )
    v.push epsilondev_xz           = Real("epsilondev_xz",           :dir => :in, :dim => [Dim()] )
    v.push epsilondev_yz           = Real("epsilondev_yz",           :dir => :in, :dim => [Dim()] )
    v.push epsilon_trace_over_3    = Real("epsilon_trace_over_3",    :dir => :in, :dim => [Dim()] )
    v.push attenuation             = Int( "ATTENUATION",             :dir => :in)
    v.push partial_phys_dispersion_only = Int( "PARTIAL_PHYS_DISPERSION_ONLY", :dir => :in)
    v.push use_3d_attenuation_arrays    = Int( "USE_3D_ATTENUATION_ARRAYS",    :dir => :in)
    v.push one_minus_sum_beta      = Real("one_minus_sum_beta",      :dir => :in, :dim => [Dim()] )
    v.push factor_common           = Real("factor_common",           :dir => :in, :dim => [Dim()] )
    v.push r_xx                    = Real("R_xx",                    :dir => :inout, :dim => [Dim()] )
    v.push r_yy                    = Real("R_yy",                    :dir => :inout, :dim => [Dim()] )
    v.push r_xy                    = Real("R_xy",                    :dir => :inout, :dim => [Dim()] )
    v.push r_xz                    = Real("R_xz",                    :dir => :inout, :dim => [Dim()] )
    v.push r_yz                    = Real("R_yz",                    :dir => :inout, :dim => [Dim()] )
    v.push alphaval                = Real("alphaval",                :dir => :in, :dim => [Dim()] )
    v.push betaval                 = Real("betaval",                 :dir => :in, :dim => [Dim()] )
    v.push gammaval                = Real("gammaval",                :dir => :in, :dim => [Dim()] )
    v.push anisotropy              = Int( "ANISOTROPY",              :dir => :in)
    v.push d_c11store              = Real("d_c11store",              :dir => :in, :dim => [Dim()] )
    v.push d_c12store              = Real("d_c12store",              :dir => :in, :dim => [Dim()] )
    v.push d_c13store              = Real("d_c13store",              :dir => :in, :dim => [Dim()] )
    v.push d_c33store              = Real("d_c33store",              :dir => :in, :dim => [Dim()] )
    v.push d_c44store              = Real("d_c44store",              :dir => :in, :dim => [Dim()] )
    v.push gravity                 = Int( "GRAVITY",                 :dir => :in)
    v.push *d_store = [d_xstore    = Real("d_xstore",                :dir => :in, :dim => [Dim()] ), d_ystore = Real("d_ystore",:dir => :in, :dim => [Dim()] ), d_zstore = Real("d_zstore",:dir => :in, :dim => [Dim()] ) ]
    v.push d_minus_gravity_table   = Real("d_minus_gravity_table",   :dir => :in, :dim => [Dim()] )
    v.push d_minus_deriv_gravity_table = Real("d_minus_deriv_gravity_table", :dir => :in, :dim => [Dim()] )
    v.push d_density_table         = Real("d_density_table",         :dir => :in, :dim => [Dim()] )
    v.push wgll_cube               = Real("wgll_cube",               :dir => :in, :dim => [Dim()] )
    v.push nspec_inner_core_strain_only = Int( "NSPEC_INNER_CORE_STRAIN_ONLY", :dir => :in)
    v.push nspec_inner_core        = Int( "NSPEC_INNER_CORE",        :dir => :in)

    ngllx        = Int("NGLLX", :const => n_gllx)
    ngll2        = Int("NGLL2", :const => n_gll2)
    ngll3        = Int("NGLL3", :const => n_gll3)
    ngll3_padded = Int("NGLL3_PADDED", :const => n_gll3_padded)
    iflag_in_fictitious_cube = Int("IFLAG_IN_FICTITIOUS_CUBE", :const => i_flag_in_fictitious_cube)

    constants = [ngllx, ngll2, ngll3, ngll3_padded]

    if textures_fields then
      d_displ_ic_tex = Real("d_displ_ic_tex", :texture => true, :dir => :in, :dim => [Dim()] )
      d_accel_ic_tex = Real("d_accel_ic_tex", :texture => true, :dir => :in, :dim => [Dim()] )
      if get_lang == CL then
        v.push(d_displ_ic_tex, d_accel_ic_tex)
        constants.push( d_displ_ic_tex.sampler, d_accel_ic_tex.sampler )
      end
    end
    if textures_constants then
      d_hprime_xx_ic_tex = Real("d_hprime_xx_ic_tex", :texture => true, :dir => :in, :dim => [Dim()] )
      if get_lang == CL then
        v.push(d_hprime_xx_ic_tex)
        constants.push( d_hprime_xx_ic_tex.sampler )
      end
    end

    p = Procedure(function_name, v, constants)
    if(get_lang == CUDA and ref) then
      @@output.print File::read("specfem3D/#{function_name}.cu")
    elsif(get_lang == CL or get_lang == CUDA) then
      if (get_lang == CL) then
        if get_default_real_size == 8 then
          @@output.puts "#pragma OPENCL EXTENSION cl_khr_fp64: enable"
          @@output.puts "#pragma OPENCL EXTENSION cl_khr_int64_base_atomics: enable"
        end
      end
      load "./atomicAdd_f.rb"
      if get_lang == CUDA
        if textures_fields then
          decl d_displ_oc_tex
          decl d_accel_oc_tex
        end
        if textures_constants then
          decl d_hprime_xx_oc_tex
        end
      end
      sub_kernel1 =  compute_element_ic_att_stress(n_gll3, n_sls)
      print sub_kernel1
      sub_kernel2 =  compute_element_ic_att_memory(n_gll3, n_gll3_padded, n_sls)
      print sub_kernel2
      sub_kernel3 =  compute_element_ic_gravity(n_gll3, r_earth_km)
      print sub_kernel3
      decl p
      decl bx = Int("bx")
      decl tx = Int("tx")
      decl k  = Int("K"), j = Int("J"), i = Int("I")
      decl active = Int("active"), offset = Int("offset")
      decl iglob  = Int("iglob")
      decl working_element = Int("working_element")
      decl l = Int("l")
      tempanl = ["x", "y", "z"].collect { |a|
        [ 1, 2, 3 ].collect { |n|
          Real("temp#{a}#{n}l")
        }
      }
      decl *(tempanl.flatten)
      decl *xil   = [ Real("xixl"),   Real("xiyl"),   Real("xizl")   ]
      decl *etal  = [ Real("etaxl"),  Real("etayl"),  Real("etazl")  ]
      decl *gammal= [ Real("gammaxl"),Real("gammayl"),Real("gammazl")]
      decl jacobianl= Real("jacobianl")
      dudl = ["x", "y", "z"].collect { |a1|
        ["x", "y", "z"].collect { |a2|
          Real("du#{a1}d#{a2}l")
        }
      }
      decl *(dudl.flatten)
      decl duxdxl_plus_duydyl = Real("duxdxl_plus_duydyl")
      decl duxdxl_plus_duzdzl = Real("duxdxl_plus_duzdzl")
      decl duydyl_plus_duzdzl = Real("duydyl_plus_duzdzl")
      decl duxdyl_plus_duydxl = Real("duxdyl_plus_duydxl")
      decl duzdxl_plus_duxdzl = Real("duzdxl_plus_duxdzl")
      decl duzdyl_plus_duydzl = Real("duzdyl_plus_duydzl")
      decl templ = Real("templ")
      decl *fac = [1,2,3].collect {|n|
        Real("fac#{n}")
      }
      decl lambdal = Real("lambdal")
      decl mul = Real("mul")
      decl lambdalplus2mul = Real("lambdalplus2mul")
      decl kappal = Real("kappal")
      decl mul_iso = Real("mul_iso")
      decl mul_aniso = Real("mul_aniso")
      sigma = ["x", "y", "z"].collect { |a1|
        ["x", "y", "z"].collect { |a2|
          Real("sigma_#{a1}#{a2}")
        }
      }
      decl *(sigma.flatten)
      decl epsilondev_xx_loc = Real("epsilondev_xx_loc")
      decl epsilondev_yy_loc = Real("epsilondev_yy_loc")
      decl epsilondev_xy_loc = Real("epsilondev_xy_loc")
      decl epsilondev_xz_loc = Real("epsilondev_xz_loc")
      decl epsilondev_yz_loc = Real("epsilondev_yz_loc")
      decl c11 = Real("c11")
      decl c12 = Real("c12")
      decl c13 = Real("c13")
      decl c33 = Real("c33")
      decl c44 = Real("c44")
      decl *sum_terms = [1,2,3].collect {|n|
        Real("sum_terms#{n}")
      }
      decl *rho_s_H = [1,2,3].collect {|n|
        Real("rho_s_H#{n}")
      }

      decl *s_dummy_loc = ["x", "y", "z"].collect { |a|
        Real("s_dummy#{a}_loc", :shared => true, :dim => [Dim(ngll3)] )
      }

      s_temp = ["x", "y", "z"].collect { |a|
        [ 1, 2, 3 ].collect { |n|
          Real("s_temp#{a}#{n}", :shared => true, :dim => [Dim(ngll3)] )
        }
      }
      decl *(s_temp.flatten)

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

        print If( d_idoubling[working_element] == iflag_in_fictitious_cube, lambda {
          print active === 0
        }, lambda {
          print iglob === d_ibool[working_element*ngll3 + tx]-1

          if textures_fields then
            (0..2).each { |indx|
              print s_dummy_loc[indx][tx] === d_displ_ic_tex[iglob*3+indx]
            }
          else
            (0..2).each { |indx|
              print s_dummy_loc[indx][tx] === d_displ[indx, iglob]
            }
          end
        })
      }
      print If(tx < ngll2) {
        if textures_constants then
          print sh_hprime_xx[tx] === d_hprime_xx_ic_tex[tx]
        else
          print sh_hprime_xx[tx] === d_hprime_xx[tx]
        end

        print sh_hprimewgll_xx[tx] === d_hprimewgll_xx[tx]
      }
      print barrier(:local)

      print If( active ) {
        (0..2).each { |indx1|
          (0..2).each { |indx2|
            print tempanl[indx1][indx2] === 0.0
          }
        }
        for_loop = For(l, 0, ngllx-1) {
          print fac[0] === sh_hprime_xx[l*ngllx + i]
          print offset === k*ngll2 + j*ngllx + l
          (0..2).each { |indx1|
            print tempanl[indx1][0] === tempanl[indx1][0] + s_dummy_loc[indx1][offset]*fac[0]
          }
          print fac[1] === sh_hprime_xx[l*ngllx + j]
          print offset === k*ngll2 + l*ngllx + i
          (0..2).each { |indx1|
            print tempanl[indx1][1] === tempanl[indx1][1] + s_dummy_loc[indx1][offset]*fac[1]
          }
          print fac[2] === sh_hprime_xx[l*ngllx + k]
          print offset === l*ngll2 + j*ngllx + i
          (0..2).each { |indx1|
            print tempanl[indx1][2] === tempanl[indx1][2] + s_dummy_loc[indx1][offset]*fac[2]
          }
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

        (0..2).each { |indx1|
          (0..2).each { |indx2|
            print dudl[indx1][indx2] === xil[indx2]*tempanl[indx1][0] + etal[indx2]*tempanl[indx1][1] + gammal[indx2]*tempanl[indx1][2]
          }
        }
        print duxdxl_plus_duydyl === dudl[0][0] + dudl[1][1]
        print duxdxl_plus_duzdzl === dudl[0][0] + dudl[2][2]
        print duydyl_plus_duzdzl === dudl[1][1] + dudl[2][2]
        print duxdyl_plus_duydxl === dudl[0][1] + dudl[1][0]
        print duzdxl_plus_duxdzl === dudl[2][0] + dudl[0][2]
        print duzdyl_plus_duydzl === dudl[2][1] + dudl[1][2]

        print If(compute_and_store_strain) {
          print templ === (dudl[0][0] + dudl[1][1] + dudl[2][2])*0.33333333333333333333333333
          print epsilondev_xx_loc === dudl[0][0] - templ
          print epsilondev_yy_loc === dudl[1][1] - templ
          print epsilondev_xy_loc === duxdyl_plus_duydxl * 0.5
          print epsilondev_xz_loc === duzdxl_plus_duxdzl * 0.5
          print epsilondev_yz_loc === duzdyl_plus_duydzl * 0.5
          print If(nspec_inner_core_strain_only == 1, lambda {
            print epsilon_trace_over_3[tx] === templ
          }, lambda {
            print epsilon_trace_over_3[tx + working_element*ngll3] === templ
          })
        }
        print kappal === d_kappav[offset]
        print mul === d_muv[offset]

        print If(attenuation, lambda {
          print If(use_3d_attenuation_arrays, lambda {
            print mul_iso  === mul * one_minus_sum_beta[tx+working_element*ngll3]
            print mul_aniso === mul * ( one_minus_sum_beta[tx+working_element*ngll3] - 1.0 )
          }, lambda {
            print mul_iso  === mul * one_minus_sum_beta[working_element]
            print mul_aniso === mul * ( one_minus_sum_beta[working_element] - 1.0 )
          })
        }, lambda {
          print mul_iso === mul
        })

        print If(anisotropy, lambda {
          print c11 === d_c11store[offset]
          print c12 === d_c12store[offset]
          print c13 === d_c13store[offset]
          print c33 === d_c33store[offset]
          print c44 === d_c44store[offset]
          print If(attenuation) {
            print c11 === c11 + mul_aniso * 1.33333333333333333333
            print c12 === c12 - mul_aniso * 0.66666666666666666666
            print c13 === c13 - mul_aniso * 0.66666666666666666666
            print c33 === c33 + mul_aniso * 1.33333333333333333333
            print c44 === c44 + mul_aniso
          }
          print sigma[0][0] === c11*dudl[0][0] + c12*dudl[1][1] + c13*dudl[2][2]
          print sigma[1][1] === c12*dudl[0][0] + c11*dudl[1][1] + c13*dudl[2][2]
          print sigma[2][2] === c13*dudl[0][0] + c13*dudl[1][1] + c33*dudl[2][2]
          print sigma[0][1] === (c11-c12)*duxdyl_plus_duydxl*0.5
          print sigma[0][2] === c44*duzdxl_plus_duxdzl
          print sigma[1][2] === c44*duzdyl_plus_duydzl
        }, lambda {
          print lambdalplus2mul === kappal + mul_iso * 1.33333333333333333333
          print lambdal === lambdalplus2mul - mul_iso * 2.0

          print sigma[0][0] === lambdalplus2mul*dudl[0][0] + lambdal*duydyl_plus_duzdzl
          print sigma[1][1] === lambdalplus2mul*dudl[1][1] + lambdal*duxdxl_plus_duzdzl
          print sigma[2][2] === lambdalplus2mul*dudl[2][2] + lambdal*duxdxl_plus_duydyl

          print sigma[0][1] === mul*duxdyl_plus_duydxl
          print sigma[0][2] === mul*duzdxl_plus_duxdzl
          print sigma[1][2] === mul*duzdyl_plus_duydzl
        })

        print If(Expression("&&", attenuation, !partial_phys_dispersion_only)) {
          print sub_kernel1.call(tx, working_element,\
                                 r_xx, r_yy, r_xy, r_xz, r_yz,\
                                 sigma[0][0].address, sigma[1][1].address, sigma[2][2].address,\
                                 sigma[0][1].address, sigma[0][2].address, sigma[1][2].address)
        }
        print sigma[1][0] === sigma[0][1]
        print sigma[2][0] === sigma[0][2]
        print sigma[2][1] === sigma[1][2]

        print jacobianl === Expression("/", 1.0, xil[0]*(etal[1]*gammal[2] - etal[2]*gammal[1])\
                                               - xil[1]*(etal[0]*gammal[2] - etal[2]*gammal[0])\
                                               + xil[2]*(etal[0]*gammal[1] - etal[1]*gammal[0]))
        print If( gravity ) {
          print sub_kernel2.call(tx, working_element,\
                                 d_ibool, d_store[0], d_store[1], d_store[2],\
                                 d_minus_gravity_table, d_minus_deriv_gravity_table, d_density_table,\
                                 wgll_cube, jacobianl,\
                                 s_dummy_loc[0], s_dummy_loc[1], s_dummy_loc[2],\
                                 sigma[0][0].address, sigma[1][1].address, sigma[2][2].address,\
                                 sigma[0][1].address, sigma[1][0].address, sigma[0][2].address,\
                                 sigma[2][0].address, sigma[1][2].address, sigma[2][1].address,\
                                 rho_s_H[0].address, rho_s_H[1].address, rho_s_H[2].address)
        }
        (0..2).each { |indx|
          print s_temp[indx][0] === jacobianl * (sigma[0][indx]*xil[0] + sigma[1][indx]*xil[1] + sigma[2][indx]*xil[2])
        }
        (0..2).each { |indx|
          print s_temp[indx][1] === jacobianl * (sigma[0][indx]*etal[0] + sigma[1][indx]*etal[1] + sigma[2][indx]*etal[2])
        }
        (0..2).each { |indx|
          print s_temp[indx][2] === jacobianl * (sigma[0][indx]*gammal[0] + sigma[1][indx]*gammal[1] + sigma[2][indx]*gammal[2])
        }
      }
      print barrier(:local)
      print If( active ){
        (0..2).each { |indx1|
          (0..2).each { |indx2|
            print tempanl[indx1][indx2] === 0.0
          }
        }
        for_loop = For(l, 0, ngllx-1) {
          print fac[0] === sh_hprimewgll_xx[i*ngllx + l]
          print offset === k*ngll2 + j*ngllx + l
          (0..2).each { |indx1|
            print tempanl[indx1][0] === tempanl[indx1][0] + s_temp[indx1][0][offset]*fac[0]
          }
          print fac[1] === sh_hprimewgll_xx[j*ngllx + l]
          print offset === k*ngll2 + l*ngllx + i
          (0..2).each { |indx1|
            print tempanl[indx1][1] === tempanl[indx1][1] + s_temp[indx1][1][offset]*fac[1]
          }
          print fac[2] === sh_hprimewgll_xx[k*ngllx + l]
          print offset === l*ngll2 + j*ngllx + i
          (0..2).each { |indx1|
            print tempanl[indx1][2] === tempanl[indx1][2] + s_temp[indx1][2][offset]*fac[2]
          }
        }
        if unroll_loops then
          for_loop.unroll
        else
          print for_loop
        end
        print fac[0] === d_wgllwgll_yz[k*ngllx+j]
        print fac[1] === d_wgllwgll_xz[k*ngllx+i]
        print fac[2] === d_wgllwgll_xy[j*ngllx+i]
        (0..2).each { |indx|
          print sum_terms[indx] === -(fac[0]*tempanl[indx][0] + fac[1]*tempanl[indx][1] + fac[2]*tempanl[indx][2])
        }

        print If( gravity ) {
          (0..2).each { |indx|
            print sum_terms[indx] === sum_terms[indx] + rho_s_H[indx]
          }
        }
        if mesh_coloring then
          if textures_fields then
            (0..2).each { |indx|
              print d_accel[indx,iglob] === d_accel_ic_tex[iglob*3+indx] + sum_terms[indx]
            }
          else
            (0..2).each { |indx|
              print d_accel[indx,iglob] === d_accel[indx,iglob] + sum_terms[indx]
            }
          end
        else
          print If( use_mesh_coloring_gpu, lambda {
            print If( nspec_inner_core > coloring_min_nspec_inner_core, lambda {
              if textures_fields then
                (0..2).each { |indx|
                  print d_accel[indx,iglob] === d_accel_ic_tex[iglob*3+indx] + sum_terms[indx]
                }
              else
                (0..2).each { |indx|
                  print d_accel[indx,iglob] === d_accel[indx,iglob] + sum_terms[indx]
                }
              end
            }, lambda{
              (0..2).each { |indx|
                print atomicAdd(d_accel+ iglob*3 +indx, sum_terms[indx])
              }
            })
          }, lambda {
            (0..2).each { |indx|
              print atomicAdd(d_accel + iglob*3 + indx, sum_terms[indx])
            }
          })
        end
        print If( Expression("&&", attenuation, !partial_phys_dispersion_only ) ) {
          print sub_kernel3.call( tx, working_element,\
                                  d_muv, factor_common,\
                                  alphaval, betaval, gammaval,\
                                  r_xx, r_yy, r_xy, r_xz, r_yz,\
                                  epsilondev_xx, epsilondev_yy, epsilondev_xy, epsilondev_xz, epsilondev_yz,\
                                  epsilondev_xx_loc, epsilondev_yy_loc, epsilondev_xy_loc, epsilondev_xz_loc, epsilondev_yz_loc,\
                                  use_3d_attenuation_arrays)
        }
        print If( compute_and_store_strain ) {
          print epsilondev_xx[tx + working_element*ngll3] === epsilondev_xx_loc
          print epsilondev_yy[tx + working_element*ngll3] === epsilondev_yy_loc
          print epsilondev_xy[tx + working_element*ngll3] === epsilondev_xy_loc
          print epsilondev_xz[tx + working_element*ngll3] === epsilondev_xz_loc
          print epsilondev_yz[tx + working_element*ngll3] === epsilondev_yz_loc
        }
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
