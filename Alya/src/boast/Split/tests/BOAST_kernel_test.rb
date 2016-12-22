gem "minitest"
require 'minitest/autorun'
require 'BOAST'
include BOAST
require '../KSplitOssRef_v1.rb'
require '../KSplitOssRef_v2.rb'
require '../KSplitOssBoast.rb'
require '../../Common/CommonArgs.rb'

class TestFor < Minitest::Test
  include CommonArgs
  def test_nests_orig_vs_boast_fortran
    epsilon = 10e-15
    set_fortran_line_length(100)
    nests = [1,2,3,4,5,6,7,8,9,10,11]

    pgaus = 4
    pnode = 1
    ndime = 3

    kernels = {}
    
    (1..2).each{|vector_size|
      k_ref_params = {:vector_length => vector_size, :preprocessor => false, :nests => nests}
      k_boast_params = {:vector_length => vector_size, :nests => nests}

      kernels[:ref] = KSplitOssRef_v2::new(k_ref_params)
      kernels[:ref].generate
      kernels[:ref].kernel.build(:FCFLAGS => "-cpp")
      
      [:included,:call,:inlined].each{|inlining|
        k_boast_params[:inline] = inlining
        [false,true].each{|unroll|
          k_boast_params[:unroll] = unroll
          (2..3).each{|dim|
            CommonArgs.init(pgaus,pnode,dim,vector_size,2)

            (1..2).each{|kfl_lumped|
              next if kfl_lumped == 1 and dim == 2
              @@kfl_lumped = kfl_lumped

              [-1,1].each{|kfl_stabi_nsi|
                @@kfl_stabi_nsi = kfl_stabi_nsi

                (1..2).each{|kfl_limit_nsi|
                  @@kfl_limit_nsi = kfl_limit_nsi
                  [FORTRAN].each{|lang|
                    set_lang(lang)
                    kernels[:boast] = KSplitBoast::new(k_boast_params,dim)
                    kernels[:boast].generate
                    kernels[:boast].kernel.build
                    kernels[:ref].kernel.run(dim,@@kfl_lumped,@@mnode,@@ntens,@@kfl_duatss,@@fact_duatss,@@kfl_stabi_nsi,
                                               @@fvins_nsi,@@fcons_nsi,@@bemol_nsi,@@kfl_regim_nsi,@@fvela_nsi,@@kfl_rmom2_nsi,
                                               @@kfl_press_nsi,@@kfl_p1ve2_nsi,@@kfl_linea_nsi,@@kfl_confi_nsi,@@nbdfp_nsi,
                                               @@kfl_sgsti_nsi,@@kfl_nota1_nsi,@@kfl_limit_nsi,@@kfl_penal_nsi,@@penal_nsi,
                                               @@kfl_bubbl_nsi,@@pnode,@@pgaus,@@gpden,@@gpvis,@@gppor,@@gpsp1,@@gpsp2,@@gpvol,
                                               @@gpsha,@@gpcar,@@gpadv,@@gpvep[0],@@gpgrp[0],@@gprhs[0],@@gprhc[0],@@gpvel,
                                               @@gpsgs,@@elvel,@@elpre,@@elbub,@@wgrgr[0],@@agrau[0],@@elauu[0],@@elaup[0],
                                               @@elapp[0],@@elapu[0],@@elrbu[0],@@elrbp[0],@@dtinv_loc,@@dtsgs,@@pbubl,
                                               @@gpsha_bub,@@gpcar_bub,@@elauq[0],@@elapq[0],@@elaqu[0],@@elaqp[0],
                                               @@elaqq[0],@@elrbq[0])

                    kernels[:boast].kernel.run(@@kfl_lumped,@@mnode,@@ntens,@@kfl_duatss,@@fact_duatss,@@kfl_stabi_nsi,
                                               @@fvins_nsi,@@fcons_nsi,@@bemol_nsi,@@kfl_regim_nsi,@@fvela_nsi,@@kfl_rmom2_nsi,
                                               @@kfl_press_nsi,@@kfl_p1ve2_nsi,@@kfl_linea_nsi,@@kfl_confi_nsi,@@nbdfp_nsi,
                                               @@kfl_sgsti_nsi,@@kfl_nota1_nsi,@@kfl_limit_nsi,@@kfl_penal_nsi,@@penal_nsi,
                                               @@kfl_bubbl_nsi,@@pnode,@@pgaus,@@gpden,@@gpvis,@@gppor,@@gpsp1,@@gpsp2,@@gpvol,
                                               @@gpsha,@@gpcar,@@gpadv,@@gpvep[1],@@gpgrp[1],@@gprhs[1],@@gprhc[1],@@gpvel,
                                               @@gpsgs,@@elvel,@@elpre,@@elbub,@@wgrgr[1],@@agrau[1],@@elauu[1],@@elaup[1],
                                               @@elapp[1],@@elapu[1],@@elrbu[1],@@elrbp[1],@@dtinv_loc,@@dtsgs,@@pbubl,
                                               @@gpsha_bub,@@gpcar_bub,@@elauq[1],@@elapq[1],@@elaqu[1],@@elaqp[1],
                                               @@elaqq[1],@@elrbq[1])

                    diff_agrau = (@@agrau[0] - @@agrau[1]).abs
                    diff_wgrgr = (@@wgrgr[0] - @@wgrgr[1]).abs
                    diff_elauu = (@@elauu[0] - @@elauu[1]).abs
                    diff_elrbu = (@@elrbu[0] - @@elrbu[1]).abs
                    diff_elapu = (@@elapu[0] - @@elapu[1]).abs
                    diff_elaqu = (@@elaqu[0] - @@elaqu[1]).abs
                    diff_elaqp = (@@elaqp[0] - @@elaqp[1]).abs
                    diff_elaqq = (@@elaqq[0] - @@elaqq[1]).abs
                    diff_elapq = (@@elapq[0] - @@elapq[1]).abs
                    diff_elauq = (@@elauq[0] - @@elauq[1]).abs
                    diff_elaup = (@@elaup[0] - @@elaup[1]).abs
                    diff_elapp = (@@elapp[0] - @@elapp[1]).abs
                    diff_elrbp = (@@elrbp[0] - @@elrbp[1]).abs
                    diff_elrbq = (@@elrbq[0] - @@elrbq[1]).abs

                    diff_gpgrp = (@@gpgrp[1] - @@gpgrp[0]).abs
                    diff_gprhs = (@@gprhs[1] - @@gprhs[0]).abs
                    diff_gpvep = (@@gpvep[1] - @@gpvep[0]).abs
                    diff_gprhc = (@@gprhc[1] - @@gprhc[0]).abs
                    
                    
                    assert( (not (diff_agrau > epsilon).to_a.flatten.include? 1))
                    assert( (not (diff_wgrgr > epsilon).to_a.flatten.include? 1))
                    assert( (not (diff_elauu > epsilon).to_a.flatten.include? 1))
                    assert( (not (diff_elrbu > epsilon).to_a.flatten.include? 1))
                    assert( (not (diff_elapu > epsilon).to_a.flatten.include? 1))
                    assert( (not (diff_elaqu > epsilon).to_a.flatten.include? 1))
                    assert( (not (diff_elaqp > epsilon).to_a.flatten.include? 1))
                    assert( (not (diff_elaqq > epsilon).to_a.flatten.include? 1))
                    assert( (not (diff_elapq > epsilon).to_a.flatten.include? 1))
                    assert( (not (diff_elauq > epsilon).to_a.flatten.include? 1))
                    assert( (not (diff_elaup > epsilon).to_a.flatten.include? 1))
                    assert( (not (diff_elapp > epsilon).to_a.flatten.include? 1))
                    assert( (not (diff_elrbp > epsilon).to_a.flatten.include? 1))
                    assert( (not (diff_elrbq > epsilon).to_a.flatten.include? 1))
                    assert( (not (diff_gpgrp > epsilon).to_a.flatten.include? 1))
                    assert( (not (diff_gprhs > epsilon).to_a.flatten.include? 1))
                    assert( (not (diff_gprhc > epsilon).to_a.flatten.include? 1))
                    assert( (not (diff_gpvep > epsilon).to_a.flatten.include? 1))
                  }
                }
              }
            }
          }
        }
      }
    }
  end

  def test_nests_orig_vs_boast_c
    epsilon = 10e-15
    set_fortran_line_length(100)
    nests = [1,2,3,4,5,6,7,8,9,10,11]

    pgaus = 8
    pnode = 8
    ndime = 3

    kernels = {}

    (1..2).each{|vector_size|
      k_ref_params = {:vector_length => vector_size, :preprocessor => false, :nests => nests}
      k_boast_params = {:vector_length => vector_size, :nests => nests}

      kernels[:ref] = KSplitOssRef_v2::new(k_ref_params)
      kernels[:ref].generate
      kernels[:ref].kernel.build(:FCFLAGS => "-cpp")

      [:included,:call,:inlined].each{|inlining|
        k_boast_params[:inline] = inlining
        [false,true].each{|unroll|
          k_boast_params[:unroll] = unroll
          (2..3).each{|dim|
            CommonArgs.init(pgaus,pnode,dim,vector_size,2)

            (1..2).each{|kfl_lumped|
              next if kfl_lumped == 1 and dim == 2
              @@kfl_lumped = kfl_lumped

              [-1,1].each{|kfl_stabi_nsi|
                @@kfl_stabi_nsi = kfl_stabi_nsi

                (1..2).each{|kfl_limit_nsi|
                  @@kfl_limit_nsi = kfl_limit_nsi
                  [C].each{|lang|
                    set_lang(lang)
                    kernels[:boast] = KSplitBoast::new(k_boast_params,dim)
                    kernels[:boast].generate
                    kernels[:boast].kernel.build
                    kernels[:ref].kernel.run(dim,@@kfl_lumped,@@mnode,@@ntens,@@kfl_duatss,@@fact_duatss,@@kfl_stabi_nsi,
                                 @@fvins_nsi,@@fcons_nsi,@@bemol_nsi,@@kfl_regim_nsi,@@fvela_nsi,@@kfl_rmom2_nsi,
                                 @@kfl_press_nsi,@@kfl_p1ve2_nsi,@@kfl_linea_nsi,@@kfl_confi_nsi,@@nbdfp_nsi,
                                 @@kfl_sgsti_nsi,@@kfl_nota1_nsi,@@kfl_limit_nsi,@@kfl_penal_nsi,@@penal_nsi,
                                 @@kfl_bubbl_nsi,@@pnode,@@pgaus,@@gpden,@@gpvis,@@gppor,@@gpsp1,@@gpsp2,@@gpvol,
                                 @@gpsha,@@gpcar,@@gpadv,@@gpvep[0],@@gpgrp[0],@@gprhs[0],@@gprhc[0],@@gpvel,
                                 @@gpsgs,@@elvel,@@elpre,@@elbub,@@wgrgr[0],@@agrau[0],@@elauu[0],@@elaup[0],
                                 @@elapp[0],@@elapu[0],@@elrbu[0],@@elrbp[0],@@dtinv_loc,@@dtsgs,@@pbubl,
                                 @@gpsha_bub,@@gpcar_bub,@@elauq[0],@@elapq[0],@@elaqu[0],@@elaqp[0],
                                 @@elaqq[0],@@elrbq[0])

                    kernels[:boast].kernel.run(@@kfl_lumped,@@mnode,@@ntens,@@kfl_duatss,@@fact_duatss,@@kfl_stabi_nsi,
                                 @@fvins_nsi,@@fcons_nsi,@@bemol_nsi,@@kfl_regim_nsi,@@fvela_nsi,@@kfl_rmom2_nsi,
                                 @@kfl_press_nsi,@@kfl_p1ve2_nsi,@@kfl_linea_nsi,@@kfl_confi_nsi,@@nbdfp_nsi,
                                 @@kfl_sgsti_nsi,@@kfl_nota1_nsi,@@kfl_limit_nsi,@@kfl_penal_nsi,@@penal_nsi,
                                 @@kfl_bubbl_nsi,@@pnode,@@pgaus,@@gpden,@@gpvis,@@gppor,@@gpsp1,@@gpsp2,@@gpvol,
                                 @@gpsha,@@gpcar,@@gpadv,@@gpvep[1],@@gpgrp[1],@@gprhs[1],@@gprhc[1],@@gpvel,
                                 @@gpsgs,@@elvel,@@elpre,@@elbub,@@wgrgr[1],@@agrau[1],@@elauu[1],@@elaup[1],
                                 @@elapp[1],@@elapu[1],@@elrbu[1],@@elrbp[1],@@dtinv_loc,@@dtsgs,@@pbubl,
                                 @@gpsha_bub,@@gpcar_bub,@@elauq[1],@@elapq[1],@@elaqu[1],@@elaqp[1],
                                 @@elaqq[1],@@elrbq[1])

                    diff_agrau = (@@agrau[0] - @@agrau[1]).abs
                    diff_wgrgr = (@@wgrgr[0] - @@wgrgr[1]).abs
                    diff_elauu = (@@elauu[0] - @@elauu[1]).abs
                    diff_elrbu = (@@elrbu[0] - @@elrbu[1]).abs
                    diff_elapu = (@@elapu[0] - @@elapu[1]).abs
                    diff_elaqu = (@@elaqu[0] - @@elaqu[1]).abs
                    diff_elaqp = (@@elaqp[0] - @@elaqp[1]).abs
                    diff_elaqq = (@@elaqq[0] - @@elaqq[1]).abs
                    diff_elapq = (@@elapq[0] - @@elapq[1]).abs
                    diff_elauq = (@@elauq[0] - @@elauq[1]).abs
                    diff_elaup = (@@elaup[0] - @@elaup[1]).abs
                    diff_elapp = (@@elapp[0] - @@elapp[1]).abs
                    diff_elrbp = (@@elrbp[0] - @@elrbp[1]).abs
                    diff_elrbq = (@@elrbq[0] - @@elrbq[1]).abs

                    diff_gpgrp = (@@gpgrp[1] - @@gpgrp[0]).abs
                    diff_gprhs = (@@gprhs[1] - @@gprhs[0]).abs
                    diff_gpvep = (@@gpvep[1] - @@gpvep[0]).abs
                    diff_gprhc = (@@gprhc[1] - @@gprhc[0]).abs
                    
                    
                    assert( (not (diff_agrau > epsilon).to_a.flatten.include? 1))
                    assert( (not (diff_wgrgr > epsilon).to_a.flatten.include? 1))
                    assert( (not (diff_elauu > epsilon).to_a.flatten.include? 1))
                    assert( (not (diff_elrbu > epsilon).to_a.flatten.include? 1))
                    assert( (not (diff_elapu > epsilon).to_a.flatten.include? 1))
                    assert( (not (diff_elaqu > epsilon).to_a.flatten.include? 1))
                    assert( (not (diff_elaqp > epsilon).to_a.flatten.include? 1))
                    assert( (not (diff_elaqq > epsilon).to_a.flatten.include? 1))
                    assert( (not (diff_elapq > epsilon).to_a.flatten.include? 1))
                    assert( (not (diff_elauq > epsilon).to_a.flatten.include? 1))
                    assert( (not (diff_elaup > epsilon).to_a.flatten.include? 1))
                    assert( (not (diff_elapp > epsilon).to_a.flatten.include? 1))
                    assert( (not (diff_elrbp > epsilon).to_a.flatten.include? 1))
                    assert( (not (diff_elrbq > epsilon).to_a.flatten.include? 1))
                    assert( (not (diff_gpgrp > epsilon).to_a.flatten.include? 1))
                    assert( (not (diff_gprhs > epsilon).to_a.flatten.include? 1))
                    assert( (not (diff_gprhc > epsilon).to_a.flatten.include? 1))
                    assert( (not (diff_gpvep > epsilon).to_a.flatten.include? 1))
                  }
                }
              }
            }
          }
        }
      }
    }
  end
end
