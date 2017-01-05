require './KSplitOssRef_v1.rb'
require './KSplitOssRef_v2.rb'
require './KSplitOssBoast.rb'
require 'narray_ffi'
require '../Common/CommonArgs.rb'

class Debug 
  include CommonArgs
  def self.run
    nests = [1,2,3,4,5,6,7,8,9,10,11]
    vector_size=2
    seed = 10

    pgaus = 8
    pnode = 8
    ndime = 3
    kernels = []

    k_orig_params = {:vector_length => vector_size, :preprocessor => false, :nests => nests}
    k_boast_params = {:vector_length => vector_size, :nests => nests, :unroll => false, :inline => :inlined}

    set_fortran_line_length(100)
    kernels[0] = KSplitOssRef_v2::new(k_orig_params)
    kernels[0].generate
    kernels[0].kernel.build(:FCFLAGS => "-cpp", :LDFLAGS => "-lgfortran")

    set_lang(C)
    kernels[1] = KSplitBoast::new(k_boast_params,ndime)
    kernels[1].generate
    puts kernels[1].kernel
    
    stats_boast = []
    stats_ref = []

    # # # 100.times{|i|
    i = 1
    CommonArgs.init(pgaus,pnode,ndime,vector_size,2,seed)
    @@kfl_lumped = 1 # 2
    @@kfl_limit_nsi = 1 # 2
    @@kfl_stabi_nsi = 1 # -1

    kernels[0].kernel.run(@@ndime,@@kfl_lumped,@@mnode,@@ntens,@@kfl_duatss,@@fact_duatss,@@kfl_stabi_nsi,
                       @@fvins_nsi,@@fcons_nsi,@@bemol_nsi,@@kfl_regim_nsi,@@fvela_nsi,@@kfl_rmom2_nsi,
                       @@kfl_press_nsi,@@kfl_p1ve2_nsi,@@kfl_linea_nsi,@@kfl_confi_nsi,@@nbdfp_nsi,
                       @@kfl_sgsti_nsi,@@kfl_nota1_nsi,@@kfl_limit_nsi,@@kfl_penal_nsi,@@penal_nsi,
                       @@kfl_bubbl_nsi,@@pnode,@@pgaus,@@gpden,@@gpvis,@@gppor,@@gpsp1,@@gpsp2,@@gpvol,
                       @@gpsha,@@gpcar,@@gpadv,@@gpvep[0],@@gpgrp[0],@@gprhs[0],@@gprhc[0],@@gpvel,
                       @@gpsgs,@@elvel,@@elpre,@@elbub,@@wgrgr[0],@@agrau[0],@@elauu[0],@@elaup[0],
                       @@elapp[0],@@elapu[0],@@elrbu[0],@@elrbp[0],@@dtinv_loc,@@dtsgs,@@pbubl,
                       @@gpsha_bub,@@gpcar_bub,@@elauq[0],@@elapq[0],@@elaqu[0],@@elaqp[0],
                       @@elaqq[0],@@elrbq[0])

    kernels[1].kernel.run(@@kfl_lumped,@@mnode,@@ntens,@@kfl_duatss,@@fact_duatss,@@kfl_stabi_nsi,
                       @@fvins_nsi,@@fcons_nsi,@@bemol_nsi,@@kfl_regim_nsi,@@fvela_nsi,@@kfl_rmom2_nsi,
                       @@kfl_press_nsi,@@kfl_p1ve2_nsi,@@kfl_linea_nsi,@@kfl_confi_nsi,@@nbdfp_nsi,
                       @@kfl_sgsti_nsi,@@kfl_nota1_nsi,@@kfl_limit_nsi,@@kfl_penal_nsi,@@penal_nsi,
                       @@kfl_bubbl_nsi,@@pnode,@@pgaus,@@gpden,@@gpvis,@@gppor,@@gpsp1,@@gpsp2,@@gpvol,
                       @@gpsha,@@gpcar,@@gpadv,@@gpvep[1],@@gpgrp[1],@@gprhs[1],@@gprhc[1],@@gpvel,
                       @@gpsgs,@@elvel,@@elpre,@@elbub,@@wgrgr[1],@@agrau[1],@@elauu[1],@@elaup[1],
                       @@elapp[1],@@elapu[1],@@elrbu[1],@@elrbp[1],@@dtinv_loc,@@dtsgs,@@pbubl,
                       @@gpsha_bub,@@gpcar_bub,@@elauq[1],@@elapq[1],@@elaqu[1],@@elaqp[1],
                       @@elaqq[1],@@elrbq[1])

    epsilon = 10e-15

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

    raise "Error: residue too big for agrau" if (diff_agrau > epsilon).to_a.flatten.include? 1
    raise "Error: residue too big for wgrgr" if (diff_wgrgr > epsilon).to_a.flatten.include? 1
    raise "Error: residue too big for elauu" if (diff_elauu > epsilon).to_a.flatten.include? 1
    raise "Error: residue too big for elrbu" if (diff_elrbu > epsilon).to_a.flatten.include? 1
    raise "Error: residue too big for elapu" if (diff_elapu > epsilon).to_a.flatten.include? 1
    raise "Error: residue too big for elaqu" if (diff_elaqu > epsilon).to_a.flatten.include? 1
    raise "Error: residue too big for elaqp" if (diff_elaqp > epsilon).to_a.flatten.include? 1
    raise "Error: residue too big for elaqq" if (diff_elaqq > epsilon).to_a.flatten.include? 1
    raise "Error: residue too big for elapq" if (diff_elapq > epsilon).to_a.flatten.include? 1
    raise "Error: residue too big for elauq" if (diff_elauq > epsilon).to_a.flatten.include? 1
    raise "Error: residue too big for elaup" if (diff_elaup > epsilon).to_a.flatten.include? 1
    raise "Error: residue too big for elapp" if (diff_elapp > epsilon).to_a.flatten.include? 1
    raise "Error: residue too big for elrbp" if (diff_elrbp > epsilon).to_a.flatten.include? 1
    raise "Error: residue too big for elrbq" if (diff_elrbq > epsilon).to_a.flatten.include? 1
    raise "Error: residue too big for gpgrp" if (diff_gpgrp > epsilon).to_a.flatten.include? 1
    raise "Error: residue too big for gprhs" if (diff_gprhs > epsilon).to_a.flatten.include? 1
    raise "Error: residue too big for gprhc" if (diff_gprhc > epsilon).to_a.flatten.include? 1
    raise "Error: residue too big for gpvep" if (diff_gpvep > epsilon).to_a.flatten.include? 1

    puts "Done"

    # }
    # t = []
    # t[1] = stats_ref.inject(0){|sum,e| sum+e[:duration]} / stats_ref.length
    # t[2] = stats_boast.inject(0){|sum,e| sum+e[:duration]} / stats_boast.length
    #  t
    # wgrgr_ref
  end
  def save_kernel(path, kernel)
    file_boast = File::new(path,"w")
    if file_boast then
      file_boast.syswrite(kernel)
    else
      puts "Failed to open #{path}"
    end
  end
end

Debug.run
