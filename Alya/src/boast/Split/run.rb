require_relative './KSplitOssRef_v1.rb'
require_relative './KSplitOssRef_v2.rb'
require_relative './KSplitOssBoast.rb'
require 'narray_ffi'
require_relative '../Common/CommonArgs.rb'
require 'yaml'
require 'pp'
require 'csv'
require 'optparse'
require_relative '../Tools/LogInfo.rb'

class Experiment
  include CommonArgs
  def self.run(output_info,dimension)
    nests = (1..11).to_a
    LogInfo.init(output_info)
    LogInfo.get_info
    vector_size=2
    seed = 10
    charac = {:pgaus => 8, :pnode => 8}
    epsilon = 10e-15

    stats = {}
    k = {}
    nests.each{|n|
      opts1 = {:vector_length => vector_size, :preprocessor => false, :nests => [n], :unroll => true, :inline => :included, :CFLAGS => "-O3"}
      opts2 = {:vector_length => vector_size, :preprocessor => false, :nests => [n], :unroll => true, :inline => :inlined, :CFLAGS => "-O3"}
      opts3 = {:vector_length => vector_size, :preprocessor => false, :nests => [n], :unroll => true, :inline => :call, :CFLAGS => "-O3"}

      set_lang(C)
      set_fortran_line_length(100)
      # k[h_ref] = KSplitOssRef_v2::new(opts)
      k[opts1] = KSplitBoast::new(opts1,dimension)
      k[opts1].generate
      LogInfo.register_kernel_info(opts1, k[opts1].kernel.to_s)
      # k[opts1].kernel.build(:FCFLAGS => "-O3")
      k[opts1].kernel.build(:CFLAGS => opts1[:CFLAGS])

      set_lang(C)
      k[opts2] = KSplitBoast::new(opts2,dimension)
      k[opts2].generate
      LogInfo.register_kernel_info(opts2, k[opts2].kernel.to_s)
      k[opts2].kernel.build(:CFLAGS => opts2[:CFLAGS])

      k[opts3] = KSplitBoast::new(opts3,dimension)
      k[opts3].generate
      LogInfo.register_kernel_info(opts3, k[opts3].kernel.to_s)
      k[opts3].kernel.build(:CFLAGS => opts3[:CFLAGS])
      
      stats[opts1] = {:time => []}
      stats[opts2] = {:time => []}
      stats[opts3] = {:time => []}
      stats[opts1][:characteristics] = charac
      stats[opts2][:characteristics] = charac
      stats[opts3][:characteristics] = charac

      CommonArgs.init(charac[:pgaus],charac[:pnode],dimension,vector_size,3,seed)
      @@kfl_lumped = 2 # 1
      @@kfl_limit_nsi = 1 # 2
      @@kfl_stabi_nsi = 1 # -1
      100.times{|i|
        k[opts1].kernel.run(@@kfl_lumped,@@mnode,@@ntens,@@kfl_duatss,@@fact_duatss,@@kfl_stabi_nsi,
               @@fvins_nsi,@@fcons_nsi,@@bemol_nsi,@@kfl_regim_nsi,@@fvela_nsi,@@kfl_rmom2_nsi,
               @@kfl_press_nsi,@@kfl_p1ve2_nsi,@@kfl_linea_nsi,@@kfl_confi_nsi,@@nbdfp_nsi,
               @@kfl_sgsti_nsi,@@kfl_nota1_nsi,@@kfl_limit_nsi,@@kfl_penal_nsi,@@penal_nsi,
               @@kfl_bubbl_nsi,@@pnode,@@pgaus,@@gpden,@@gpvis,@@gppor,@@gpsp1,@@gpsp2,@@gpvol,
               @@gpsha,@@gpcar,@@gpadv,@@gpvep[0],@@gpgrp[0],@@gprhs[0],@@gprhc[0],@@gpvel,
               @@gpsgs,@@elvel,@@elpre,@@elbub,@@wgrgr[0],@@agrau[0],@@elauu[0],@@elaup[0],
               @@elapp[0],@@elapu[0],@@elrbu[0],@@elrbp[0],@@dtinv_loc,@@dtsgs,@@pbubl,
               @@gpsha_bub,@@gpcar_bub,@@elauq[0],@@elapq[0],@@elaqu[0],@@elaqp[0],
               @@elaqq[0],@@elrbq[0])

        k[opts2].kernel.run(@@kfl_lumped,@@mnode,@@ntens,@@kfl_duatss,@@fact_duatss,@@kfl_stabi_nsi,
               @@fvins_nsi,@@fcons_nsi,@@bemol_nsi,@@kfl_regim_nsi,@@fvela_nsi,@@kfl_rmom2_nsi,
               @@kfl_press_nsi,@@kfl_p1ve2_nsi,@@kfl_linea_nsi,@@kfl_confi_nsi,@@nbdfp_nsi,
               @@kfl_sgsti_nsi,@@kfl_nota1_nsi,@@kfl_limit_nsi,@@kfl_penal_nsi,@@penal_nsi,
               @@kfl_bubbl_nsi,@@pnode,@@pgaus,@@gpden,@@gpvis,@@gppor,@@gpsp1,@@gpsp2,@@gpvol,
               @@gpsha,@@gpcar,@@gpadv,@@gpvep[1],@@gpgrp[1],@@gprhs[1],@@gprhc[1],@@gpvel,
               @@gpsgs,@@elvel,@@elpre,@@elbub,@@wgrgr[1],@@agrau[1],@@elauu[1],@@elaup[1],
               @@elapp[1],@@elapu[1],@@elrbu[1],@@elrbp[1],@@dtinv_loc,@@dtsgs,@@pbubl,
               @@gpsha_bub,@@gpcar_bub,@@elauq[1],@@elapq[1],@@elaqu[1],@@elaqp[1],
               @@elaqq[1],@@elrbq[1])

        k[opts3].kernel.run(@@kfl_lumped,@@mnode,@@ntens,@@kfl_duatss,@@fact_duatss,@@kfl_stabi_nsi,
               @@fvins_nsi,@@fcons_nsi,@@bemol_nsi,@@kfl_regim_nsi,@@fvela_nsi,@@kfl_rmom2_nsi,
               @@kfl_press_nsi,@@kfl_p1ve2_nsi,@@kfl_linea_nsi,@@kfl_confi_nsi,@@nbdfp_nsi,
               @@kfl_sgsti_nsi,@@kfl_nota1_nsi,@@kfl_limit_nsi,@@kfl_penal_nsi,@@penal_nsi,
               @@kfl_bubbl_nsi,@@pnode,@@pgaus,@@gpden,@@gpvis,@@gppor,@@gpsp1,@@gpsp2,@@gpvol,
               @@gpsha,@@gpcar,@@gpadv,@@gpvep[2],@@gpgrp[2],@@gprhs[2],@@gprhc[2],@@gpvel,
               @@gpsgs,@@elvel,@@elpre,@@elbub,@@wgrgr[2],@@agrau[2],@@elauu[2],@@elaup[2],
               @@elapp[2],@@elapu[2],@@elrbu[2],@@elrbp[2],@@dtinv_loc,@@dtsgs,@@pbubl,
               @@gpsha_bub,@@gpcar_bub,@@elauq[2],@@elapq[2],@@elaqu[2],@@elaqp[2],
               @@elaqq[2],@@elrbq[2])

      }
      
      1000.times{|i|
        CommonArgs.init(charac[:pgaus],charac[:pnode],dimension,vector_size,3,seed)
        @@kfl_lumped = 2 # 1
        @@kfl_limit_nsi = 1 # 2
        @@kfl_stabi_nsi = 1 # -1
        
        stats[opts1][:time][i] = k[opts1].kernel.run(@@kfl_lumped,@@mnode,@@ntens,@@kfl_duatss,@@fact_duatss,@@kfl_stabi_nsi,
               @@fvins_nsi,@@fcons_nsi,@@bemol_nsi,@@kfl_regim_nsi,@@fvela_nsi,@@kfl_rmom2_nsi,
               @@kfl_press_nsi,@@kfl_p1ve2_nsi,@@kfl_linea_nsi,@@kfl_confi_nsi,@@nbdfp_nsi,
               @@kfl_sgsti_nsi,@@kfl_nota1_nsi,@@kfl_limit_nsi,@@kfl_penal_nsi,@@penal_nsi,
               @@kfl_bubbl_nsi,@@pnode,@@pgaus,@@gpden,@@gpvis,@@gppor,@@gpsp1,@@gpsp2,@@gpvol,
               @@gpsha,@@gpcar,@@gpadv,@@gpvep[0],@@gpgrp[0],@@gprhs[0],@@gprhc[0],@@gpvel,
               @@gpsgs,@@elvel,@@elpre,@@elbub,@@wgrgr[0],@@agrau[0],@@elauu[0],@@elaup[0],
               @@elapp[0],@@elapu[0],@@elrbu[0],@@elrbp[0],@@dtinv_loc,@@dtsgs,@@pbubl,
               @@gpsha_bub,@@gpcar_bub,@@elauq[0],@@elapq[0],@@elaqu[0],@@elaqp[0],
               @@elaqq[0],@@elrbq[0])[:duration]

        stats[opts2][:time][i] = k[opts2].kernel.run(@@kfl_lumped,@@mnode,@@ntens,@@kfl_duatss,@@fact_duatss,@@kfl_stabi_nsi,
               @@fvins_nsi,@@fcons_nsi,@@bemol_nsi,@@kfl_regim_nsi,@@fvela_nsi,@@kfl_rmom2_nsi,
               @@kfl_press_nsi,@@kfl_p1ve2_nsi,@@kfl_linea_nsi,@@kfl_confi_nsi,@@nbdfp_nsi,
               @@kfl_sgsti_nsi,@@kfl_nota1_nsi,@@kfl_limit_nsi,@@kfl_penal_nsi,@@penal_nsi,
               @@kfl_bubbl_nsi,@@pnode,@@pgaus,@@gpden,@@gpvis,@@gppor,@@gpsp1,@@gpsp2,@@gpvol,
               @@gpsha,@@gpcar,@@gpadv,@@gpvep[1],@@gpgrp[1],@@gprhs[1],@@gprhc[1],@@gpvel,
               @@gpsgs,@@elvel,@@elpre,@@elbub,@@wgrgr[1],@@agrau[1],@@elauu[1],@@elaup[1],
               @@elapp[1],@@elapu[1],@@elrbu[1],@@elrbp[1],@@dtinv_loc,@@dtsgs,@@pbubl,
               @@gpsha_bub,@@gpcar_bub,@@elauq[1],@@elapq[1],@@elaqu[1],@@elaqp[1],
               @@elaqq[1],@@elrbq[1])[:duration]

        stats[opts3][:time][i] = k[opts3].kernel.run(@@kfl_lumped,@@mnode,@@ntens,@@kfl_duatss,@@fact_duatss,@@kfl_stabi_nsi,
               @@fvins_nsi,@@fcons_nsi,@@bemol_nsi,@@kfl_regim_nsi,@@fvela_nsi,@@kfl_rmom2_nsi,
               @@kfl_press_nsi,@@kfl_p1ve2_nsi,@@kfl_linea_nsi,@@kfl_confi_nsi,@@nbdfp_nsi,
               @@kfl_sgsti_nsi,@@kfl_nota1_nsi,@@kfl_limit_nsi,@@kfl_penal_nsi,@@penal_nsi,
               @@kfl_bubbl_nsi,@@pnode,@@pgaus,@@gpden,@@gpvis,@@gppor,@@gpsp1,@@gpsp2,@@gpvol,
               @@gpsha,@@gpcar,@@gpadv,@@gpvep[2],@@gpgrp[2],@@gprhs[2],@@gprhc[2],@@gpvel,
               @@gpsgs,@@elvel,@@elpre,@@elbub,@@wgrgr[2],@@agrau[2],@@elauu[2],@@elaup[2],
               @@elapp[2],@@elapu[2],@@elrbu[2],@@elrbp[2],@@dtinv_loc,@@dtsgs,@@pbubl,
               @@gpsha_bub,@@gpcar_bub,@@elauq[2],@@elapq[2],@@elaqu[2],@@elaqp[2],
               @@elaqq[2],@@elrbq[2])[:duration]

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
      }
    }
    LogInfo.dump_info
    puts "Done"
    return stats,k
  end
end

options = {}

opt_parser = OptionParser.new { |opts|
  opts.banner = "Usage: run.rb --[options]=[value]"

  opts.on("-dVAL", "--data=VAL", "Specify the path where to store the data in yaml format") { |n|
    options[:data_path] = n
  }

  opts.on("-iVAL", "--info=VAL", "Specify the path where to store the information such as the generated source, compiler info, etc... in yaml format") { |n|
    options[:info_path] = n
  }

  opts.on("--dimension=VAL", "Specify the number of dimension") { |n|
    options[:dimension] = n.to_i
  }

  opts.on("-h", "--help", "Prints this help") {
    puts opts
    exit
  }
}.parse!

stats, infos = Experiment.run(options[:info_path],options[:dimension])

if options[:data_path] then
  File::open( options[:data_path], "w") { |f|
    f.print YAML::dump(stats)
  }
end
