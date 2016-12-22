
require_relative '../Common/Parameters.rb'
class KSplitOss
  include Parameters
  attr_reader :kernel

  def declare_parameters(ndime = nil)
    Parameters.initialize(@opts[:vector_length],ndime)
    @args = []
    if ndime.nil? then
      @args.push @ndime       = $ndime
    else
      @ndime                  = ndime
    end
    @args.push @kfl_lumped    = $kfl_lumped
    @args.push @mnode         = $mnode
    @args.push @ntens         = $ntens
    @args.push @kfl_duatss    = $kfl_duatss
    @args.push @fact_duatss   = $fact_duatss
    @args.push @kfl_stabi_nsi = $kfl_stabi_nsi
    @args.push @fvins_nsi     = $fvins_nsi
    @args.push @fcons_nsi     = $fcons_nsi
    @args.push @bemol_nsi     = $bemol_nsi
    @args.push @kfl_regim_nsi = $kfl_regim_nsi
    @args.push @fvela_nsi     = $fvela_nsi
    @args.push @kfl_rmom2_nsi = $kfl_rmom2_nsi
    @args.push @kfl_press_nsi = $kfl_press_nsi
    @args.push @kfl_p1ve2_nsi = $kfl_p1ve2_nsi
    @args.push @kfl_linea_nsi = $kfl_linea_nsi
    @args.push @kfl_confi_nsi = $kfl_confi_nsi
    @args.push @nbdfp_nsi     = $nbdfp_nsi
    @args.push @kfl_sgsti_nsi = $kfl_sgsti_nsi
    @args.push @kfl_nota1_nsi = $kfl_nota1_nsi
    @args.push @kfl_limit_nsi = $kfl_limit_nsi
    @args.push @kfl_penal_nsi = $kfl_penal_nsi
    @args.push @penal_nsi     = $penal_nsi
    @args.push @kfl_bubbl_nsi = $kfl_bubbl_nsi
  
    @args.push @pnode     = $pnode
    @args.push @pgaus     = $pgaus
    @args.push @gpden     = $gpden
    @args.push @gpvis     = $gpvis
    @args.push @gppor     = $gppor
    @args.push @gpsp1     = $gpsp1
    @args.push @gpsp2     = $gpsp2
    @args.push @gpvol     = $gpvol
    @args.push @gpsha     = $gpsha
    @args.push @gpcar     = $gpcar
    @args.push @gpadv     = $gpadv
    @args.push @gpvep     = $gpvep
    @args.push @gpgrp     = $gpgrp
    @args.push @gprhs     = $gprhs
    @args.push @gprhc     = $gprhc
    @args.push @gpvel     = $gpvel
    @args.push @gpsgs     = $gpsgs
    @args.push @elvel     = $elvel
    @args.push @elpre     = $elpre
    @args.push @elbub     = $elbub
  
    @args.push @wgrgr     = $wgrgr
    @args.push @agrau     = $agrau
  
    # Matrices
    @args.push @elauu     = $elauu
    @args.push @elaup     = $elaup
    @args.push @elapp     = $elapp
    @args.push @elapu     = $elapu
    @args.push @elrbu     = $elrbu
    @args.push @elrbp     = $elrbp
    # Others
    @args.push @dtinv_loc = $dtinv_loc
    @args.push @dtsgs     = $dtsgs
    @args.push @pbubl     = $pbubl
    @args.push @gpsha_bub = $gpsha_bub
    @args.push @gpcar_bub = $gpcar_bub
    # Enrichement Element matrices
    @args.push @elauq     = $elauq
    @args.push @elapq     = $elapq
    @args.push @elaqu     = $elaqu
    @args.push @elaqp     = $elaqp
    @args.push @elaqq     = $elaqq
    @args.push @elrbq     = $elrbq
  end

  def initialize(options,ndime = nil)
    @opts = {:vector_length => 1, :preprocessor => false, :nests => (1..10).to_a, :unroll => false, :inline => :included}
    @opts.update(options)
    
    declare_parameters(ndime)
  end

    def declare_procedure(name, functions = nil)
      return Procedure(name, @args, :functions => functions)
    end
end
