require '../Common/Parameters.rb'

class KASGS
  include Parameters
  attr_reader :kernel

  def initialize(options)
    @opts = {:vector_length => 1, :preprocessor => false, :nests => (1..10).to_a, :unroll => false, :inline => :included}
    @opts.update(options)
    
    Parameters.initialize(@opts[:vector_length])

    @args = []
    @args.push @ndime              = $ndime
    @args.push @kfl_lumped         = $kfl_lumped
    @args.push @mnode              = $mnode
    @args.push @ntens              = $ntens
    @args.push @kfl_duatss         = $kfl_duatss
    @args.push @fact_duatss        = $fact_duatss
    @args.push @kfl_stabi_nsi      = $kfl_stabi_nsi
    @args.push @fvins_nsi          = $fvins_nsi
    @args.push @fcons_nsi          = $fcons_nsi
    @args.push @bemol_nsi          = $bemol_nsi
    @args.push @kfl_regim_nsi      = $kfl_regim_nsi
    @args.push @fvela_nsi          = $fvela_nsi
    @args.push @kfl_rmom2_nsi      = $kfl_rmom2_nsi
    @args.push @kfl_press_nsi      = $kfl_press_nsi
    @args.push @kfl_p1ve2_nsi      = $kfl_p1ve2_nsi
    @args.push @kfl_linea_nsi      = $kfl_linea_nsi
    @args.push @pabdf_nsi          = $pabdf_nsi
    @args.push @kfl_confi_nsi      = $kfl_confi_nsi
    @args.push @nbdfp_nsi          = $nbdfp_nsi
    @args.push @kfl_sgsti_nsi      = $kfl_sgsti_nsi
    @args.push @kfl_nota1_nsi      = $kfl_nota1_nsi
    @args.push @kfl_limit_nsi      = $kfl_limit_nsi
    @args.push @kfl_penal_nsi      = $kfl_penal_nsi
    @args.push @penal_nsi          = $penal_nsi
    @args.push @kfl_bubbl_nsi      = $kfl_bubbl_nsi

    @args.push @pnode              = $pnode
    @args.push @pgaus              = $pgaus
    @args.push @gpden              = $gpden
    @args.push @gpvis              = $gpvis
    @args.push @gppor              = $gppor
    @args.push @gpgvi              = $gpgvi
    @args.push @gpsp1              = $gpsp1
    @args.push @gptt1              = $gptt1
    @args.push @gpsp2              = $gpsp2
    @args.push @gptt2              = $gptt2
    @args.push @gpvol              = $gpvol
    @args.push @gpsha              = $gpsha
    @args.push @gpcar              = $gpcar
    @args.push @gpadv              = $gpadv
    @args.push @gpvep              = $gpvep
    @args.push @gplap              = $gplap
    @args.push @gphes              = $gphes
    @args.push @gprhs              = $gprhs
    @args.push @gprhs_sgs          = $gprhs_sgs
    @args.push @gprhc              = $gprhc
    @args.push @gprh2              = $gprh2
    @args.push @rmomu              = $rmomu
    @args.push @rcont              = $rcont
    @args.push @elpre              = $elpre
    @args.push @elbub              = $elbub

    # Element matrices
    @args.push @elauu              = $elauu
    @args.push @elaup              = $elaup
    @args.push @elapp              = $elapp
    @args.push @elapu              = $elapu
    @args.push @elrbu              = $elrbu
    @args.push @elrbp              = $elrbp
    @args.push @rmom2              = $rmom2
    @args.push @gpst1              = $gpst1
    @args.push @gpgve              = $gpgve
    @args.push @elvel              = $elvel
    @args.push @ellum              = $ellum
    @args.push @dtinv_loc          = $dtinv_loc
    @args.push @pbubl              = $pbubl
    @args.push @gpsha_bub          = $gpsha_bub                    # Ne
    @args.push @gpcar_bub          = $gpcar_bub              # dNe/dxi
    @args.push @gppre              = $gppre

    # Enrichement Element matrices
    @args.push @elauq              = $elauq
    @args.push @elapq              = $elapq
    @args.push @elaqu              = $elaqu
    @args.push @elaqp              = $elaqp
    @args.push @elaqq              = $elaqq
    @args.push @elrbq              = $elrbq
    
    @args.push @wgrgr              = $wgrgr
  end
  def declare_procedure(functions = nil)
    return Procedure("nsi_element_assembly_asgs_oss", @args, :functions => functions)
  end
end
