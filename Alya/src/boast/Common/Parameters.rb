
require 'BOAST'
include BOAST

module Parameters
  $ndime         = Int("ndime", :dir => :in )
  $kfl_lumped    = Int("kfl_lumped",    :dir => :in)
  $mnode         = Int("mnode",         :dir => :in)
  $ntens         = Int("ntens",         :dir => :in)
  $kfl_duatss    = Int("kfl_duatss",    :dir => :in)
  $fact_duatss   = Int("fact_duatss",   :dir => :in)
  $kfl_stabi_nsi = Int("kfl_stabi_nsi", :dir => :in)
  $fvins_nsi     = Real("fvins_nsi",    :dir => :in)
  $fcons_nsi     = Real("fcons_nsi",    :dir => :in)
  $bemol_nsi     = Real("bemol_nsi",    :dir => :in)
  $kfl_regim_nsi = Int("kfl_regim_nsi", :dir => :in)
  $fvela_nsi     = Real("fvela_nsi",    :dir => :in, :dim => [Dim(3)])
  $kfl_rmom2_nsi = Int("kfl_rmom2_nsi", :dir => :in)
  $kfl_press_nsi = Int("kfl_press_nsi", :dir => :in)
  $kfl_p1ve2_nsi = Int("kfl_p1ve2_nsi", :dir => :in)
  $kfl_linea_nsi = Int("kfl_linea_nsi", :dir => :in)
  $pabdf_nsi     = Real("pabdf_nsi",    :dir => :in, :dim => [Dim()])
  $kfl_confi_nsi = Int("kfl_confi_nsi", :dir => :in)
  $nbdfp_nsi     = Int("nbdfp_nsi",     :dir => :in)
  $kfl_sgsti_nsi = Int("kfl_sgsti_nsi", :dir => :in)
  $kfl_nota1_nsi = Int("kfl_nota1_nsi", :dir => :in)
  $kfl_limit_nsi = Int("kfl_limit_nsi", :dir => :in)
  $kfl_penal_nsi = Int("kfl_penal_nsi", :dir => :in)
  $penal_nsi     = Real("penal_nsi",    :dir => :in)
  $kfl_bubbl_nsi = Int("kfl_bubbl_nsi", :dir => :in)

  $pnode     = Int("pnode",             :dir => :in)
  $pgaus     = Int("pgaus",             :dir => :in)

  $inode     = Int("inode")
  $jnode     = Int("jnode")
  $jdime     = Int("jdime")
  $idofv     = Int("idofv")
  $ivect     = Int("ivect")
  $igaus     = Int("igaus")
  $idime     = Int("idime")
  $jdofv     = Int("jdofv")
  $itime     = Int("itime")
  
  def self.initialize(vector_length, ndime = nil)
    allocate = get_lang == C ? true : false
    ndime_var  = ndime.nil? ? $ndime : ndime

    $gpden     = Real("gpden",     :dir => :in,     :vector_length => vector_length, :dim => [Dim($pgaus)])
    $gpvis     = Real("gpvis",     :dir => :in,     :vector_length => vector_length, :dim => [Dim($pgaus)])              
    $gppor     = Real("gppor",     :dir => :in,     :vector_length => vector_length, :dim => [Dim($pgaus)])
    $gpgvi     = Real("gpgvi",     :dir => :in,     :vector_length => vector_length, :dim => [Dim(ndime_var),Dim($pgaus)])
    $gptt1     = Real("gptt1",     :dir => :in,     :vector_length => vector_length, :dim => [Dim($pgaus)])
    $gptt2     = Real("gptt2",     :dir => :in,     :vector_length => vector_length, :dim => [Dim($pgaus)])
    $gplap     = Real("gplap",     :dir => :in,     :vector_length => vector_length, :dim => [Dim($pnode),Dim($pgaus)])
    $gphes     = Real("gphes",     :dir => :in,     :vector_length => vector_length, :dim => [Dim($ntens),Dim($mnode),Dim($pgaus)])
    $gprhs_sgs = Real("gprhs_sgs", :dir => :inout,  :vector_length => vector_length, :dim => [Dim(ndime_var),Dim($pgaus)])
    $gprh2     = Real("gprh2",     :dir => :in,     :vector_length => vector_length, :dim => [Dim($pgaus)])
    $rmomu     = Real("rmomu",     :dir => :in,     :vector_length => vector_length, :dim => [Dim($pnode),Dim($pgaus)])
    $rcont     = Real("rcont",     :dir => :in,     :vector_length => vector_length, :dim => [Dim(ndime_var),Dim($pnode),Dim($pgaus)])
    $gpsp1     = Real("gpsp1",     :dir => :in,     :vector_length => vector_length, :dim => [Dim($pgaus)])
    $gpsp2     = Real("gpsp2",     :dir => :in,     :vector_length => vector_length, :dim => [Dim($pgaus)])
    $gpvol     = Real("gpvol",     :dir => :in,     :vector_length => vector_length, :dim => [Dim($pgaus)])
    $gpsha     = Real("gpsha",     :dir => :in,     :vector_length => vector_length, :dim => [Dim($pnode),Dim($pgaus)])
    $gpcar     = Real("gpcar",     :dir => :in,     :vector_length => vector_length, :dim => [Dim(ndime_var),Dim($mnode),Dim($pgaus)])
    $gpadv     = Real("gpadv",     :dir => :in,     :vector_length => vector_length, :dim => [Dim(ndime_var),Dim($pgaus)])
    $gpvep     = Real("gpvep",     :dir => :inout , :vector_length => vector_length, :dim => [Dim(ndime_var),Dim($pgaus)])
    $gpgrp     = Real("gpgrp",     :dir => :inout , :vector_length => vector_length, :dim => [Dim(ndime_var),Dim($pgaus)])
    $gprhs     = Real("gprhs",     :dir => :inout , :vector_length => vector_length, :dim => [Dim(ndime_var),Dim($pgaus)])
    $gprhc     = Real("gprhc",     :dir => :inout , :vector_length => vector_length, :dim => [Dim($pgaus)])
    $gpvel     = Real("gpvel",     :dir => :in,     :vector_length => vector_length, :dim => [Dim(ndime_var),Dim($pgaus),Dim()])
    $gpsgs     = Real("gpsgs",     :dir => :in,     :vector_length => vector_length, :dim => [Dim(ndime_var),Dim($pgaus),Dim()])
    $elvel     = Real("elvel",     :dir => :in,     :vector_length => vector_length, :dim => [Dim(ndime_var),Dim($pnode),Dim()])
    $elpre     = Real("elpre",     :dir => :in,     :vector_length => vector_length, :dim => [Dim($pnode),Dim()])
    $elbub     = Real("elbub",     :dir => :in,     :vector_length => vector_length, :dim => [Dim(1)])

    $wgrgr     = Real("wgrgr",     :dir => :out,    :vector_length => vector_length, :dim => [Dim($pnode),Dim($pnode),Dim($pgaus)])
    $agrau     = Real("agrau",     :dir => :out,    :vector_length => vector_length, :dim => [Dim($pnode),Dim($pgaus)])

    # Matrices
    $elauu = Real("elauu", :dir => :out, :vector_length => vector_length, :dim => [Dim($pnode*ndime_var),Dim($pnode*ndime_var)])
    $elaup = Real("elaup", :dir => :out, :vector_length => vector_length, :dim => [Dim($pnode*ndime_var),Dim($pnode)])
    $elapp = Real("elapp", :dir => :out, :vector_length => vector_length, :dim => [Dim($pnode),Dim($pnode)])
    $elapu = Real("elapu", :dir => :out, :vector_length => vector_length, :dim => [Dim($pnode),Dim($pnode*ndime_var)])
    $elrbu = Real("elrbu", :dir => :out, :vector_length => vector_length, :dim => [Dim(ndime_var),Dim($pnode)])
    $elrbp = Real("elrbp", :dir => :out, :vector_length => vector_length, :dim => [Dim($pnode)])
    $rmom2 = Real("rmom2", :dir => :in,  :vector_length => vector_length, :dim => [Dim(ndime_var),Dim(ndime_var),Dim($pnode),Dim($pgaus)])
    $gpst1 = Real("gpst1", :dir => :in,  :vector_length => vector_length, :dim => [Dim($pgaus)])
    $gpgve = Real("gpgve", :dir => :in,  :vector_length => vector_length, :dim => [Dim(ndime_var),Dim(ndime_var),Dim($pgaus)])
    $ellum = Real("ellum", :dir => :out, :vector_length => vector_length, :dim => [Dim($pnode)])
    $gppre = Real("gppre", :dir => :in,  :vector_length => vector_length, :dim => [Dim($pgaus)]) 

    # Others
    $dtinv_loc = Real("dtinv_loc", :dir => :in, :vector_length => vector_length, :dim => [Dim(1)])
    $dtsgs     = Real("dtsgs",     :dir => :in, :vector_length => vector_length, :dim => [Dim(1)])
    $pbubl     = Int("pbubl",      :dir => :in, :vector_length => vector_length, :dim => [Dim(1)])
    $gpsha_bub = Real("gpsha_bub", :dir => :in, :vector_length => vector_length, :dim => [Dim($pgaus)])
    $gpcar_bub = Real("gpcar_bub", :dir => :in, :vector_length => vector_length, :dim => [Dim(ndime_var),Dim($pgaus)])

    # Enrichement Element matrices
    $elauq     = Real("elauq",     :dir => :out, :vector_length => vector_length, :dim => [Dim($pnode*ndime_var),Dim(1)])
    $elapq     = Real("elapq",     :dir => :out, :vector_length => vector_length, :dim => [Dim($pnode),Dim(1)])
    $elaqu     = Real("elaqu",     :dir => :out, :vector_length => vector_length, :dim => [Dim(1),Dim($pnode*ndime_var)])
    $elaqp     = Real("elaqp",     :dir => :out, :vector_length => vector_length, :dim => [Dim(1),Dim($pnode)])
    $elaqq     = Real("elaqq",     :dir => :out, :vector_length => vector_length, :dim => [Dim(1),Dim(1)])
    $elrbq     = Real("elrbq",     :dir => :out, :vector_length => vector_length, :dim => [Dim(1)])

    # Locals
    $ivect  = Int("ivect")
    $kdime  = Int("kdime")
    $xvis2  = Real("xvis2")
    $xvisc  = Real("xvisc")
    $one_rp = Real("one_rp")
    $p1ve2  = Real("p1ve2",        :vector_length => vector_length, :dim => [Dim(ndime_var),Dim(ndime_var),Dim($pnode),Dim($pgaus)], :allocate => allocate)
    $p1vec  = Real("p1vec",        :vector_length => vector_length, :dim => [Dim($pnode),Dim($pgaus)], :allocate => allocate)
    $p2vec  = Real("p2vec",        :vector_length => vector_length, :dim => [Dim(ndime_var),Dim($pnode),Dim($pgaus)], :allocate => allocate)
    $p2sca  = Real("p2sca",        :vector_length => vector_length, :dim => [Dim($pnode),Dim($pgaus)], :allocate => allocate)
    $wgrvi  = Real("wgrvi",        :vector_length => vector_length, :dim => [Dim($pnode),Dim($pgaus)], :allocate => allocate)

    $factx     = Real("factx",     :vector_length => vector_length)
    $facty     = Real("facty",     :vector_length => vector_length)
    $facx1     = Real("facx1",     :vector_length => vector_length)
    $facy1     = Real("facy1",     :vector_length => vector_length)
    $ugraN     = Real("ugraN",     :vector_length => vector_length)
    $gramugraN = Real("gramugraN", :vector_length => vector_length)
    $penal     = Real("penal",     :vector_length => vector_length)
    $gprhh     = Real("gprhh",     :vector_length => vector_length, :dim => [Dim(ndime_var),Dim($pgaus)], :allocate => allocate)
    $taupr     = Real("taupr",     :vector_length => vector_length, :dim => [Dim($pgaus)], :allocate => allocate)
    $gpveo     = Real("gpveo",     :vector_length => vector_length, :dim => [Dim(ndime_var)], :allocate => allocate)
    $gpcar1ji  = Real("gpcar1ji",  :vector_length => vector_length)
    $gpcar2ji  = Real("gpcar2ji",  :vector_length => vector_length)
    $gpcar3ji  = Real("gpcar3ji",  :vector_length => vector_length)
    $p2sca_bub = Real("p2sca_bub", :vector_length => vector_length, :dim => [Dim($pgaus)], :allocate => allocate)
    $p2vec_bub = Real("p2vec_bub", :vector_length => vector_length, :dim => [Dim(ndime_var),Dim($pgaus)], :allocate => allocate)

    for i in 1..3 do
      for j in 1..3 do 
        eval '$elauu#{i}#{j} = Real("elauu#{i}#{j}", :vector_length => vector_length, :dim => [Dim(4*((pnode+3)/4)), Dim($pnode)], :allocate => allocate)'
      end
    end

    for i in 1..4 do
      eval '$factvec#{i} = Real("factvec#{i}", :vector_length => vector_length, :dim=> [Dim(4*((pnode+3)/4))], :allocate => allocate)'
    end

    $gpsp1_p   = Real("gpsp1_p",   :vector_length => vector_length, :dim => [Dim($pgaus)], :allocate => allocate)
    $gpsp1_v   = Real("gpsp1_v",   :vector_length => vector_length, :dim => [Dim($pgaus)], :allocate => allocate)
    $gpsp2_v   = Real("gpsp2_v",   :vector_length => vector_length, :dim => [Dim($pgaus)], :allocate => allocate)
    $c1        = Real("c1",        :vector_length => vector_length, :dim => [Dim(1)],      :allocate => allocate)
    $c2        = Real("c2",        :vector_length => vector_length, :dim => [Dim(1)],      :allocate => allocate)
    $c3        = Real("c3",        :vector_length => vector_length, :dim => [Dim(1)],      :allocate => allocate)
    $c4        = Real("c4",        :vector_length => vector_length, :dim => [Dim(1)],      :allocate => allocate)
    $alpha     = Real("alpha",     :vector_length => vector_length, :dim => [Dim(1)],      :allocate => allocate)
    $beta      = Real("beta",      :vector_length => vector_length, :dim => [Dim(1)],      :allocate => allocate)

    $gpveo     = Real("gpveo",     :vector_length => vector_length, :dim => [Dim(3)],      :allocate => allocate)
    $fact1_p   = Real("fact1_p",   :vector_length => vector_length, :dim => [Dim(1)],      :allocate => allocate)
    $dtinv_mod = Real("dtinv_mod", :vector_length => vector_length, :dim => [Dim(1)],      :allocate => allocate)
    $fact      = Real('fact',      :vector_length => vector_length, :dim => [Dim(9)],      :allocate => allocate)

    $idof      = Int("idof", :dim => [Dim(3)], :allocate => allocate)
    $jdof      = Int("jdof", :dim => [Dim(3)], :allocate => allocate)
  end

  def self.copy(arg, direction = nil)
    if direction == arg.direction or direction.nil? then
      return arg
    else
      (copy = arg.clone).direction = direction
      return copy
    end
  end
end
