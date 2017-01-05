module CommonArgs
  def self.init(pgaus,pnode,ndime,vector_length,nb_instances,seed = nil)
    ANArray.srand(seed) unless seed.nil?
    @@pgaus = pgaus
    @@pnode = pnode
    @@ndime = ndime
 
    @@mnode     = 8
    @@ntens     = 3*@@ndime-3
    @@ncomp_nsi = 3

    @@kfl_lumped    = 0
    @@kfl_duatss    = 0
    @@fact_duatss   = 0.0
    @@kfl_stabi_nsi = 2
    @@fvins_nsi     = 1.0

    @@fcons_nsi     = 0.0
    @@bemol_nsi     = 0.0
    @@kfl_regim_nsi = 0
    @@fvela_nsi     = ANArray.float(8*vector_length,3).random!
    @@kfl_rmom2_nsi = 0
    @@kfl_press_nsi = 1
    @@kfl_p1ve2_nsi = 0
    @@kfl_linea_nsi = 2
    @@pabdf_nsi     = 1.0
    @@kfl_confi_nsi = 0
    @@nbdfp_nsi     = @@ncomp_nsi-1
    @@kfl_sgsti_nsi = 1
    @@kfl_nota1_nsi = 0
    @@kfl_limit_nsi = 0
    @@kfl_penal_nsi = 0
    @@penal_nsi     = 0.0
    @@kfl_bubbl_nsi = 0

    @@gpden     = ANArray.float(8*vector_length,vector_length,@@pgaus).random!
    @@gpvis     = ANArray.float(8*vector_length,vector_length,@@pgaus).random!
    @@gppor     = ANArray.float(8*vector_length,vector_length,@@pgaus).random! 
    @@gpgvi     = ANArray.float(8*vector_length,vector_length,@@ndime,@@pgaus).random!
    @@gpsp1     = ANArray.float(8*vector_length,vector_length,@@pgaus).random! 
    @@gptt1     = ANArray.float(8*vector_length,vector_length,@@pgaus).random! 
    @@gpsp2     = ANArray.float(8*vector_length,vector_length,@@pgaus).random!
    @@gptt2     = ANArray.float(8*vector_length,vector_length,@@pgaus).random! 
    @@gpvol     = ANArray.float(8*vector_length,vector_length,@@pgaus).random! 
    @@gpsha     = ANArray.float(8*vector_length,vector_length,@@pnode,@@pgaus).random! 
    @@gpcar     =  ANArray.float(8*vector_length,vector_length,@@ndime,@@mnode,@@pgaus).random!
    @@gplap     = ANArray.float(8*vector_length,vector_length,@@pnode,@@pgaus).random! 
    @@gphes     = ANArray.float(8*vector_length,vector_length,@@ntens,@@mnode,@@pgaus).random! 
    @@gpadv     = ANArray.float(8*vector_length,vector_length,@@ndime,@@pgaus).random!
    @@rmomu     =  ANArray.float(8*vector_length,vector_length,@@pnode,@@pgaus).random!
    @@rcont     = ANArray.float(8*vector_length,vector_length,@@ndime,@@pnode,@@pgaus).random!
    @@elpre     = ANArray.float(8*vector_length,vector_length,@@pnode,10).random!
    @@elbub     = ANArray.float(8*vector_length,vector_length).random!
    @@rmom2     = ANArray.float(8*vector_length,vector_length,@@ndime,@@ndime,@@pnode,@@pgaus).random!
    @@gpst1     = ANArray.float(8*vector_length,vector_length,@@pgaus).random!
    @@gpgve     = ANArray.float(8*vector_length,vector_length,@@ndime,@@ndime,@@pgaus).random!
    @@gprh2     = ANArray.float(8*vector_length,vector_length,@@pgaus).random!
    @@elvel = ANArray.float(8*vector_length,vector_length,@@ndime,@@pnode,10).random!
    @@ellum = ANArray.float(8*vector_length,vector_length,@@pnode).random!
    @@dtinv_loc = ANArray.float(8*vector_length,vector_length).random!
    @@dtsgs = ANArray.float(8*vector_length,vector_length).random!

    @@pbubl = ANArray.int(4*vector_length,vector_length).random!
    @@gpsha_bub = ANArray.float(8*vector_length,vector_length,@@pgaus).random!
    @@gpcar_bub = ANArray.float(8*vector_length,vector_length,@@ndime,@@pgaus).random!
    @@gppre = ANArray.float(8*vector_length,vector_length,@@pgaus).random!
    @@gpvel = ANArray.float(8*vector_length,vector_length,@@ndime,@@pgaus,10).random!
    @@gpsgs = ANArray.float(8*vector_length,vector_length,@@ndime,@@pgaus,10).random!

    @@wgrgr = []
    @@agrau = []

    @@elauu = []
    @@elaup = []
    @@elapp = []
    @@elapu = []
    @@elrbu = []
    @@elrbp = []

    @@elauq = [] 
    @@elapq = []
    @@elaqu = []
    @@elaqp = []
    @@elaqq = []
    @@elrbq = []

    @@gpgrp = []
    @@gprhs = []
    @@gpvep = []
    @@gprhc = []
    @@gprhs_sgs = []

    nb_instances.times{|i|
      @@wgrgr.push ANArray.float(8*vector_length,vector_length,@@pnode,@@pnode,@@pgaus)
      @@agrau.push ANArray.float(8*vector_length,vector_length,@@pnode,@@pgaus)

      @@elauu.push ANArray.float(8*vector_length,vector_length,@@pnode*@@ndime,@@pnode*@@ndime)
      @@elaup.push ANArray.float(8*vector_length,vector_length,@@pnode*@@ndime,@@pnode)
      @@elapp.push ANArray.float(8*vector_length,vector_length,@@pnode,@@pnode)
      @@elapu.push ANArray.float(8*vector_length,vector_length,@@pnode,@@pnode*@@ndime)
      @@elrbu.push ANArray.float(8*vector_length,vector_length,@@ndime,@@pnode)
      @@elrbp.push ANArray.float(8*vector_length,vector_length,@@pnode)

      @@elauq.push ANArray.float(8*vector_length,vector_length,@@pnode*@@ndime,1)
      @@elapq.push ANArray.float(8*vector_length,vector_length,@@pnode,1)
      @@elaqu.push ANArray.float(8*vector_length,vector_length,1,@@pnode*@@ndime)
      @@elaqp.push ANArray.float(8*vector_length,vector_length,1,@@pnode)
      @@elaqq.push ANArray.float(8*vector_length,vector_length,1,1)
      @@elrbq.push ANArray.float(8*vector_length,vector_length,1)

      case i
      when 0
        @@gpgrp.push ANArray.float(8*vector_length,vector_length,@@ndime,@@pgaus).random! 
        @@gprhs.push ANArray.float(8*vector_length,vector_length,@@ndime,@@pgaus).random!
        @@gpvep.push ANArray.float(8*vector_length,vector_length,@@ndime,@@pgaus).random! 
        @@gprhc.push ANArray.float(8*vector_length,vector_length,@@pgaus).random! 
        @@gprhs_sgs.push ANArray.float(8*vector_length,vector_length,@@ndime,@@pgaus).random!

      else
        @@gpgrp.push @@gpgrp[0].clone
        @@gprhs.push @@gprhs[0].clone
        @@gpvep.push @@gpvep[0].clone
        @@gprhc.push @@gprhc[0].clone
        @@gprhs_sgs.push @@gprhs_sgs[0].clone
      end
    }
  end
end
