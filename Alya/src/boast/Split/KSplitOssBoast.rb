
require_relative '../Common/subroutine.rb'
require_relative './KSplitOss.rb'

class Nest1 < Subroutine
  def initialize(options,functions = nil)
    super("nest1",options,functions)
    # ins
    @mnode = Parameters.copy($mnode,:in)
    @pnode = Parameters.copy($pnode,:in)
    @pgaus = Parameters.copy($pgaus,:in)
    @gpden = Parameters.copy($gpden,:in)
    @gpcar = Parameters.copy($gpcar,:in)
    @gpadv = Parameters.copy($gpadv,:in)

    # inouts
    @wgrgr = Parameters.copy($wgrgr,:inout)
    @agrau = Parameters.copy($agrau,:inout)

    @args.push @mnode,@pnode,@pgaus,@gpden,@gpcar,@gpadv,@wgrgr,@agrau
  end
 
  def generate
    # locals
    inode = $inode
    jnode = $jnode
    igaus = $igaus

    block = lambda{
      decl inode,jnode,igaus unless @usage == :included

      exp1 = ""
      exp2 = ""
      @ndime.times{|i|
        exp1 += "@gpadv[#{i+1},igaus]*@gpcar[#{i+1},inode,igaus]"
        exp2 += "@gpcar[#{i+1},inode,igaus]*@gpcar[#{i+1},jnode,igaus]"
        exp1 +=  " + " if i+1 < @ndime
        exp2 +=  " + " if i+1 < @ndime
      }
      
      form_agrau = "pr @agrau[inode,igaus] === @gpden[igaus] * (#{exp1})"
      form_wgrgr = "pr @wgrgr[inode,jnode,igaus] === #{exp2}"
      
      pr For(igaus,1,@pgaus){
        pr For(inode,1,@pnode){
          eval form_agrau
          pr For(jnode,1,@pnode){                    
            eval form_wgrgr
          }
        }
      }
    }
    
    construct(block)
  end
end

class Nest2 < Subroutine
  def initialize(options, functions = nil)
    super("nest2",options,functions)

    # ins
    @mnode = Parameters.copy($mnode,:in)
    @pnode = Parameters.copy($pnode,:in)
    @pgaus = Parameters.copy($pgaus,:in)
    @gpden = Parameters.copy($gpden,:in)
    @gpcar = Parameters.copy($gpcar,:in)
    @gpsp2 = Parameters.copy($gpsp2,:in)
    @gpvol = Parameters.copy($gpvol,:in)
    @gpvis = Parameters.copy($gpvis,:in)
    @dtinv_loc = Parameters.copy($dtinv_loc,:in)
    @gppor = Parameters.copy($gppor,:in)
    @gpsha = Parameters.copy($gpsha,:in)
    @gpsp1_v = Parameters.copy($gpsp1_v,:in)
    @wgrgr = Parameters.copy($wgrgr,:in)
    @agrau = Parameters.copy($agrau,:in)

    #inouts
    @elauu = Parameters.copy($elauu,:inout)

    @args.push @mnode,@pnode,@pgaus,@gpden,@gpcar,@gpsp2,@gpvol,@gpvis,@dtinv_loc,
                @gppor,@gpsha,@gpsp1_v,@wgrgr,@agrau,@elauu
  end
  def generate
    # locals
    fact = Parameters.copy($fact)
    igaus = Parameters.copy($igaus)
    inode = Parameters.copy($inode)
    jnode = Parameters.copy($jnode)
    idof = Parameters.copy($idof)
    jdof = Parameters.copy($jdof)
    idime = Parameters.copy($idime)
    jdime = Parameters.copy($jdime)

    for1 = For(igaus,1,@pgaus){ 
      pr fact[1] === @gpsp2[igaus] * @gpvol[igaus]
      pr fact[7] === @gpvis[igaus] * @gpvol[igaus]
      pr fact[8] === @gpsp1_v[igaus] * @gpvol[igaus]
      pr fact[9] ===  @gpden[igaus] * @functions[:pabdf_nsi].call(1) * @dtinv_loc[1] + @gppor[igaus]
    
      pr For(inode,1,@pnode){
        pr For(idime,1,@ndime){
          pr idof[1] === (inode-1)*@ndime+idime
          pr fact[2] === fact[1] * @gpcar[idime,inode,igaus]      
          pr For(jnode,1,@pnode){
            pr For(jdime,1,@ndime){
              pr jdof[jdime]              === (jnode-1)*@ndime+jdime
              pr @elauu[idof[1],jdof[jdime]] === @elauu[idof[1],jdof[jdime]] + fact[2] * @gpcar[jdime,jnode,igaus]                   
            }.unroll(@unroll)
    
            pr jdof[1]              === (jnode-1)*@ndime+idime
            pr fact[5]              === @gpsha[inode,igaus] * @gpvol[igaus]
            pr fact[6]              === fact[5] * ( @agrau[jnode,igaus] + fact[9] * @gpsha[jnode,igaus] ) + fact[7] * @wgrgr[inode,jnode,igaus] + fact[8] *   @agrau[jnode,igaus] * @agrau[inode,igaus]
            pr @elauu[idof[1],jdof[1]] === @elauu[idof[1],jdof[1]] + fact[6]
    
          }
        }
      }
    }
    
    main_block = lambda{
      decl fact,igaus,inode,jnode,idof,jdof,idime,jdime unless @usage == :included
      pr for1
    }
          

    construct(main_block)
  end
end

class Nest3 < Subroutine
    def initialize(options, functions = nil)
      super("nest3",options,functions)
      
      # ins
      @pgaus = Parameters.copy($pgaus,:in)
      @pnode = Parameters.copy($pnode,:in)
      @mnode = Parameters.copy($mnode,:in)
      @gpvis = Parameters.copy($gpvis,:in)
      @gpvol = Parameters.copy($gpvol,:in)
      @gpcar = Parameters.copy($gpcar,:in)
      @fvins_nsi = Parameters.copy($fvins_nsi,:in)

      # inouts
      @elauu = Parameters.copy($elauu,:inout)

      @args.push @pgaus,@pnode,@mnode,@gpvis,@gpvol,@gpcar,@fvins_nsi,@elauu
    end
  def generate
    # locals
    igaus = $igaus
    idime = $idime
    jdime = $jdime
    inode = $inode
    jnode = $jnode
    fact = $fact
    idofv = $idofv
    jdofv = $jdofv
    
    for_inode = For(inode,1,@pnode)
    for_jnode = For(jnode,1,@pnode)
    for_idime = For(idime,1,@ndime)
    for_jdime = lambda{ |dim| return For(jdime,1,dim) }
    
    exp_fact1 = fact[2] === @gpvis[igaus] * @gpvol[igaus] * @gpcar[idime,jnode,igaus]
    exp_idofv = idofv === (inode-1)*@ndime + idime
    
    block_elauu_1 = lambda{
      pr jdofv === (jnode-1)*@ndime + jdime
      pr @elauu[idofv,jdofv] === @elauu[idofv,jdofv] + fact[2] * @gpcar[jdime,inode,igaus]
    }
    
    exp_fact1_2 = fact[2] === -2.0.to_var / 3.0.to_var * @gpvis[igaus] * @gpvol[igaus] * @gpcar[idime,inode,igaus]
    
    block_elauu_2 = lambda{
      pr jdofv === (jnode-1)*@ndime + jdime
      pr @elauu[idofv,jdofv] === @elauu[idofv,jdofv] + fact[2] * @gpcar[jdime,jnode,igaus]
    }

    main_block = lambda{ 
      decl igaus,idime,jdime,inode,jnode,fact,idofv,jdofv unless @usage == :included
    
      pr If(@fvins_nsi > 0.9){
        for1 = For(igaus,1,@pgaus){ 
    
          opn for_inode
          opn for_idime
          pr exp_idofv
          opn for_jnode
          pr exp_fact1
          pr for_jdime.call(@ndime).unroll(@unroll), &block_elauu_1
          close for_jnode
          pr If(@fvins_nsi == 2.0){
            pr exp_fact1_2
            opn for_jnode
            pr for_jdime.call(@ndime).unroll(@unroll), &block_elauu_2
            close for_jnode
          }
          close for_idime
          close for_inode
        }
        
        pr for1
      }
    }              

    construct(main_block)
  end
end

class Nest4 < Subroutine
  def initialize(options, functions = nil)
    super("nest4",options,functions)

    # ins
    @pgaus = Parameters.copy($pgaus,:in)
    @pnode = Parameters.copy($pnode,:in)
    @gpvol = Parameters.copy($gpvol,:in)
    @gpden = Parameters.copy($gpden,:in)
    @dtinv_mod = Parameters.copy($dtinv_mod,:in)
    @gpsha = Parameters.copy($gpsha,:in)
    @elvel = Parameters.copy($elvel,:in)
    @kfl_lumped = Parameters.copy($kfl_lumped,:in)

    # inouts
    @gpveo = Parameters.copy($gpveo,:inout)
    @elrbu = Parameters.copy($elrbu,:inout)
    @elauu = Parameters.copy($elauu,:inout)

    @args.push @pgaus,@pnode,@gpvol,@gpden,@dtinv_mod,@gpsha,@elvel,@kfl_lumped,@gpveo,@elrbu,@elauu
  end
  def generate
    # locals 
    idime = $idime
    jdime = $jdime
    igaus = $igaus
    inode = $inode 
    jnode = $jnode
    fact = $fact
    idof = $idof
    jdof = $jdof

    exp_fact0 = fact[1] === @gpvol[igaus] * @gpden[igaus] * @dtinv_mod[1]            
    exp_fact1 = fact[2] === fact[1] * @gpsha[inode,igaus]
    exp_elrbu = @elrbu[idime,inode] === @elrbu[idime,inode] + fact[2] * @elvel[idime,inode,2]
    for_elauu1 = For(idime,1,@ndime){
      pr idof[idime] === (inode-1)*@ndime+idime
      pr @elauu[idof[idime],idof[idime]] === @elauu[idof[idime],idof[idime]] + fact[2]
    }
    for_elrbu2 = For(idime,1,@ndime){
      pr exp_elrbu
    }
    
    for_elauu2 = For(jdime,1,3){
      pr jdof[jdime] === (jnode-1)*@ndime+jdime
      pr @elauu[idof[jdime],jdof[jdime]] === @elauu[idof[jdime],jdof[jdime]] - fact[2] * @gpsha[jnode,igaus]
    }
    
    for_elrbu1 = For(idime,1,3){
      pr @elrbu[idime,inode] ===  @elrbu[idime,inode] - fact[2] * @gpveo[idime]
      pr exp_elrbu
    }              
        block_kfl_lumped_1 = lambda{
          pr For(igaus,1,@pgaus){
            if get_lang == FORTRAN then
              pr @gpveo === 0.0.to_var
            elsif get_lang == C
              code =<<EOF
              memset(gpveo,0,sizeof(#{@gpveo.type.decl}) * #{@gpveo.dimension[0]});
EOF
              get_output.print code
            end
            pr For(inode,1,@pnode){
              pr For(idime,1,3){
                pr @gpveo[idime] === @gpveo[idime] + @elvel[idime,inode,2] * @gpsha[inode,igaus] 
              }.unroll(@unroll)

            }
            pr exp_fact0
            pr For(inode,1,@pnode){
              pr exp_fact1
              if @unroll then
                pr for_elauu1.unroll
                pr for_elrbu1.unroll
              else
                pr for_elauu1
                pr for_elrbu1
              end
              pr For(jnode,1,@pnode){
                pr for_elauu2.unroll(@unroll)
              }
            }
          }
        }
    block_kfl_lumped_2 = lambda{
        pr For(igaus,1,@pgaus){
          pr exp_fact0
          pr For(inode,1,@pnode){
            pr exp_fact1
            pr for_elauu1.unroll(@unroll)
            pr for_elrbu2.unroll(@unroll)
          }
        }
      }
    if_kfl_lumped_1 = @ndime == 2 ? lambda{} : block_kfl_lumped_1
    if_kfl_lumped_2 = block_kfl_lumped_2
    
    main_block = lambda{ 
      decl idime,jdime,igaus,inode,jnode,fact,idof,jdof unless @usage == :included
      pr If(@kfl_lumped == 1 => if_kfl_lumped_1, @kfl_lumped == 2 => if_kfl_lumped_2) 
    }      
    
    construct(main_block)
  end
end

class Nest5 < Subroutine
  def initialize(options, functions = nil)
    super("nest5",options,functions)
    
    # ins
    @mnode = Parameters.copy($mnode,:in)
    @pgaus = Parameters.copy($pgaus,:in)
    @pnode = Parameters.copy($pnode,:in)
    @gpvol = Parameters.copy($gpvol,:in)
    @gpsha = Parameters.copy($gpsha,:in)
    @gpcar = Parameters.copy($gpcar,:in)

    # inouts
    @elapu = Parameters.copy($elapu,:inout)
    @elaup = Parameters.copy($elaup,:inout)

    @args.push @mnode,@pgaus,@pnode,@gpvol,@gpsha,@gpcar,@elapu,@elaup
  end

  def generate
    # locals  
    igaus = $igaus
    inode = $inode
    jnode = $jnode
    idime = $idime
    jdime = $jdime
    idof = $idof
    fact = $fact

    for1 = For(igaus,1,@pgaus){
      pr For(inode,1,@pnode){
        pr For(idime,1,@ndime){
          pr idof[idime] === (inode-1)*@ndime + idime
        }.unroll(@unroll)
    
        pr For(jnode,1,@pnode){
          pr fact[1] === @gpvol[igaus] * @gpsha[jnode,igaus]
          pr For(jdime,1,@ndime){
            pr fact[jdime+1] === fact[1] * @gpcar[jdime,inode,igaus]
            pr @elapu[jnode,idof[jdime]] === @elapu[jnode,idof[jdime]] + fact[jdime+1]
          }.unroll(@unroll)
    
          pr For(jdime,1,@ndime){
            pr @elaup[idof[jdime],jnode] === @elaup[idof[jdime],jnode] - fact[jdime+1]
          }.unroll(@unroll)
        }
      }
    }
    
    
    main_block = lambda{
      decl igaus,inode,jnode,idime,jdime,idof,fact unless @usage == :included
      pr for1
    }
    
    construct(main_block)
  end
end

class Nest6 < Subroutine
  def initialize(options, functions = nil)
    super("nest6",options,functions)
    
    # ins
    @kfl_stabi_nsi = Parameters.copy($kfl_stabi_nsi)
    @pgaus = Parameters.copy($pgaus)
    @pnode = Parameters.copy($pnode)
    @gpsp1_p = Parameters.copy($gpsp1_p)
    @wgrgr = Parameters.copy($wgrgr)
    @gpvol = Parameters.copy($gpvol)

    # inouts
    @elapp = Parameters.copy($elapp)

    @args.push @kfl_stabi_nsi,@pgaus,@pnode,@gpsp1_p,@wgrgr,@gpvol,@elapp
  end
  def generate
    # locals
    igaus = $igaus
    inode = $inode
    jnode = $jnode
    fact = $fact
    
    main_block = lambda{
      decl igaus,inode,jnode,fact unless @usage == :included
      pr If(@kfl_stabi_nsi != -1){
        pr For(igaus,1,@pgaus){
          pr For(inode,1,@pnode){
            pr For(jnode,inode+1,@pnode){
              pr fact[2]             === @gpsp1_p[igaus] * @wgrgr[jnode,inode,igaus] * @gpvol[igaus]
              pr @elapp[jnode,inode] === @elapp[jnode,inode] + fact[2]
              pr @elapp[inode,jnode] === @elapp[inode,jnode] + fact[2]
            }
            pr fact[2]             === @gpsp1_p[igaus] * @wgrgr[inode,inode,igaus] * @gpvol[igaus]
            pr @elapp[inode,inode] === @elapp[inode,inode] + fact[2]
          }
        }
      }
    }      

    construct(main_block)
  end
end

class Nest7 < Subroutine
  def initialize(options, functions = nil)
    super("nest7",options,functions)

    # ins
    @pnode = Parameters.copy($pnode)
    @pgaus = Parameters.copy($pgaus)
    @penal_nsi = Parameters.copy($penal_nsi)
    @gpvol = Parameters.copy($gpvol)
    @gpsha = Parameters.copy($gpsha)
    @elpre = Parameters.copy($elpre)

    # inouts
    @elapp = Parameters.copy($elapp)
    @elrbp = Parameters.copy($elrbp)

    @args.push @pnode,@pgaus,@penal_nsi,@gpvol,@gpsha,@elpre,@elapp,@elrbp
  end  
  def generate
    # locals
    fact = $fact
    igaus = $igaus
    inode = $inode

    main_block = lambda {
      decl fact,igaus,inode unless @usage == :included
      pr For(igaus,1,@pgaus){
        pr fact[2] === @penal_nsi * @gpvol[igaus]
        pr For(inode,1,@pnode){
          pr @elapp[inode,inode] === @elapp[inode,inode] + fact[2] * @gpsha[inode,igaus]
          pr @elrbp[inode]       === @elrbp[inode]       + fact[2] * @gpsha[inode,igaus] * @elpre[inode,1]
        }
      }
    }

    construct(main_block)
  end
end

class Nest8 < Subroutine
  def initialize(options, functions = nil)
    super("nest8",options,functions)
    
    # ins
    @pgaus = Parameters.copy($pgaus)
    @pnode = Parameters.copy($pnode)
    @elvel = Parameters.copy($elvel)
    @agrau = Parameters.copy($agrau)
    @gpsp1 = Parameters.copy($gpsp1)
    @kfl_limit_nsi = Parameters.copy($kfl_limit_nsi)
    @kfl_stabi_nsi = Parameters.copy($kfl_stabi_nsi)

    # inouts
    @gpvep = Parameters.copy($gpvep)

    @args.push @pgaus,@pnode,@elvel,@agrau,@gpsp1,@kfl_limit_nsi,@kfl_stabi_nsi,@gpvep
  end
  def generate
    # local
    igaus = $igaus
    idime = $idime
    inode = $inode
    c1 = $c1
    c2 = $c2
    c3 = $c3
    c4 = $c4
    beta = $beta
    alpha = $alpha

    register_funccall("epsilon") if get_lang == FORTRAN
    
    tanh = lambda{|x|
      if get_lang == FORTRAN or @vector_length < 2 then
        return Tanh(x)
      else
        return @functions[:tanh].call(x)
      end
    }
    
    min = lambda{|x,y|
      if get_lang == FORTRAN then 
        return Min( x, y )
      else
        if @vector_length > 1 then
          return @functions[:min].call(x,y)
        else
          return Ternary(x < y, x, y)
        end
      end
    }
              for_igaus = For(igaus,1,@pgaus){
                if get_lang == FORTRAN then
                  pr c1[1] === 0.0.to_var 
                  pr c2[1] === 0.0.to_var 
                  pr c3[1] === 0.0.to_var 
                elsif get_lang == C
                  code =<<EOF
                  memset(c1,0,sizeof(c1));
                  memset(c2,0,sizeof(c2));
                  memset(c3,0,sizeof(c3));
EOF
                  get_output.print code
                end
                pr For(idime,1,@ndime){
                  if get_lang == FORTRAN then
                    pr c4[1] === 0.0.to_var
                  elsif get_lang == C
                    code =<<EOF
                    memset(c4,0,sizeof(c4));
EOF
                    get_output.print code
                  end
      pr For(inode,1,@pnode){
        pr c4[1] === c4[1] + @agrau[inode,igaus] * @elvel[idime,inode,1]
      }
      pr c4[1] === @gpsp1[igaus] * c4[1]
      # Exponential intrinsics only works with ICC
      # pr c1[1] === c1[1] + ( @gpvep[idime,igaus] - c4[1] )**2
      pr c1[1] === c1[1] + ( @gpvep[idime,igaus] - c4[1] ) * ( @gpvep[idime,igaus] - c4[1] )
      pr c3[1] === c3[1] + @gpvep[idime,igaus] * @gpvep[idime,igaus]
      pr c2[1] === c2[1] + c4[1] * c4[1]
    }
    
    pr c3[1]   === Sqrt( c2[1] ) + Sqrt( c3[1] )  
    pr c1[1]   === Sqrt( c1[1] )
                if get_lang == FORTRAN then
                  pr beta[1] === c1[1] / ( c3[1] + epsilon(1.0.to_var) )
                elsif get_lang == C
                  code =<<EOF
                  beta[0] = c1[0] / ( c3[0] + DBL_EPSILON );
EOF
                  get_output.print code
                end
    pr If(@kfl_limit_nsi == 1 => lambda{
            pr alpha[1] === min.call( Set(1.0.to_var, alpha), 2.0.to_var * ( 1.0.to_var  - beta[1] ) )
          }, @kfl_limit_nsi == 2 => lambda{
            pr alpha[1] === 0.5.to_var * ( tanh.call( 20.0.to_var * ( beta[1] - 0.8.to_var ) ) + 1.0.to_var )
          })              
                pr For(idime,1,@ndime){
                  pr @gpvep[idime,igaus] === alpha[1] * @gpvep[idime,igaus]
                }
              }

              if get_lang == C then
                code =<<EOF
                memset(gpvep,0,sizeof(#{@gpvep.type.decl}) * #{@ndime} * pgaus);
EOF
              end

              main_block = lambda{
                decl igaus,idime,inode,c1,c2,c3,c4,beta,alpha unless @usage == :included

                pr If( @kfl_stabi_nsi == -1 => lambda{
                         if get_lang == FORTRAN then
                           pr @gpvep === 0.0
                         elsif get_lang == C
                           get_output.print code
                         end           
                       }, @kfl_limit_nsi == -1 => lambda{
                         if get_lang == FORTRAN then
                           pr @gpvep === 0.0
                         elsif get_lang == C
                           get_output.print code
                         end           
                       }, @kfl_limit_nsi > 0 => lambda{ pr for_igaus } )
              }

    construct(main_block)
  end
end

class Nest9 < Subroutine
  def initialize(options, functions = nil)
    super("nest9",options,functions)
    
    # ins
    @pgaus = Parameters.copy($pgaus)
    @gpsp1_p = Parameters.copy($gpsp1_p)
    @gprhs = Parameters.copy($gprhs)
    @gpden = Parameters.copy($gpden)
    @dtsgs = Parameters.copy($dtsgs)
    @gpsp1_v = Parameters.copy($gpsp1_v)
    @gpsgs = Parameters.copy($gpsgs)
    @kfl_sgsti_nsi = Parameters.copy($kfl_sgsti_nsi)
    @kfl_stabi_nsi = Parameters.copy($kfl_stabi_nsi)

    # inouts
    @gpgrp = Parameters.copy($gpgrp)
    @gpvep = Parameters.copy($gpvep)

    @args.push @pgaus,@gpsp1_p,@gprhs,@gpden,@dtsgs,@gpsp1_v,@gpsgs,@kfl_sgsti_nsi,@kfl_stabi_nsi,@gpgrp,@gpvep
  end
  def generate
    # local
    igaus = $igaus
    idime = $idime
    fact = $fact
    fact1_p = $fact1_p

              for_igaus = []
              for_igaus[0] = For(igaus,1,@pgaus){
                pr For(idime,1,@ndime){
                  pr @gpgrp[idime,igaus] === @gpgrp[idime,igaus] + @gpsp1_p[igaus] * @gprhs[idime,igaus]
                }.unroll(@unroll)
              }
              
              for_igaus[1] = For(igaus,1,@pgaus){
                pr fact[2]     === @gpden[igaus] * @dtsgs[1] * @gpsp1_v[igaus]
                pr fact1_p[1]  === @gpden[igaus] * @dtsgs[1] * @gpsp1_p[igaus]
                pr For(idime,1,@ndime){
                  pr @gpvep[idime,igaus] === @gpvep[idime,igaus] + fact[2]    * @gpsgs[idime,igaus,2]
                  pr @gpgrp[idime,igaus] === @gpgrp[idime,igaus] + fact1_p[1] * @gpsgs[idime,igaus,2]
                }.unroll(@unroll)
              }

              block = lambda{
                pr for_igaus[0]
                pr If(@kfl_sgsti_nsi == 1){
                  pr for_igaus[1]
                }
              }

              main_block = lambda{
                decl igaus,idime,fact,fact1_p unless @usage == :included
                pr If(@kfl_stabi_nsi == -1 => lambda{ 
                        if get_lang == FORTRAN then
                          pr @gpgrp === 0.0 
                        elsif get_lang == C
                          code =<<EOF
                          memset(gpgrp,0,sizeof(#{@gpgrp.type.decl}) * #{@ndime} * pgaus);
EOF
                        end
                        get_output.print code
                      }, :else => block)
              }
      
    
    construct(main_block)
  end
end

class Nest10 < Subroutine
  def initialize(options, functions = nil)
    super("nest10",options,functions)
    
    # ins
    @mnode = Parameters.copy($mnode)
    @pgaus = Parameters.copy($pgaus)
    @gpden = Parameters.copy($gpden)
    @dtinv_mod = Parameters.copy($dtinv_mod)
    @nbdfp_nsi = Parameters.copy($nbdfp_nsi)
    @gpvel = Parameters.copy($gpvel)
    @pnode = Parameters.copy($pnode)
    @gpvol = Parameters.copy($gpvol)
    @gpsha = Parameters.copy($gpsha)
    @agrau = Parameters.copy($agrau)
    @gprhc = Parameters.copy($gprhc)
    @gpcar = Parameters.copy($gpcar)
    @gpvep = Parameters.copy($gpvep)
    @gpgrp = Parameters.copy($gpgrp)

    # inouts
    @gprhs = Parameters.copy($gprhs)
    @elrbu = Parameters.copy($elrbu)
    @elrbp = Parameters.copy($elrbp)

    @args.push @mnode,@pgaus,@gpden,@dtinv_mod,@nbdfp_nsi,@gpvel,@pnode,@gpvol,@gpsha,@agrau,@gprhc,@gpcar,@gpvep,@gpgrp,@gprhs,@elrbu,@elrbp
  end

  def generate
    # locals
    igaus = $igaus
    fact = $fact
    idime = $idime
    itime = $itime
    inode = $inode

      for_igaus = For(igaus,1,@pgaus){
        pr fact[5] === @gpden[igaus] * @dtinv_mod[1]
        pr For(itime,2,@nbdfp_nsi){
          pr For(idime,1,@ndime){
            pr @gprhs[idime,igaus] === @gprhs[idime,igaus] - @functions[:pabdf_nsi].call(itime) * fact[5] * @gpvel[idime,igaus,itime]
          }.unroll(@unroll)
        }
        pr For(inode,1,@pnode){
          pr fact[2] === @gpvol[igaus] * @gpsha[inode,igaus]  # ( f + rho*u^n/dt , v )
          pr fact[4] === @gpvol[igaus] * @agrau[inode,igaus]  # ( rho * a.grad(v) , P1' ) 
          pr For(idime,1,@ndime){
            pr @elrbu[idime,inode] === @elrbu[idime,inode] + fact[2] * @gprhs[idime,igaus] + fact[4] * @gpvep[idime,igaus]
          }.unroll(@unroll)
          pr @elrbp[inode] === @elrbp[inode] + @gpvol[igaus] * @gpsha[inode,igaus] * @gprhc[igaus]  # ( rhs, q )
          pr For(idime,1,@ndime){
            pr @elrbp[inode] === @elrbp[inode] + @gpvol[igaus] * @gpcar[idime,inode,igaus] * @gpgrp[idime,igaus] # ( P2' , grad(q) ) 
          }.unroll(@unroll)
        }
      }
    
    main_block = lambda{
      decl igaus,fact,idime,itime,inode unless @usage == :included
      pr for_igaus
    }

    construct(main_block)
    end
end

class Nest11 < Subroutine
            def initialize(options, functions = nil)
              super("nest11",options,functions)
              
              # ins
              @mnode = Parameters.copy($mnode)
              @pgaus = Parameters.copy($pgaus)
              @pnode = Parameters.copy($pnode)
              @gpcar = Parameters.copy($gpcar)
              @gpvol = Parameters.copy($gpvol)
              @gpsha = Parameters.copy($gpsha)
              @gpcar_bub = Parameters.copy($gpcar_bub)
              @pbubl = Parameters.copy($pbubl)
              @kfl_stabi_nsi = Parameters.copy($kfl_stabi_nsi)
              @kfl_press_nsi = Parameters.copy($kfl_press_nsi)
              @penal_nsi = Parameters.copy($penal_nsi)
              @elbub = Parameters.copy($elbub)
              @gprhc = Parameters.copy($gprhc)
              @gpsha_bub = Parameters.copy($gpsha_bub)

              # inouts
              @elauq = Parameters.copy($elauq)
              @elaqu = Parameters.copy($elaqu)
              @elapq = Parameters.copy($elapq)
              @elaqp = Parameters.copy($elaqp)
              @elaqq = Parameters.copy($elaqq)
              @elrbq = Parameters.copy($elrbq)

              @args.push @mnode,@pgaus,@pnode,@gpcar,@gpvol,@gpsha,@gpcar_bub,@pbubl,@kfl_stabi_nsi,@kfl_press_nsi,@penal_nsi,@elbub,@gprhc,@gpsha_bub,@elauq,@elaqu,@elapq,@elaqp,@elaqq,@elrbq
            end

            def generate
              # locals
              fact = $fact
              igaus  = $igaus
              idof = $idof
              idime = $idime
              inode = $inode

              register_funccall("stop") if get_lang == FORTRAN
              register_funccall("maxval") if get_lang == FORTRAN

              exp_elauq = []
              exp_elauq[0] = "@elauq[idof[@ndime],1] === @elauq[idof[@ndime],1] - fact[2] * @gpcar[idime,inode,igaus]"
              exp_elauq[1] = "@elauq[idof[@ndime],1] === @elauq[idof[@ndime],1] + @gpvol[igaus] * @gpsha[inode,igaus] * @gpcar_bub[idime,igaus]"
              
              block_igaus = lambda{|exp|
                return For(igaus,1,@pgaus){
                  pr fact[2] === @gpvol[igaus] * @gpsha_bub[igaus]
                  pr For(inode,1,@pnode){
                    pr For(idime,1,@ndime){
                      pr idof[@ndime] === (inode-1)*@ndime + idime
                      pr eval exp
                      pr @elaqu[1,idof[@ndime]] === @elaqu[1,idof[@ndime]] + fact[2] * @gpcar[idime,inode,igaus]
                    }.unroll(@unroll)
                  }
                }
              }
              
              if_kfl_press_nsi_1 = block_igaus.call(exp_elauq[0])
              if_kfl_press_nsi_else = block_igaus.call(exp_elauq[1])

              call_maxval = lambda{|x|
                if get_lang == FORTRAN
                  return maxval(x)
                else
                  return @functions[:maxval].call(x)
                end
              }

              main_block = lambda{
                decl fact,igaus,idof,idime,inode unless @usage == :included
                pr If(call_maxval.call(@pbubl) == 1){
                  pr If(@kfl_stabi_nsi != -1){
                    if get_lang == FORTRAN then
                      stop
                    end
                  }
                  
                  if get_lang == FORTRAN then
                    pr @elauq === 0.0
                    pr @elapq === 0.0
                    pr @elaqu === 0.0
                    pr @elaqp === 0.0
                    pr @elaqq === 0.0
                    pr @elrbq === 0.0
                  elsif get_lang == C
                    code =<<EOF
                    memset(elauq, 0, sizeof(#{@elauq.type.decl}) * pnode * #{@ndime});
                    memset(elapq, 0, sizeof(#{@elapq.type.decl}) * pnode);
                    memset(elaqu, 0, sizeof(#{@elaqu.type.decl}) * pnode * #{@ndime});
                    memset(elaqp, 0, sizeof(#{@elaqp.type.decl}) * pnode);
                    memset(elaqq, 0, sizeof(#{@elaqq.type.decl}));
                    memset(elrbq, 0, sizeof(#{@elrbq.type.decl}));
EOF
                    get_output.print code
                  end
                  pr If(@kfl_press_nsi == 1 => lambda{ pr if_kfl_press_nsi_1}, :else => lambda{ pr if_kfl_press_nsi_else })
                  
                  pr For(igaus,1,@pgaus){
                    pr @elaqq[1,1] === @elaqq[1,1] + @gpvol[igaus] * @gpsha_bub[igaus] * @penal_nsi
                    pr @elrbq[1]   === @elrbq[1]   + @gpvol[igaus] * @gpsha_bub[igaus] * @penal_nsi * @elbub[1]
                    pr @elrbq[1]   === @elrbq[1]   + @gpvol[igaus] * @gpsha_bub[igaus] * @gprhc[igaus]
                  }
                }
              }
              
              construct(main_block)
            end
          end

class KSplitBoast < KSplitOss
  def generate
    register_funccall("maxval") if get_lang == FORTRAN
    register_funccall("epsilon") if get_lang == FORTRAN
    register_funccall("stop") if get_lang == FORTRAN

    x1 = Int("x", :dir => :in)
    y1 = Real("y")
    @pabdf_nsi = Procedure("pabdf_nsi",[x1], :return => y1){
      pr y1 === 1.0
    }          

    unless get_lang == FORTRAN then
      x2 = Int("x", :dir => :in, :vector_length => @opts[:vector_length], :dim => [Dim(1)])
      y2 = Int("y")
      @p_maxval = Procedure("maxval",[x2], :return => y2){
        if @opts[:vector_length] > 1 then
          decl i = Int("i")
          decl a = Int("a", :dim => [Dim(@opts[:vector_length])], :allocate => true)
    
          pr a[1] === x2.dereference
          pr y2 === a[1]
    
          pr For(i,2,@opts[:vector_length]){
            pr If(a[i] > y2){
              pr y2 === a[i]
            }
          }
        else
          pr y2 === x2[1]
        end
      }            
    end          
    unless get_lang == FORTRAN then
      x3 = Real("x", :dir => :in, :vector_length => @opts[:vector_length])
      y3 = Real("y", :vector_length => @opts[:vector_length])
    
      @p_tanh = Procedure("Tanh", [x3], :return => y3){
        decl i = Int("i")
        decl tmp1 = Real("tmp1", :dim => [Dim(@opts[:vector_length])], :allocate => true)
        decl tmp2 = Real("tmp2", :dim => [Dim(@opts[:vector_length])], :allocate => true)
        
        pr tmp1[1] === x3
    
        pr For(i,1,@opts[:vector_length]){
          pr tmp2[i] === Tanh(tmp1[i])
        }
        pr y3 === tmp2[1] 
      }
    end
    unless get_lang == FORTRAN then
      x4 = Real("x", :dir => :in, :vector_length => @opts[:vector_length])
      y4 = Real("y", :vector_length => @opts[:vector_length])
    
      @p_sqrt = Procedure("Sqrt", [x4], :return => y4){
        decl i = Int("i")
        decl tmp1 = Real("tmp1", :dim => [Dim(@opts[:vector_length])], :allocate => true)
        decl tmp2 = Real("tmp2", :dim => [Dim(@opts[:vector_length])], :allocate => true)
        
        pr tmp1[1] === x4
    
        pr For(i,1,@opts[:vector_length]){
          pr tmp2[i] === Sqrt(tmp1[i])
        }
        pr y4 === tmp2[1] 
      }
    end
    unless get_lang == FORTRAN then
      x5 = Real("x", :dir => :in, :vector_length => @opts[:vector_length])
      y5 = Real("y", :dir => :in, :vector_length => @opts[:vector_length])
      z5 = Real("z", :vector_length => @opts[:vector_length])
    
      @p_min = Procedure("Min", [x5,y5], :return => z5){
        decl i = Int("i")
        decl tmp1 = Real("tmp1", :dim => [Dim(@opts[:vector_length])], :allocate => true)
        decl tmp2 = Real("tmp2", :dim => [Dim(@opts[:vector_length])], :allocate => true)
        decl tmp3 = Real("tmp3", :dim => [Dim(@opts[:vector_length])], :allocate => true)
        
        pr tmp1[1] === x5
        pr tmp2[1] === y5
    
        pr For(i,1,@opts[:vector_length]){
          pr tmp3[i] === Ternary(tmp1[i] < tmp2[i], tmp1[i], tmp2[i])
        }
        pr z5 === tmp3[1] 
      }
    end              

    @opts[:ndime] = @ndime
    nests = {}
    @opts[:nests].each{|i|
      case i
      when 2,10
        functions = {:pabdf_nsi => @pabdf_nsi}
      when 8
        functions = get_lang == FORTRAN ? nil :{:tanh => @p_tanh, :min => @p_min}
      when 11
        functions = get_lang == FORTRAN ? nil : {:maxval => @p_maxval}
      else
        functions = nil
      end
      eval "nests[:nest#{i}] = Nest#{i}::new(@opts,functions)"
      eval "nests[:nest#{i}].generate"
    }
  
    includes = ["immintrin.h"]
    includes.push("string.h", "math.h", "float.h") if get_lang == C
    @kernel = CKernel::new( :includes => includes )
    functions = [@pabdf_nsi]
    functions.push @p_maxval unless get_lang == FORTRAN
    @kernel.procedure = declare_procedure("nsi_element_assembly_split_oss_ndime#{@ndime}",functions.flatten)

    # Printing subroutines           
    pr @pabdf_nsi
    unless get_lang == FORTRAN then
      pr @p_maxval
      pr @p_min
      pr @p_tanh
    end
    nests.values.each{|n| pr n.code } unless @opts[:inline] == :included

    # Main procedure body
    opn @kernel.procedure
      # Local arrays
      # Local arrays
      @gpsp1_p   = $gpsp1_p
      @gpsp1_v   = $gpsp1_v
      @gpsp2_v   = $gpsp2_v
      @c1        = $c1
      @c2        = $c2
      @c3        = $c3
      @c4        = $c4
      @alpha     = $alpha
      @beta      = $beta
      
      @gpveo     = $gpveo
      @fact1_p   = $fact1_p
      @dtinv_mod = $dtinv_mod
      @inode     = $inode
      @jnode     = $jnode
      @jdime     = $jdime
      @idofv     = $idofv
      @ivect     = $ivect
      @igaus     = $igaus
      @idime     = $idime
      @jdofv     = $jdofv
      @itime     = $itime
      
      @fact = $fact
      
      @idof = $idof
      @jdof = $jdof
      
      decl @gpsp1_p,@gpsp1_v,@gpsp2_v,@c1,@c2,@c3,@c4,@alpha,@beta,
           @gpveo,@fact1_p,@dtinv_mod,@inode,@jnode,@jdime,@idofv,
           @ivect,@igaus,@idime,@jdofv,@itime,@fact,@idof,@jdof 
    
if get_lang == FORTRAN
  pr @dtinv_mod === @dtinv_loc
  pr @gpsp1_p === @gpsp1
  pr @gpsp1_v === @gpsp1
  pr @gpsp2_v === @gpsp2

  pr If( @kfl_nota1_nsi == 1 => lambda{ pr @gpsp1_v === 0.0 } )

  pr @elrbp === 0.0
  pr @elrbu === 0.0
  pr @elapp === 0.0
  pr @elauu === 0.0
  pr @elaup === 0.0
  pr @elapu === 0.0
elsif get_lang == C
  code =<<EOF
  memcpy(dtinv_mod, dtinv_loc, sizeof(dtinv_mod));
  memcpy(gpsp1_p, gpsp1, sizeof(gpsp1_p));
  memcpy(gpsp1_v, gpsp1, sizeof(gpsp1_v));
  memcpy(gpsp2_v, gpsp2, sizeof(gpsp2_v));

  if (kfl_nota1_nsi == 1) memset(gpsp1_v, 0, sizeof(gpsp1_v));

  memset(elrbp, 0, sizeof(#{@elrbp.type.decl}) * pnode);
  memset(elrbu, 0, sizeof(#{@elrbu.type.decl}) * #{@ndime} * pnode);
  memset(elapp, 0, sizeof(#{@elapp.type.decl}) * pnode * pnode);
  memset(elauu, 0, sizeof(#{@elauu.type.decl}) * pnode * #{@ndime} * pnode * #{@ndime});
  memset(elaup, 0, sizeof(#{@elaup.type.decl}) * pnode * #{@ndime} * pnode);
  memset(elapu, 0, sizeof(#{@elapu.type.decl}) * pnode * pnode * #{@ndime});
EOF
  get_output.print code
end

      # Either call subroutine or paste their code
      nests.values.each{|n| n.call} 

    close @kernel.procedure
  end        
end
