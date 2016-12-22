require_relative "./KASGS.rb"

class KASGSBoast < KASGS
  def set_arr(arr,value)
    if get_lang == FORTRAN then
      pr arr === value
    else
      get_output.print <<EOF
      memset(#{arr}, 0, sizeof(#{arr.type.decl}) * #{arr.dimension.join("*")});
EOF
    end
  end

  def copy_arr(dst,src)
    if get_lang == FORTRAN then
      pr dst === src
    else
      get_output.print <<EOF
      memcpy(#{dst}, #{src}, sizeof(#{dst.type.decl}) * #{dst.dimension.join("*")});
EOF
    end
  end

  def decl_max

  end

  def decl_local
    @inode = $inode
    @jnode = $jnode
    @idofv = $idofv
    @idof = $idof
    @jdime = $jdime
    @ivect = $ivect
    @igaus = $igaus
    @idime = $idime
    @jdofv = $jdofv
    @kdime = $kdime
    @xvis2 = $xvis2
    @xvisc = $xvisc
    @one_rp = $one_rp
    @p1ve2 = $p1ve2
    @p1vec = $p1vec
    @p2vec = $p2vec
    @p2sca = $p2sca
    @wgrvi = $wgrvi
    @fact = $fact
    @factx = $factx
    @facx1 = $facx1
    @ugraN = $ugraN
    @gramugraN = $gramugraN
    @penal = $penal
    @gprhh = $gprhh
    @taupr = $taupr
    @gpveo = $gpveo
    @gpcar1ji = $gpcar1ji
    @gpcar2ji = $gpcar2ji
    @gpcar3ji = $gpcar3ji
    @p2sca_bub = $p2sca_bub
    @p2vec_bub = $p2vec_bub
    @elauu_a = $elauu_a
    @factvec = $factvec

    decl @inode,@jnode,@idofv,@idof,@jdime,@ivect,@igaus,@idime,
         @jdofv,@kdime,@xvis2,@xvisc,@one_rp,@p1ve2,@p1vec,@p2vec,
         @p2sca,@wgrvi,@fact,@factx,@facx1,@ugraN,@gramugraN,
         @penal,@gprhh,@taupr,@gpveo,@gpcar1ji,@gpcar2ji,@gpcar3ji,
         @p2sca_bub,@p2vec_bub
    
    @elauu_a.each{|e1|
      e1.each{|e2|
        decl e2
      }
    }

    @factvec.each{|e|
      decl e
    }
  end

  def generate
    @openacc = false

    includes = ["immintrin.h"]
    includes.push("string.h", "math.h", "float.h") if get_lang == C
    @kernel = CKernel::new( :includes => includes )
    @kernel.procedure = declare_procedure

    opn @kernel.procedure
      decl_local
      #----------------------------------------------------------------------
      #
      # Local variables
      #
      #----------------------------------------------------------------------
      #
      # Viscous term factor
      #
      pr @xvis2 === 0.0.to_var
      pr If( @fvins_nsi > 0.9.to_var => lambda{
               pr @xvisc === 1.0.to_var
               pr If( @fvins_nsi > 1.9.to_var) {pr @xvis2 === 2.0.to_var/3.0.to_var}
               
             },
             :else => lambda{pr @xvisc === 0.0.to_var} )
      #
      # Do not stabilize convection, reaction, coriolis
      #
      pr If( @kfl_nota1_nsi == 1 => lambda{pr @one_rp === 0.0.to_var}, 
             :else => lambda{pr @one_rp === 1.0.to_var})

      for_openacc = For($ivect,1,@opts[:vector_length]) if @openacc

      opn for_openacc if @openacc
        # initialization
          set_arr(@elauu, 0.0.to_var)
          set_arr(@elaup, 0.0.to_var)
          set_arr(@elapu, 0.0.to_var)
          set_arr(@elapp, 0.0.to_var)
          set_arr(@elrbu, 0.0.to_var)
          set_arr(@elrbp, 0.0.to_var)

          3.times{|i|
            3.times{|j|
              set_arr(@elauu_a[i][j], 0.0.to_var)
            }
          }

        #
        # Incompressible vs low Mach
        #

        # Pb writting a max equivalent with BOAST that can 
        # return Intrinsics
        pr If( @kfl_regim_nsi == 3 => lambda{
                 set(@taupr, 1.0.to_var)
               },
               :else => lambda{
                 # set(@taupr, 1.0.to_var / max(@gpst1,zeror))
               })
      
      close for_openacc if @openacc
    close @kernel.procedure
  end        
end
