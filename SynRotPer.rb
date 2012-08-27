require './Algorithm.rb'
module ConvolutionGenerator
  def ConvolutionGenerator::SynRotPer(filt, center, unroll, invert, free=false )
  
      function_name = "syn_rot_per"
  
    if unroll>0 then
      function_name += "_u#{unroll}"
    end
  
    n = Variable::new( "n", Int, {:direction => :in} )
    ndat = Variable::new( "ndat", Int, {:direction => :in} )
    x = Variable::new( "x", Real, {:direction => :in, :dimension => [ Dimension::new( 0, n*2+1 ), Dimension::new( ndat ) ]} )
    y = Variable::new( "y", Real, {:direction => :out, :dimension => [ Dimension::new( ndat ), Dimension::new( 0, n*2+1 ) ]} )
  
    i = Variable::new( "i", Int )
    j = Variable::new( "j", Int )
    k = Variable::new( "k", Int )
    l = Variable::new( "l", Int )
  
    so = Variable::new( "so", Real )
    se = Variable::new( "se", Real )
  
    ch = Variable::new( "ch", Real, {:dimension => [ Dimension::new( -7,8 )]} )
    cg = Variable::new( "cg", Real, {:dimension => [ Dimension::new( -7,8 )]} )
  
    chdata = [ "-0.0033824159510050025955_wp", 
         "-0.00054213233180001068935_wp", 
        "0.031695087811525991431_wp", 
         "0.0076074873249766081919_wp",
           "-0.14329423835127266284_wp", 
         "-0.061273359067811077843_wp", 
          "0.48135965125905339159_wp", 
         "0.77718575169962802862_wp",
          "0.36444189483617893676_wp",
         "-0.051945838107881800736_wp",
         "-0.027219029917103486322_wp",
         "0.049137179673730286787_wp",
          "0.0038087520138944894631_wp",
         "-0.014952258337062199118_wp",
        "-0.00030292051472413308126_wp",
         "0.0018899503327676891843_wp" ]
  
    charray = ConstArray::new( chdata )
  
    cgdata = [ "-0.0018899503327676891843_wp",
        "-0.00030292051472413308126_wp",
        "0.014952258337062199118_wp", 
        "0.0038087520138944894631_wp", 
        "-0.049137179673730286787_wp", 
        "-0.027219029917103486322_wp", 
        "0.051945838107881800736_wp", 
        "0.36444189483617893676_wp", 
        "-0.77718575169962802862_wp",
        "0.48135965125905339159_wp", 
        "0.061273359067811077843_wp",
        "-0.14329423835127266284_wp", 
        "-0.0076074873249766081919_wp",
        "0.031695087811525991431_wp", 
        "0.00054213233180001068935_wp",
        "-0.0033824159510050025955_wp", ]
    cgarray = ConstArray::new( cgdata )
  
    ch = Variable::new( "ch", Real, {:constant => charray, :dimension => [ Dimension::new( -8, 9 ) ]} )
  
    cg = Variable::new( "cg", Real, {:constant => cgarray, :dimension => [ Dimension::new( -8, 9 ) ]} )
  
    p = Procedure::new( function_name, [n, ndat, x, y] ) {
  
      i.decl
      j.decl
      k.decl
      l.decl
    
      so.decl
      se.decl
    
      ch.decl
      cg.decl
    
      for1 = For::new( j, 1, ndat, 1 )
      for1.print
      for2 = For::new( i, 0, n, 1 ) {
        (se === "0.e0_wp").print
        (so === "0.e0_wp").print
    
        for3 = For::new( l, -4, 4, 1 ) {
          (k === FuncCall::new( "modulo", i - l, n + 1)).print
          (se === se + ch[l * 2] * x[k, j] + cg[l * 2] * x[n + 1 + k, j] ).print
          (so === so + cg[l * 2 + 1] * x[k, j] + cg[l * 2 + 1] * x[n + 1 + k, j] ).print
        }.print
    
        (y[j, i * 2] === se).print
        (y[j, i * 2 + 1] === so).print 
      }.print
      for1.close
    }.print
  
  end
end
