require 'BOAST'
include BOAST

set_lang(CL)

register_funccall("MAX")
register_funccall("MIN")

set_array_start(0)

height = Int("height")
width = Int("width")
pdst = Int("pdst", :dir => :in,  :signed => false, :size => 1, :dim => [ Dim(3), Dim(width), Dim(height) ] )
psrc = Int("psrc", :dir => :out, :signed => false, :size => 1, :dim => [ Dim(3), Dim(width), Dim(height) ] )

get_output.puts <<EOF
#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))
EOF
p = Procedure("math", [pdst, psrc, width, height]) {
  decl y = Int("y")
  decl x = Int("x")
  decl xBoundary = Int("xBoundary")
  decl yBoundary = Int("yBoundary")
  pr xBoundary === width - 2
  pr yBoundary === height - 2
  
  pr y === get_global_id(0)
  pr x === get_global_id(1)

  pr If( Expression( "||", x >= xBoundary, y >= yBoundary ), lambda {
    (0..2).each { |ind|
      pr pdst[ind, x, y] === psrc[ind, x, y]
    }
  }, lambda {
    acolors = ["b", "g", "r"].collect { |c|
      Int(c+"Color")
    }
    decl *acolors
    acolors.each { |c|
      pr c === 0
    }
    (0..2).each { |ind|
      pr acolors[ind] === acolors[ind] - psrc[ind, x, y]     - psrc[ind, x + 1, y]         - psrc[ind, x + 2, y] \
                                       - psrc[ind, x, y + 1] + psrc[ind, x + 1, y + 1] * 9 - psrc[ind, x + 2, y + 1] \
                                       - psrc[ind, x, y + 2] - psrc[ind, x + 1, y + 2]     - psrc[ind, x + 2, y + 2]
    }
    colors = ["blue", "green", "red"].collect { |c|
      Int(c, :size => 1, :signed => false)
    }
    decl *colors
    (0..2).each { |ind|
      pr colors[ind] === MAX(MIN(acolors[ind], 255), 0)
    }
    (0..2).each { |ind|
      pr pdst[ind, x, y] === colors[ind]
    }
  })

}

pr p
