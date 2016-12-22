require 'yaml'
require 'pp'
require 'csv'

input = ARGV[0]
output = ARGV[1]

header = [:vector_length,:preprocessor,:nest,:unroll,:usage,:CFLAGS,:pgaus,:pnode,:time]
body = []

h = YAML::load(File::open(input).read)
h.each{|k1,v1|
  v1[:time].each{|e1|
    t = []
    k1.each{|k2,v2|
      t.push v2.class == Array ? v2[0] : v2
    }
    t.push v1[:characteristics][:pgaus]
    t.push v1[:characteristics][:pnode]
    t.push(e1)
    body.push(t)
  }
}

CSV.open(output, "w"){ |f|
  f << header
  body.each{ |e|
    f << e
  }
}

