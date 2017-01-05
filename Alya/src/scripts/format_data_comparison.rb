require 'yaml'
require 'pp'
require 'csv'

input = ARGV[0]
output = ARGV[1]

header = [:kernel,:nest,:time]
body = []

h = YAML::load(File::open(input).read)
h.each{|k1,v1|
  v1.each{|e|
    t = []
    k1.each{|k2,v2|
      t.push(v2)
    }
    t.push(e)
    body.push(t)
  }
}

CSV.open(output, "w"){ |f|
  f << header
  body.each{ |e|
    f << e
  }
}

