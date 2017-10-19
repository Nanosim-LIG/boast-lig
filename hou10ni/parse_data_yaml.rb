require 'yaml'
require 'pp'
require 'csv'


# USAGE : ruby parse_data_yaml.rb ../experiments/../Data.yaml test.cvs
input = ARGV[0]
output = ARGV[1]

#header = [:kernel,:optim_main, :optim_nested, :omp_num_threads,:time]
header = []
body = []

h = YAML::load(File::open(input).read)
h.each{|k1,v1|
  #v1[:time].each{|e|
  	e = v1[:time].min
    t = []
    t.push k1[:kernel]
    t.push k1[:optim_main]
    t.push k1[:optim_nested]
    t.push k1[:omp_num_threads]
		t.push k1[:FCFLAGS]
    t.push(e)
    body.push(t)
  #}
}

CSV.open(output, "w"){ |f|
  f << header
  body.each{ |e|
    f << e
  }
}
