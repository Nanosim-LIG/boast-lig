require 'yaml'

f = File::open('intrinsicsX86.yaml')
hash = YAML::load(f.read)

h = {}

hash["intrinsics_list"]["intrinsic"].each { |intrinsic|
  name = intrinsic.delete("name")
  h[name] = intrinsic["CPUID"]
}

puts YAML::dump(h)
