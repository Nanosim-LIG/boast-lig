require 'yaml'
require 'active_support/core_ext/hash/conversions'
f = File::open('intrinsicsX86.yaml')
hash = YAML::load(f.read)
hash["intrinsics_list"]["intrinsic"].each { |intrinsic|
#  puts intrinsic["name"]
  list = []
  if intrinsic["CPUID"] then
    [intrinsic["CPUID"]].flatten.each { |cpuid|
      list += cpuid.split("/")
    }
    intrinsic["CPUID"]=list
  end
  list = []
  if intrinsic["parameter"] then
#    puts intrinsic["parameter"].inspect
    intrinsic["parameter"] = [intrinsic["parameter"]].flatten
  end
} 
puts YAML::dump hash
