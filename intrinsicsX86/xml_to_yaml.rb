require 'yaml'
require 'active_support/core_ext/hash/conversions'
f = File::open('data-3.3.12.xml')
hash = Hash.from_xml(f.read)
hash["intrinsics_list"]["intrinsic"].each { |intrinsic|
  list = []
  if intrinsic["CPUID"] then
    [intrinsic["CPUID"]].flatten.each { |cpuid|
      list += cpuid.split("/")
    }
    intrinsic["CPUID"]=list
  end
  list = []
  if intrinsic["parameter"] then
    intrinsic["parameter"] = [intrinsic["parameter"]].flatten
  end
}
puts YAML::dump hash
