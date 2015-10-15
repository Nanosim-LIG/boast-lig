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
  intrinsic["parameter"] = [intrinsic["parameter"]].flatten if intrinsic["parameter"]
  intrinsic["type"] = [intrinsic["type"]].flatten if intrinsic["type"]
  intrinsic["category"] = [intrinsic["category"]].flatten if intrinsic["category"]
}
puts YAML::dump hash
