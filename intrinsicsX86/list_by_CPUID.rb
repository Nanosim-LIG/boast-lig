require 'nokogiri'
require 'yaml'

intrinsics_hash = Hash::new { |h,k| h[k] = [] }

doc = File.open('data-3.3.12.xml') { |f| Nokogiri::XML(f) }
intrinsics = doc.xpath("//intrinsic")
intrinsics.each { |intrinsic|
  intrinsic.xpath("CPUID").each { |cpuid|
    intrinsics_hash[cpuid.child.to_s].push(intrinsic['name'].to_s)
  }
}
puts YAML::dump(intrinsics_hash)
