require 'BOAST'
include BOAST
require 'narray'

class Numeric
  def clamp( min, max )
    [[self, max].min, min].max
  end
end

def rndup( val, div)
  return (val%div) == 0 ? val : val + div - (val%div)
end

set_array_start(0)
register_funccall("clamp")

def laplacian_ref( input, output, width, height )
  (1...(height-1)).each { |y_i|
    (3...(3*width-3)).each { |x_i|
      output[x_i, y_i] = ( -input[x_i - 3, y_i - 1] - input[x_i, y_i - 1] - input[x_i + 3, y_i - 1]\
                          - input[x_i - 3, y_i]     + input[x_i, y_i] * 9 - input[x_i + 3, y_i]\
                          - input[x_i - 3, y_i + 1] - input[x_i, y_i + 1] - input[x_i + 3, y_i + 1] ).clamp(0,255)
    }
  }
end

def laplacian_c_ref
  k = CKernel::new
  height = Int("height", :dir => :in)
  width = Int("width", :dir => :in)
  w = Int("w")
  pdst = Int("pdst", :dir => :out,  :signed => false, :size => 1, :dim => [ Dim(w), Dim(height) ] )
  psrc = Int("psrc", :dir => :in, :signed => false, :size => 1, :dim => [ Dim(w), Dim(height) ] )
  p = Procedure("math", [psrc, pdst, width, height]) {
    decl i = Int("i")
    decl j = Int("j")
    decl tmp = Int("tmp")
    decl w
    pr w === width * 3
    pr For(j, 1, height-2) {
      pr For(i, 3, w-4) {
         
        pr tmp === ( -psrc[i - 3, j - 1] - psrc[i, j - 1] - psrc[i + 3, j - 1]\
                    - psrc[i - 3, j]     + psrc[i, j] * 9 - psrc[i + 3, j]\
                    - psrc[i - 3, j + 1] - psrc[i, j + 1] - psrc[i + 3, j + 1] )
        pr pdst[i,j] === Ternary(tmp < 0, 0, Ternary(tmp>255, 255, tmp))
      }
    }
  }
  pr p
  k.procedure = p
  return k
end

def laplacian(x_component_number = 1, vector_length=1, y_component_number = 1, temporary_size = 4, vector_recompute = false, load_overlap = false)
 
  k = CKernel::new
  height = Int("height", :dir => :in)
  width = Int("width", :dir => :in)
  w = Int("w")
  pdst = Int("pdst", :dir => :out,  :signed => false, :size => 1, :dim => [ Dim(w), Dim(height) ] )
  psrc = Int("psrc", :dir => :in, :signed => false, :size => 1, :dim => [ Dim(w), Dim(height) ] )
  
  get_output.puts <<EOF
#define vload1(offset, p) *(p + offset)
#define vstore1(d, offset, p) *(p + offset) = d
EOF
  p = Procedure("math", [psrc, pdst, width, height]) {
    decl y = Int("y")
    decl x = Int("x")
    decl w
    
    pr x === get_global_id(0) * x_component_number 
    pr y === get_global_id(1) * y_component_number
    pr w === width * 3

    vector_number = (x_component_number.to_f/vector_length).ceil
    vector_number = vector_number < 1 ?  1 : vector_number
    total_x_size = vector_recompute ? vector_number * vector_length : x_component_number

    x_offset = total_x_size + 3
    y_offset = y_component_number + 1

    pr x === Ternary(x < 3, 3, Ternary( x > w      - x_offset, w      - x_offset, x ) )
    pr y === Ternary(y < 1, 1, Ternary( y > height - y_offset, height - y_offset, y ) )

    temp_type = "#{Int("dummy", :size => temporary_size).type.decl}"

    if not load_overlap then
      total_load_window = total_x_size + 6
      tempload = []
      ranges = []
      temp_vec_length = vector_length
      temp_load_window = total_load_window
      begin
        if temp_vec_length <= temp_load_window then
          tempload.push( Int("tempload#{total_load_window-temp_load_window}_#{total_load_window-temp_load_window+temp_vec_length-1}", :size => 1, :vector_length => temp_vec_length, :signed => false) )
          ranges.push ((total_load_window-temp_load_window)..(total_load_window-temp_load_window+temp_vec_length-1))
          temp_load_window -= temp_vec_length
        else
          temp_vec_length /= 2
        end
      end while temp_load_window > 0
      decl *(tempload)
    else
      tempnn = (0..2).collect { |v_i|
        (0...vector_number).collect { |x_i|
          (0...(y_component_number+2)).collect { |y_i|
            Int("temp#{x_i}#{v_i}#{y_i}", :size => 1, :vector_length => vector_length, :signed => false)
          }
        }
      }
      decl *(tempnn.flatten)
    end
    resnn = (0...(vector_number)).collect { |v_i|
      (0...(y_component_number)).collect { |y_i|
        Int("res#{v_i}#{y_i}", :size => 1, :vector_length => vector_length, :signed => false)
      }
    }
    decl *(resnn.flatten)

    tempcnn = (0..2).collect { |v_i|
      (0...vector_number).collect { |x_i|
        (0...(y_component_number+2)).collect { |y_i|
          Int("tempc#{x_i}#{v_i}#{y_i}", :size => temporary_size, :vector_length => vector_length)
        }
      }
    }
    decl *(tempcnn.flatten)
    rescnn = (0...vector_number).collect { |v_i|
      (0...y_component_number).collect { |y_i|
        Int("resc#{v_i}#{y_i}", :size => temporary_size, :vector_length => vector_length)
      }
    }
    decl *(rescnn.flatten)

    (0...(y_component_number+2)).each { |y_i|
      if not load_overlap then
        load_start = -3
        temp_vec_length = vector_length
        temp_load_window = total_load_window
        vec_indx = 0
        begin
          if temp_vec_length <= temp_load_window then
            pr tempload[vec_indx] === FuncCall("vload#{temp_vec_length}",0, psrc[x + load_start, y + (y_i - 1)].address)
            temp_load_window -= temp_vec_length
            load_start += temp_vec_length
            vec_indx += 1
          else
            temp_vec_length /= 2
          end
        end while temp_load_window > 0
        (0..2).each { |x_i|
          (0...vector_number).each { |v_i|
            start = v_i * vector_length + x_i * 3
            end_indx = start + vector_length - 1
            merge_expr = []
            begin
              vec_indx = ranges.find_index { |r| r.include?(start) }
              if vec_indx then
                vec = tempload[vec_indx]
                r = ranges[vec_indx]
                range = (start-r.begin)..(end_indx < r.end ? end_indx - r.begin : r.end - r.begin)
                #beware valid component idexing length are 1,2,3,4,8,16
                if not [1, 2, 3, 4, 8, 16].include?(range.end - range.begin + 1) then
                  new_size = [1, 2, 3, 4, 8, 16].rindex { |e| e < (range.end - range.begin + 1) }
                  range = range.begin..(range.begin+new_size-1)
                end
                merge_expr.push( "#{vec.components(range)}" )
                start += range.end - range.begin + 1
              else # in the case vectors have dummy elements...
                (end_indx - start + 1).times {
                  merge_expr.push( "0" )
                }
                start = end_indx + 1
              end
            end while start <= end_indx
            pr tempcnn[x_i][v_i][y_i] === FuncCall("convert_#{tempcnn[0][0][0].type.decl}", "(#{tempload.first.type.decl})(#{merge_expr.join(",")})" )
          }
        }
      else
        (0..2).each { |x_i|
          (0...vector_number).each { |v_i|
            pr tempnn[x_i][v_i][y_i] === FuncCall("vload#{vector_length}",0, psrc[x + v_i * vector_length + 3 * (x_i - 1), y + (y_i - 1)].address)
            pr tempcnn[x_i][v_i][y_i] === FuncCall("convert_#{tempcnn[0][0][0].type.decl}", tempnn[x_i][v_i][y_i])
          }
        }
      end
    }
    (0...vector_number).each { |v_i|
      (0...y_component_number).each { |y_i|
        pr rescnn[v_i][y_i] === - tempcnn[0][v_i][y_i]     - tempcnn[1][v_i][y_i]                         - tempcnn[2][v_i][y_i]\
                                - tempcnn[0][v_i][y_i + 1] + tempcnn[1][v_i][y_i + 1] * "(#{temp_type})9" - tempcnn[2][v_i][y_i + 1]\
                                - tempcnn[0][v_i][y_i + 2] - tempcnn[1][v_i][y_i + 2]                     - tempcnn[2][v_i][y_i + 2]
        pr resnn[v_i][y_i] === FuncCall("convert_#{resnn[0][0].type.decl}", clamp(rescnn[v_i][y_i],"(#{temp_type})0","(#{temp_type})255"))
      }
    }

    (0...(y_component_number)).each { |y_i|
      remaining_elem = total_x_size
      (0...vector_number).each { |v_i|
        if remaining_elem >= vector_length then
          pr FuncCall("vstore#{vector_length}", resnn[v_i][y_i], 0, pdst[x + v_i * vector_length, y + y_i].address)
          remaining_elem -= vector_length
        else
          temp_vec_length = vector_length
          begin
            temp_vec_length = temp_vec_length/2
            elem_indexes = 0
            if remaining_elem >= temp_vec_length then
              pr FuncCall("vstore#{temp_vec_length}", resnn[v_i][y_i].components(elem_indexes...(elem_indexes+temp_vec_length)), 0, pdst[x + (v_i * vector_length + elem_indexes), y + y_i].address)
              elem_indexes += temp_vec_length
              remaining_elem -= temp_vec_length
            end
          end while remaining_elem > 0
        end
      }
    }
  }
  pr p
  k.procedure = p
  return k
end
#puts laplacian
#puts laplacian(2)
#puts laplacian(2,2)
#puts laplacian(4,2)
#puts laplacian(3,2)

sizes = [[768, 432], [2560, 1600], [2048, 2048], [5760, 3240], [7680, 4320]]
inputs = []
refs = []
results = []
width = 1024
height = 512

set_lang(C)

k = laplacian_c_ref
puts k
sizes.each { |width, height|
  input = NArray.byte(width*3,height+1).random(256)
  output_ref = NArray.byte(width*3,height)

  k.run(input, output_ref, width, height)
  inputs.push(input)
  refs.push(output_ref[3..-4,1..-2])
  results.push( [] )
}

set_lang(CL)

#(1..16).each { |x_component_number|
[1,2,4,8,16].each { |x_component_number|
  [1,2,4,8,16].reject{ |v| v > x_component_number }.each { |vector_length| # = 1, vector_length=1, y_component_number = 1, temporary_size = 4
    (1..4).each { |y_component_number|
      [2,4].each { |temporary_size|
        [true, false].each { |vector_recompute|
          [true, false].each { |load_overlap|
            id = "x_component_number: #{x_component_number}, vector_length: #{vector_length}, y_component_number: #{y_component_number}, temporary_size: #{temporary_size}, vector_recompute: #{vector_recompute}, load_overlap: #{load_overlap}"
            puts id
            k = laplacian(x_component_number, vector_length, y_component_number, temporary_size, vector_recompute, load_overlap)
            puts k
            sizes.each_index { |indx|
              width, height = sizes[indx]
              puts "#{width} x #{height}"
              output = NArray.byte(width*3,height)
              output.random(256)
              durations=[]
              (0..3).each {
                stats = k.run(inputs[indx], output, width, height, :global_work_size => [rndup((width*3/x_component_number.to_f).ceil,32), (height/y_component_number.to_f).ceil, 1], :local_work_size => [32, 1, 1])
                durations.push stats[:duration]
              }
              puts durations.min
          
#              diff = ( refs[indx] - output[3..-4,1..-2] ).abs
#              i = 0
#              diff.each { |elem|
#                #puts elem
#                i += 1
#                raise "Warning: residue too big: #{elem} #{i%3},#{(i / 3 ) % (width-2)},#{i / 3 / (width - 2)}" if elem > 0
#              }
              results[indx].push( [id, durations.min] )
            }
          }
        }
      }
    }
  }
}
puts "Best candidates:"
results.each_index { |indx|
  width, height = sizes[indx]
  puts "#{width} x #{height}"
  results[indx].sort! { |x,y| x[1] <=> y[1] }
  puts results[indx][0]
}
