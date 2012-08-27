require './Algorithm.rb'
module ConvolutionGenerator
  def ConvolutionGenerator::SimpleLoop(unroll)
    i = Variable::new('i', Int)
    sum = Variable::new('sum', Int)
    data = Variable::new("data",Int,{:dimension => [ Dimension::new(10)]})
    For::new(i, 0, 9, unroll) {
        1.upto(unroll) { |index|
          (sum === sum + data[i + index]).print
        }
    }.print
  end
end

